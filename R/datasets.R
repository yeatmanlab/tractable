#' Parse an S3 URI
#'
#' @param uri An AWS S3 URI
#'
#' @return A named list with bucket and object elements
#' @export
#'
#' @examples
#' parse.s3.uri("s3://bucket-name/a/path/to/an/object.txt")
parse.s3.uri <- function(uri) {
  s3_uri_regex <- stringr::regex("
    ^s3://  # initial S3 prefix
    (?<bucket>(?!(xn--|-s3alias$))[a-z0-9][a-z0-9-]{1,61}[a-z0-9]) # bucket name
    /  # separating /
    (?<object>[a-zA-Z0-9!-_.*'()][a-zA-Z0-9!-_.*'()/]+)$  # object name
    ", comments = TRUE)

  m <- stringr::str_match(uri, s3_uri_regex)[1, ]
  bucket <- getElement(m, "bucket")
  object <- getElement(m, "object")

  return(c("bucket" = bucket, "object" = object))
}

#' Create a merged AFQ/phenotype dataframe
#'
#' @param nodes_csv path to a nodes file
#' @param pheno_csv path to a phenotypic file
#' @param index specification of the column used for merging
#' @param index.nodes specification of the column used for merging
#' @param index.pheno specification of the column used for merging
#' @param dwi_metrics which diffusion metrics should be retained
#' @param factor_cols which columns should be treated as factors
#' @param pheno_cols which columns to include from pheno file
#' @param ... arguments to be passed to read.csv
#'
#' @export
#'
#' @examples
#' \dontrun{
#'   df_afq <- read.afq.files(nodes_csv = "a/path/to/nodes.csv",
#'                            pheno_csv = "a/path/to/pheno.csv")
#' }
read.afq.files <- function(nodes_csv, pheno_csv, index = "subjectID", index.nodes = index, index.pheno = index, dwi_metrics = NULL, factor_cols = NULL, pheno_cols = NULL, ...) {
  if (grepl("^s3://", nodes_csv)) {
    # Read from S3
    uri <- parse.s3.uri(nodes_csv)
    nodes_df <- aws.s3::s3read_using(FUN=utils::read.csv,
                                     bucket = uri[["bucket"]],
                                     object = uri[["object"]],
                                     ... = ...)
  } else {
    # Read from local or URL
    nodes_df <- utils::read.csv(nodes_csv, ... = ...)
  }

  # Select only the user's desired DWI metrics
  if (!is.null(dwi_metrics)) {
    node_cols <- unique(c(index.nodes,
                          "siteID",
                          "tractID",
                          "nodeID",
                          dwi_metrics))
    nodes_df <- dplyr::select(nodes_df, dplyr::any_of(node_cols))
  }

  # Treat some pheno columns as factors
  if (is.null(factor_cols)) {
    colClasses <- NA
  } else {
    colClasses <- rep("factor", length(factor_cols))
    names(colClasses) <- factor_cols
  }

  # Determine separator from file suffix
  if (grepl("\\.tsv$", pheno_csv)) {
    sep <- "\t"
  } else {
    sep <- ","
  }

  if (grepl("^s3://", pheno_csv)) {
    # Read from S3
    uri <- parse.s3.uri(pheno_csv)
    pheno_df <- aws.s3::s3read_using(FUN=utils::read.csv,
                                     bucket = uri[["bucket"]],
                                     object = uri[["object"]],
                                     sep = sep,
                                     ... = ...)
  } else {
    # Read from local or URL
    pheno_df <- utils::read.csv(pheno_csv,
                         check.names = FALSE,
                         colClasses=colClasses,
                         sep = sep,
                         ... = ...)
  }

  # Select only the user supplied pheno columns
  if (!is.null(pheno_cols)) {
    pheno_df <- dplyr::select(pheno_df, unique(c(index.pheno, pheno_cols)))
  }

  # Before merging, add "sub-" prefix to participant IDs if not already there
  pheno_df[[index.pheno]] <- trimws(as.character(pheno_df[[index.pheno]]))
  nodes_df[[index.nodes]] <- trimws(as.character(nodes_df[[index.nodes]]))

  pheno_df[[index.pheno]] <- ifelse(grepl("^sub-", pheno_df[[index.pheno]]),
                                    pheno_df[[index.pheno]],
                                    paste0("sub-", pheno_df[[index.pheno]]))
  nodes_df[[index.nodes]] <- ifelse(grepl("^sub-", nodes_df[[index.nodes]]),
                                    nodes_df[[index.nodes]],
                                    paste0("sub-", nodes_df[[index.nodes]]))

  # Merge and return
  merged_df <- merge(nodes_df, pheno_df, by.x = index.nodes, by.y = index.pheno, incomparables = NA)
  return(merged_df)
}

#' Load tract profiles from Sarica et al.
#'
#' @param ... arguments to be passed to read.csv
#'
#' @return A merged dataframe with data from Sarica et al.
#' @export
#'
#' @examples
#' df_sarica <- read.afq.sarica()
read.afq.sarica <- function(...) {
  url.nodes <- "https://github.com/yeatmanlab/Sarica_2017/raw/gh-pages/data/nodes.csv"
  url.pheno <- "https://github.com/yeatmanlab/Sarica_2017/raw/gh-pages/data/subjects.csv"
  df <- read.afq.files(nodes_csv = url.nodes,
                       pheno_csv = url.pheno,
                       dwi_metrics = c("fa", "md"),
                       factor_cols = c("class", "gender"),
                       pheno_cols = c("age", "class", "gender"),
                       ... = ...)
  return(df)
}

#' Load tract profiles from Yeatman et al.
#'
#' @param ... arguments to be passed to read.csv
#'
#' @return A merged dataframe with data from Yeatman et al.
#' @export
#'
#' @examples
#' df_weston_havens <- read.afq.weston.havens()
read.afq.weston.havens <- function(...) {
  url.nodes <- "https://yeatmanlab.github.io/AFQBrowser-demo/data/nodes.csv"
  url.pheno <- "https://yeatmanlab.github.io/AFQBrowser-demo/data/subjects.csv"
  df <- read.afq.files(nodes_csv = url.nodes,
                       pheno_csv = url.pheno,
                       dwi_metrics = c("fa", "md"),
                       factor_cols = c("Gender"),
                       pheno_cols = c("Age", "Gender", "IQ"),
                       ... = ...)
  return(df)
}

#' Load tract profiles from the Healthy Brain Network dataset
#'
#' @param truncate if TRUE, truncate the data to 49 rows. default = FALSE
#' @param ... arguments to be passed to read.csv
#'
#' @return A merged dataframe with data from HBN
#' @export
#'
#' @examples
#' \dontrun{
#'   df_hbn <- read.afq.hbn()
#' }
read.afq.hbn <- function(truncate = FALSE, ...) {
  if (truncate) {
    url.nodes <- "s3://fcp-indi/data/Projects/HBN/BIDS_curated/derivatives/afq/.truncated_tract_profiles.csv"
    url.pheno <- "s3://fcp-indi/data/Projects/HBN/BIDS_curated/derivatives/afq/.truncated_participants.tsv"
  } else { # nocov start
    url.nodes <- "s3://fcp-indi/data/Projects/HBN/BIDS_curated/derivatives/afq/combined_tract_profiles.csv"
    url.pheno <- "s3://fcp-indi/data/Projects/HBN/BIDS_curated/derivatives/qsiprep/participants.tsv"
  } # nocov end

  df <- read.afq.files(nodes_csv = url.nodes,
                       pheno_csv = url.pheno,
                       index.pheno = "subject_id",
                       dwi_metrics = c("dki_fa", "dki_md"),
                       factor_cols = c("sex"),
                       pheno_cols = c("age", "sex"),
                       ... = ...)
  return(df)
}
