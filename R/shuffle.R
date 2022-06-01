#' Shuffle an AFQ dataframe
#'
#' This function shuffles participants' age, group, and sex,
#' thereby destroying correlations between the participants'
#' tract profiles and phenotypic data.
#'
#' @param input_df The input AFQ dataframe
#' @param dwi_metric The diffusion MRI metric (e.g. "FA", "MD")
#'
#' @return A shuffled AFQ dataframe
#' @export
#'
#' @examples
#' \dontrun{
#' df_afq <- read.csv("/path/to/afq/output.csv")
#' df_shuffle <- shuffle_df(df_afq, "dti_fa")
#' }
shuffle_df <- function(input_df, dwi_metric) {
  # Spread the input dataframe to one row per participant
  col_names <- colnames(input_df)
  wide_df <- tidyr::pivot_wider(
    input_df,
    names_from = "nodeID",
    values_from = dwi_metric
  )

  # Then shuffle participants' age and group
  wide_df$group <- sample(wide_df$group, length(wide_df$group))
  wide_df$age <- sample(wide_df$age, length(wide_df$age))

  # Randomly assign the subjects' sex
  wide_df$sex <- sample(c(0, 1), length(wide_df$sex), replace=TRUE)

  # Gather back to long format (one row per node)
  output_df <- tidyr::pivot_longer(
    wide_df,
    -dplyr::any_of(c("subjectID", "age", "sex", "tractID", "group")),
    names_to = "nodeID",
    values_to = dwi_metric
  )

  output_df <- dplyr::select(output_df, dplyr::all_of(col_names))
  output_df$nodeID <- as.integer(output_df$nodeID)

  return(output_df)
}

#' Simulate the null distribution using permutation testing
#'
#' @param df_tract AFQ Dataframe of node metric values for single tract
#' @param n_permutations Number of permutation tests to perform
#' @param dwi_metric The diffusion metric to model (e.g. "FA", "MD")
#' @param covariates List of strings of GAM covariates, not including
#'     the smoothing terms over nodes and the random effect due to subjectID.
#'     This list can also include smoothing terms.
#' @param k Dimension of the basis used to represent the node smoothing term
#' @param family Distribution to use for the gam. Must be 'gamma' or 'beta'
#' @param formula Optional explicit formula to use for the GAM. If provided,
#'     this will override the dynamically generated formula build from the
#'     target, covariate, and k inputs. Default = NULL.
#'
#' @return Dataframe with permutation test coefficients
#' @export
#'
#' @examples
#' \dontrun{
#' df_afq <- read.csv("/path/to/afq/output.csv")
#' df_tract <- df_afq[which(df_afq$tractID == tract), ]
#' permutation_coefs <- fit_gam(df_afq,
#'                              dwi_metric = "dti_fa",
#'                              covariates = list("group", "sex"),
#'                              family = "gamma",
#'                              k = 40,
#'                              n_permutations = 1000)
#' }
permutation_test <- function(df_tract,
                             n_permutations,
                             dwi_metric,
                             covariates = NULL,
                             k = NULL,
                             family = NULL,
                             formula = NULL) {
  coefs = vector(mode="list", length=n_permutations)

  pb <- progress::progress_bar$new(total = n_permutations)
  for (idx in 1:n_permutations) {
    pb$tick()
    df_shuffle <- shuffle_df(input_df = df_tract, dwi_metric = dwi_metric)
    gam_shuffle <- fit_gam(df_tract = df_shuffle,
                           target = dwi_metric,
                           covariates = covariates,
                           formula = formula,
                           k = k,
                           family = family)

    coefs[[idx]] <- gam_shuffle$coefficients[["group1"]]
  }

  group_coefs <- unlist(coefs)
  df_permutation_test = data.frame(group_coefs)

  return(df_permutation_test)
}

