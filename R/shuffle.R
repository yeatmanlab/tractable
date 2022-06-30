#' Shuffle an AFQ dataframe
#'
#' This function shuffles participants' age, group, and sex,
#' thereby destroying correlations between the participants'
#' tract profiles and phenotypic data.
#'
#' @param input_df The input AFQ dataframe
#' @param dwi_metric The diffusion MRI metric (e.g. "FA", "MD")
#' @param group.by The grouping variable used to group nodeID smoothing terms
#' @param participant.id The name of the column that encodes participant ID
#' @param shuffle_vars List of strings of column names that should be shuffled
#' @param sample_uniform Boolean flag. If TRUE, shuffling should sample
#'     uniformly from the unique values in the columns. If FALSE, shuffling
#'     will shuffle without replacement.
#'
#' @return A shuffled AFQ dataframe
#' @export
#'
#' @examples
#' \dontrun{
#' df_afq <- read.csv("/path/to/afq/output.csv")
#' df_shuffle <- shuffle_df(df_afq, "dti_fa")
#' }
shuffle_df <- function(input_df, dwi_metric, group.by = "group", participant.id = "subjectID", shuffle_vars = NULL, sample_uniform = FALSE) {
  # Spread the input dataframe to one row per participant
  col_names <- colnames(input_df)
  wide_df <- tidyr::pivot_wider(
    input_df,
    names_from = "nodeID",
    values_from = dwi_metric
  )

  if (is.null(shuffle_vars)) {
    shuffle_vars <- col_names[-which(col_names %in% c("nodeID", "tractID", dwi_metric, participant.id))]
  }

  # Then shuffle participants' shuffle_vars and the grouping variable
  for (svar in unique(c(shuffle_vars, group.by))) {
    if (sample_uniform) {
      # Sample uniformly from the unique values
      wide_df[[svar]] <- sample(unique(wide_df[[svar]]), length(wide_df[[svar]]), replace=TRUE)
    } else {
      # Shuffle the existing values
      wide_df[[svar]] <- sample(wide_df[[svar]], length(wide_df[[svar]]))
    }
  }

  # Gather back to long format (one row per node)
  output_df <- tidyr::pivot_longer(
    wide_df,
    -dplyr::any_of(c(participant.id, shuffle_vars, "tractID", group.by)),
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
#' @param tract AFQ tract name
#' @param group.by The grouping variable used to group nodeID smoothing terms
#' @param participant.id The name of the column that encodes participant ID
#' @param sample_uniform Boolean flag. If TRUE, shuffling should sample
#'     uniformly from the unique values in the columns. If FALSE, shuffling
#'     will shuffle without replacement.
#' @param covariates List of strings of GAM covariates, not including
#'     the smoothing terms over nodes and the random effect due to subjectID.
#' @param smooth_terms Smoothing terms, not including
#'     the smoothing terms over nodes and the random effect due to subjectID.
#' @param k Dimension of the basis used to represent the node smoothing term
#' @param family Distribution to use for the gam. Must be 'gamma' or 'beta'
#' @param formula Optional explicit formula to use for the GAM. If provided,
#'     this will override the dynamically generated formula build from the
#'     target, covariate, and k inputs. Default = NULL.
#' @param factor_a First group factor, string
#' @param factor_b Second group factor, string
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
                             tract,
                             group.by = "group",
                             participant.id = "subjectID",
                             sample_uniform = FALSE,
                             covariates = NULL,
                             smooth_terms = NULL,
                             k = NULL,
                             family = NULL,
                             formula = NULL,
                             factor_a = NULL,
                             factor_b = NULL) {
  if (group.by %in% covariates) {
    coefs <- vector(mode = "list", length = n_permutations)
  }
  node_diffs <- data.frame(nodeID = numeric(0),
                           est = numeric(0),
                           permIdx = numeric(0))

  pb <- progress::progress_bar$new(total = n_permutations)
  for (idx in 1:n_permutations) {
    pb$tick()
    df_shuffle <- shuffle_df(input_df = df_tract,
                             dwi_metric = dwi_metric,
                             group.by = group.by,
                             participant.id = participant.id,
                             shuffle_vars = covariates,
                             sample_uniform = sample_uniform)

    gam_shuffle <- fit_gam(df_tract = df_shuffle,
                           target = dwi_metric,
                           covariates = covariates,
                           smooth_terms = smooth_terms,
                           group.by = group.by,
                           participant.id = participant.id,
                           formula = formula,
                           k = k,
                           family = family)

    if (group.by %in% covariates) {
      coef_name <- grep(paste0("^", group.by),
                        names(gam_shuffle$coefficients),
                        value = TRUE)

      coefs[[idx]] <- gam_shuffle$coefficients[[coef_name]]
    }

    if (!is.null(factor_a) & !is.null(factor_b)) {
      df_pair <- spline_diff(gam_model = gam_shuffle,
                             tract = tract,
                             group.by = group.by,
                             factor_a,
                             factor_b,
                             save_output = FALSE,
                             out_dir = NULL)

      df_pair$permIdx <- idx
      df_pair <- dplyr::select(df_pair, c("nodeID", "est", "permIdx"))
      node_diffs <- rbind(node_diffs, df_pair)
    }
  }

  df_permutation_test <- tidyr::pivot_wider(
    node_diffs,
    names_from = "nodeID",
    names_prefix = "node_",
    values_from = "est",
    id_cols = "permIdx"
  )

  if (group.by %in% covariates) {
    df_permutation_test$group_coefs <- unlist(coefs)
  }

  return(df_permutation_test)
}

