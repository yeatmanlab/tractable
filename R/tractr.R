#' Analyze group differences in a single dMRI tract profile using GAMs
#'
#' @param df_afq Input AFQ dataframe. If NULL, this function will load data
#'     using read.afq.data and the additional arguments in ...
#' @param dwi_metric The diffusion metric to model (e.g. "FA", "MD")
#' @param participant_id The name of the column that encodes participant ID
#' @param group_by The grouping variable used to group nodeID smoothing terms
#' @param covariates List of strings of GAM covariates,
#'     not including the smoothing terms over nodes and the random effect due
#'     to subjectID.
#' @param smooth_terms Smoothing terms, not including
#'     the smoothing terms over nodes and the random effect due to subjectID.
#' @param comp_list List of factor strings for pairwise comparison.
#'     These should be values from the group column in the AFQ CSV file.
#' @param out_dir Directory in which to save text output and plots
#' @param resampling_technique Type of resampling technique, either "boostrap,"
#'     "permutation," or NULL. If NULL, no resampling test is performed. If
#'     "bootstrap," use bootstrap resampling to estimate confidence intervals.
#'     If "permutation," use permutation resampling to estimate the null
#'     distribution.
#' @param n_samples Number of bootstrap or permutation samples
#' @param k Dimension of the basis used to represent the node smoothing term,
#'     If k = 'auto', this function will attempt to find the best value
#' @param family Distribution to use for the gam. Must be either 'gamma',
#'     'beta', or 'auto'. If 'auto', this function will select the best fit
#'     between beta and gamma distributions.
#' @param ... Arguments to pass to read.afq.files
#'
#' @export
#'
#' @examples
#' \dontrun{
#' sarica <- read.afq.sarica()
#' sarica$group <- factor(sarica$class)
#' sarica$subjectID <- unclass(factor(sarica$subjectID))
#' tractr_bwas(df_afq = sarica,
#'             out_dir = out_dir,
#'             dwi_metric= "fa",
#'             participant_id = "subjectID",
#'             group_by = "group",
#'             covariates = c("age","group"),
#'             comp_list = c("ALS", "CTRL"),
#'             resampling_technique = "bootstrap",
#'             n_samples = 100)
#' }
tractr_bwas <- function(df_afq = NULL,
                        dwi_metric,
                        out_dir,
                        participant_id = "subjectID",
                        group_by = "group",
                        covariates = c(group_by),
                        smooth_terms = NULL,
                        comp_list = unique(df_afq[[group_by]]),
                        resampling_technique = NULL,
                        n_samples = 100,
                        k = "auto",
                        family = "auto",
                        ...) {
  if (is.null(df_afq)) {
    df_afq <- read.afq.files(..., index = participant_id, dwi_metrics = c(dwi_metric))
  }

  tracts <- unique(df_afq$tractID)
  pb <- progress::progress_bar$new(total = length(tracts))

  for (tract in tracts) {
    pb$tick()
    print(tract)
    tractr_single_bundle(df_afq = df_afq,
                         tract = tract,
                         dwi_metric = dwi_metric,
                         out_dir = out_dir,
                         participant_id = participant_id,
                         group_by = group_by,
                         covariates = covariates,
                         smooth_terms = smooth_terms,
                         comp_list = comp_list,
                         resampling_technique = resampling_technique,
                         n_samples = n_samples,
                         k = k,
                         family = family,
                         ...)
  }
}

#' Analyze group differences in a single dMRI tract profile using GAMs
#'
#' @param df_afq Input AFQ dataframe. If NULL, this function will load data
#'     using read.afq.data and the additional arguments in ...
#' @param tract Abbreviated tract name, e.g., "CST_L" or "OR"
#' @param dwi_metric The diffusion metric to model (e.g. "FA", "MD")
#' @param participant_id The name of the column that encodes participant ID
#' @param group_by The grouping variable used to group nodeID smoothing terms
#' @param covariates List of strings of GAM covariates,
#'     not including the smoothing terms over nodes and the random effect due
#'     to subjectID.
#' @param smooth_terms Smoothing terms, not including
#'     the smoothing terms over nodes and the random effect due to subjectID.
#' @param k Dimension of the basis used to represent the node smoothing term,
#'     If k = 'auto', this function will attempt to find the best value
#' @param family Distribution to use for the gam. Must be either 'gamma',
#'     'beta', or 'auto'. If 'auto', this function will select the best fit
#'     between beta and gamma distributions.
#' @param ... Arguments to pass to fit_gam
#'
#' @export
#'
#' @examples
#' sarica <- read.afq.sarica()
#' sarica$group <- factor(sarica$class)
#' sarica$subjectID <- unclass(factor(sarica$subjectID))
#' tractr_single_bundle(df_afq = sarica,
#'                      tract = "Right Corticospinal",
#'                      participant_id = "subjectID",
#'                      group_by = "group",
#'                      covariates = c("age","group"),
#'                      dwi_metric = "fa")
tractr_single_bundle <- function(df_afq,
                                 tract,
                                 dwi_metric,
                                 participant_id = "subjectID",
                                 group_by = "group",
                                 covariates = c(group_by),
                                 smooth_terms = NULL,
                                 k = "auto",
                                 family = "auto",
                                 ...) {
  selected <- select_bundle(
    df_afq = df_afq,
    tract = tract,
    dwi_metric = dwi_metric,
    covariates = covariates,
    participant_id = participant_id,
    group_by = group_by)

  df_tract <- selected$df_tract
  tract_names <- selected$tract_names

  gam_fit <- fit_gam(df_tract = df_tract,
                     target = dwi_metric,
                     covariates = covariates,
                     smooth_terms = smooth_terms,
                     group_by = group_by,
                     participant_id = participant_id,
                     k = k,
                     family = family,
                     ... = ...)

  return(gam_fit)
}
