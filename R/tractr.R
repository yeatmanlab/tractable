#' Analyze group differences in a single dMRI tract profile using GAMs
#'
#' @param df_afq Input AFQ dataframe. If NULL, this function will load data
#'     using read.afq.data and the additional arguments in ...
#' @param dwi_metric The diffusion metric to model (e.g. "FA", "MD")
#' @param participant.id The name of the column that encodes participant ID
#' @param group.by The grouping variable used to group nodeID smoothing terms
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
#'             participant.id = "subjectID",
#'             group.by = "group",
#'             covariates = c("age","group"),
#'             comp_list = c("ALS", "CTRL"),
#'             resampling_technique = "bootstrap",
#'             n_samples = 100)
#' }
tractr_bwas <- function(df_afq = NULL,
                        dwi_metric,
                        out_dir,
                        participant.id = "subjectID",
                        group.by = "group",
                        covariates = c(group.by),
                        smooth_terms = NULL,
                        comp_list = unique(df_afq[[group.by]]),
                        resampling_technique = NULL,
                        n_samples = 100,
                        k = "auto",
                        family = "auto",
                        ...) {
  if (is.null(df_afq)) {
    df_afq <- read.afq.files(..., index = participant.id, dwi_metrics = c(dwi_metric))
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
                         participant.id = participant.id,
                         group.by = group.by,
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
#' @param participant.id The name of the column that encodes participant ID
#' @param group.by The grouping variable used to group nodeID smoothing terms
#' @param covariates List of strings of GAM covariates,
#'     not including the smoothing terms over nodes and the random effect due
#'     to subjectID.
#' @param smooth_terms Smoothing terms, not including
#'     the smoothing terms over nodes and the random effect due to subjectID.
#' @param comp_list List of factor strings for pairwise comparison.
#'     These should be values from the group column in the AFQ CSV file.
#' @param out_dir Directory in which to save text output and plots
#' @param resampling_technique Type of resampling technique, either "boostrap,"
#'     "permute," or NULL. If NULL, no resampling test is performed. If
#'     "bootstrap," use bootstrap resampling to estimate confidence intervals.
#'     If "permute," use permutation resampling to estimate the null
#'     distribution.
#' @param n_samples Number of bootstrap or permutation samples
#' @param k Dimension of the basis used to represent the node smoothing term,
#'     If k = 'auto', this function will attempt to find the best value
#' @param family Distribution to use for the gam. Must be either 'gamma',
#'     'beta', or 'auto'. If 'auto', this function will select the best fit
#'     between beta and gamma distributions.
#' @param sim.ci Logical: Using simultaneous confidence intervals or not
#'   (default set to FALSE). The implementation of simultaneous CIs follows
#'   Gavin Simpson's blog of December 15, 2016:
#'   http://www.fromthebottomoftheheap.net/2016/12/15/simultaneous-interval-revisited/.
#'   This interval is calculated from simulations based. Please specify a seed
#'   (e.g., set.seed(123)) for reproducable results. Note: in contrast with
#'   Gavin Simpson's code, here the Bayesian posterior covariance matrix of
#'   the parameters is uncertainty corrected (unconditional=TRUE) to reflect
#'   the uncertainty on the estimation of smoothness parameters.
#' @param ... Arguments to pass to read.afq.files
#'
#' @export
#'
#' @examples
#' \dontrun{
#' sarica <- read.afq.sarica()
#' sarica$group <- factor(sarica$class)
#' sarica$subjectID <- unclass(factor(sarica$subjectID))
#' tractr_single_bundle(df_afq = sarica,
#'                      out_dir = ".",
#'                      tract = "Right Corticospinal",
#'                      participant.id = "subjectID",
#'                      group.by = "group",
#'                      covariates = c("age","group"),
#'                      dwi_metric = "fa",
#'                      comp_list = c("ALS", "CTRL"),
#'                      resampling_technique = "bootstrap",
#'                      n_samples = 100)
#' }
tractr_single_bundle <- function(df_afq = NULL,
                                 tract,
                                 dwi_metric,
                                 out_dir,
                                 participant.id = "subjectID",
                                 group.by = "group",
                                 covariates = c(group.by),
                                 smooth_terms = NULL,
                                 comp_list = unique(df_afq[[group.by]]),
                                 resampling_technique = NULL,
                                 n_samples = 100,
                                 k = "auto",
                                 family = "auto",
                                 sim.ci = FALSE,
                                 ...) {
  # Create output directories
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  plot_dir <- file.path(out_dir, "plots")
  dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

  stats_dir <- file.path(out_dir, "stats")
  dir.create(stats_dir, showWarnings = FALSE, recursive = TRUE)

  selected <- select_bundle(
    df_afq,
    tract,
    dwi_metric,
    participant.id,
    covariates,
    group.by)

  df_tract <- selected$df_tract
  tract_names <- selected$tract_names

  gam_fit <- fit_gam(df_tract = df_tract,
                     target = dwi_metric,
                     covariates = covariates,
                     smooth_terms = smooth_terms,
                     group.by = group.by,
                     participant.id = participant.id,
                     k = k,
                     family = family,
                     tract_name = tract,
                     out_dir = stats_dir,
                     save_output = TRUE)

  for (this_tract in tract_names) {
    this_df <- df_tract[which(df_tract$tractID == this_tract), ]

    if (tract == "all") {
      this_comp_list <- grep(paste0("^", this_tract),
                             unique(df_tract[[group.by]]),
                             value = TRUE)
    } else {
      this_comp_list <- comp_list
    }

    if (!is.null(resampling_technique)) {
      permute <- resampling_technique == "permute"

      df_resampling <- sampling_test(df_tract = this_df,
                                     n_samples = n_samples,
                                     dwi_metric = dwi_metric,
                                     tract = this_tract,
                                     group.by = group.by,
                                     participant.id = participant.id,
                                     covariates = covariates,
                                     sample_uniform = TRUE,
                                     family = gam_fit$family,
                                     formula = gam_fit$formula,
                                     factor_a = this_comp_list[1],
                                     factor_b = this_comp_list[2],
                                     permute = permute)

      filename <- paste0(resampling_technique,
                         "_",
                         sub(" ", "_", this_tract),
                         ".csv")
      utils::write.csv(df_resampling,
                       file.path(stats_dir, filename),
                       row.names = FALSE)

      # if (group.by %in% covariates) {
      #   coef_name <- grep(paste0("^", group.by),
      #                     names(gam_fit$coefficients),
      #                     value = TRUE)
      #   observed_coef = gam_fit$coefficients[[coef_name]]
      #   group_p_value = sum(
      #     abs(df_resampling$group_coefs) >= abs(observed_coef)
      #   ) / length(df_resampling$group_coefs)
      #
      #   print(paste0("Bootstrapped group p value = ", group_p_value))
      # }
    }

    plot_gam_splines(gam_model = gam_fit,
                     tract = this_tract,
                     df_tract = this_df,
                     dwi_metric = dwi_metric,
                     covariates = covariates,
                     group.by = group.by,
                     participant.id = participant.id,
                     out_dir = plot_dir)

    if (!is.null(group.by)) {
      df_diff <- spline_diff(gam_model = gam_fit,
                             tract = tract,
                             group.by = group.by,
                             factor_a = this_comp_list[1],
                             factor_b = this_comp_list[2],
                             out_dir = plot_dir,
                             save_output = FALSE,
                             sim.ci = sim.ci)

      filename <- paste0("spline_diff_",
                         sub(" ", "_", this_tract),
                         ".csv")
      utils::write.csv(df_diff,
                       file.path(stats_dir, filename),
                       row.names = FALSE)
    }
  }

  return(gam_fit)
}
