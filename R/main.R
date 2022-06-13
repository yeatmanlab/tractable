#' Analyze group differences in dMRI tract profiles data using GAMs
#'
#' @param df_afq Input AFQ dataframe. If NULL, this function will load data
#'     using load.afq.data and the additional arguments in ...
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
#' @param permute If TRUE, estimate the null distribution of the overall
#'     group coefficient using a permutation test.
#' @param n_permute Number of permutations to simulate the null distribution
#' @param ... Arguments to pass to load.afq.files
#'
#' @export
#'
#' @examples
#' \dontrun{
#' df_afq = read.csv("/path/to/afq/output.csv")
#' tractr(df_afq = df_afq,
#'        tract = "CST_R",
#'        dwi_metric = "dti_fa",
#'        covariates = c("sex", "group"),
#'        comp_list = c("0", "1"),
#'        out_dir = ".",
#'        permute = TRUE,
#'        n_permute = 100)
#' }
tractr <- function(df_afq = NULL,
                   tract,
                   dwi_metric,
                   out_dir,
                   participant.id = "subjectID",
                   group.by = "group",
                   covariates = c(group.by),
                   smooth_terms = NULL,
                   comp_list = unique(df_afq[[group.by]]),
                   permute = FALSE,
                   n_permute = 100,
                   ...) {
  # Create output directories
  dir.create(out_dir, showWarnings = TRUE, recursive = TRUE)

  plot_dir <- file.path(out_dir, "plots")
  dir.create(plot_dir, showWarnings = TRUE, recursive = TRUE)

  stats_dir <- file.path(out_dir, "stats")
  dir.create(stats_dir, showWarnings = TRUE, recursive = TRUE)

  if (is.null(df_afq)) {
    df_afq <- load.afq.files(..., index = participant.id, dwi_metrics = c(dwi_metric))
  }

  df_afq[[group.by]] <- factor(df_afq[[group.by]])
  df_afq[[participant.id]] <- unclass(factor(df_afq[[participant.id]]))

  cols <- unique(c(participant.id,
                   "nodeID",
                   "tractID",
                   dwi_metric,
                   group.by,
                   covariates))
  df_afq <- dplyr::select(df_afq, dplyr::any_of(cols))
  df_tract <- df_afq[which(df_afq$tractID == tract), ]
  df_tract <- stats::na.omit(df_tract)

  gam_fit <- fit_gam(df_tract = df_tract,
                     target = dwi_metric,
                     covariates = covariates,
                     smooth_terms = smooth_terms,
                     group.by = group.by,
                     participant.id = participant.id,
                     k = "auto",
                     family = "auto",
                     tract_name = tract,
                     out_dir = stats_dir,
                     save_output = TRUE)

  if (permute) {
    df_perm <- permutation_test(df_tract = df_tract,
                                n_permutations = n_permute,
                                dwi_metric = dwi_metric,
                                group.by = group.by,
                                participant.id = participant.id,
                                covariates = covariates,
                                sample_uniform = TRUE,
                                family = gam_fit$family,
                                formula = gam_fit$formula)

    coef_name <- grep(paste0("^", group.by),
                      names(gam_fit$coefficients),
                      value = TRUE)
    observed_coef = gam_fit$coefficients[[coef_name]]
    group_p_value = sum(
      abs(df_perm$group_coefs) >= observed_coef
    ) / length(df_perm$group_coefs)

    print(paste0("Bootstrapped group p value = ", group_p_value))
  }

  plot_gam_splines(gam_model = gam_fit,
                   tract = tract,
                   df_tract = df_tract,
                   dwi_metric = dwi_metric,
                   covariates = covariates,
                   group.by = group.by,
                   participant.id = participant.id,
                   out_dir = plot_dir)

  df_diff <- spline_diff(gam_model = gam_fit,
                         tract = tract,
                         group.by = group.by,
                         factor_a = comp_list[1],
                         factor_b = comp_list[2],
                         out_dir = plot_dir)

  utils::write.csv(df_diff, file.path(stats_dir, "spline_diff.csv"))
}
