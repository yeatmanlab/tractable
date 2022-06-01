#' Analyze group differences in dMRI tract profiles data using GAMs
#'
#' @param afq.csv Input CSV file with AFQ outputs
#' @param tract Abbreviated tract name, e.g., "CST_L" or "OR"
#' @param dwi_metric The diffusion metric to model (e.g. "FA", "MD")
#' @param covariates List of strings of GAM covariates,
#'     not including the smoothing terms over nodes and the random effect due
#'     to subjectID. This list can also include smoothing terms.
#' @param comp_list List of factor strings for pairwise comparison.
#'     These should be values from the group column in the AFQ CSV file.
#' @param out_dir Directory in which to save text output and plots
#' @param permute.test If TRUE, estimate the null distribution of the overall
#'     group coefficient using a permutation test.
#' @param permute.n Number of permutations to simulate the null distribution
#'
#' @export
#'
#' @examples
#' \dontrun{
#' tractr(afq.csv = "/path/to/afq/output.csv",
#'        tract = "CST_R",
#'        dwi_metric = "dti_fa",
#'        covariates = c("sex", "group"),
#'        comp_list = c("0", "1"),
#'        out_dir = ".",
#'        permute.test = TRUE,
#'        permute.n = 100)
#' }
tractr <- function(afq.csv,
                   tract,
                   dwi_metric,
                   covariates,
                   comp_list,
                   out_dir,
                   permute.test = FALSE,
                   permute.n = 100) {
  # Create output directories
  dir.create(out_dir, showWarnings = TRUE, recursive = TRUE)

  plot_dir <- file.path(out_dir, "plots")
  dir.create(plot_dir, showWarnings = TRUE, recursive = TRUE)

  stats_dir <- file.path(out_dir, "stats")
  dir.create(stats_dir, showWarnings = TRUE, recursive = TRUE)

  df_afq = utils::read.csv(afq.csv)

  df_afq$group <- factor(df_afq$group)
  df_afq$subjectID <- unclass(factor(df_afq$subjectID))

  df_afq <- dplyr::select(df_afq, c(
    "subjectID", "nodeID", dwi_metric, "age", "sex", "tractID", "group"
  ))

  df_tract <- df_afq[which(df_afq$tractID == tract), ]

  gam_fit <- fit_gam(df_tract = df_tract,
                     target = dwi_metric,
                     covariates = covariates,
                     k = "auto",
                     family = "auto",
                     tract_name = tract,
                     out_dir = stats_dir,
                     save_output = TRUE)

  df_perm <- permutation_test(df_tract = df_tract,
                              n_permutations = permute.n,
                              dwi_metric = dwi_metric,
                              family = gam_fit$family,
                              formula = gam_fit$formula)

  observed_coef = gam_fit$coefficients[["group1"]]
  group_p_value = sum(
    abs(df_perm$group_coefs) >= observed_coef
  ) / length(df_perm$group_coefs)

  plot_gam_splines(gam_model = gam_fit,
                   tract = tract,
                   df_tract = df_tract,
                   dwi_metric = dwi_metric,
                   out_dir = plot_dir)

  df_diff <- spline_diff(gam_model = gam_fit,
                         tract = tract,
                         factor_a = comp_list[1],
                         factor_b = comp_list[2],
                         out_dir = plot_dir)

  utils::write.csv(df_diff, file.path(stats_dir, "spline_diff.csv"))
}
