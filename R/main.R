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
#' @param permute If TRUE, estimate the null distribution of the overall
#'     group coefficient using a permutation test.
#' @param n_permute Number of permutations to simulate the null distribution
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
#' df_afq = read.csv("/path/to/afq/output.csv")
#' tractr.bwas(df_afq = df_afq,
#'             tract = "CST_R",
#'             dwi_metric = "dti_fa",
#'             covariates = c("sex", "group"),
#'             comp_list = c("0", "1"),
#'             out_dir = ".",
#'             permute = TRUE,
#'             n_permute = 100)
#' }
tractr_bwas <- function(df_afq = NULL,
                        dwi_metric,
                        out_dir,
                        participant.id = "subjectID",
                        group.by = "group",
                        covariates = c(group.by),
                        smooth_terms = NULL,
                        comp_list = unique(df_afq[[group.by]]),
                        permute = FALSE,
                        n_permute = 100,
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
                         permute = permute,
                         n_permute = n_permute,
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
#' @param permute If TRUE, estimate the null distribution of the overall
#'     group coefficient using a permutation test.
#' @param n_permute Number of permutations to simulate the null distribution
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
#' df_afq = read.csv("/path/to/afq/output.csv")
#' tractr_single_bundle(df_afq = df_afq,
#'                      tract = "CST_R",
#'                      dwi_metric = "dti_fa",
#'                      covariates = c("sex", "group"),
#'                      comp_list = c("0", "1"),
#'                      out_dir = ".",
#'                      permute = TRUE,
#'                      n_permute = 100)
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
                                 permute = FALSE,
                                 n_permute = 100,
                                 k = "auto",
                                 family = "auto",
                                 ...) {
  # Create output directories
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  plot_dir <- file.path(out_dir, "plots")
  dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

  stats_dir <- file.path(out_dir, "stats")
  dir.create(stats_dir, showWarnings = FALSE, recursive = TRUE)

  if (is.null(df_afq)) {
    df_afq <- read.afq.files(..., index = participant.id, dwi_metrics = c(dwi_metric))
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

  if (tract != "all") {
    df_tract <- df_afq[which(df_afq$tractID == tract), ]
    df_tract <- stats::na.omit(df_tract)
    tract_names <- c(tract)
  } else {
    df_tract <- df_afq
    old_group.by <- group.by
    group.by <- paste0("tractID", group.by)
    df_tract[[group.by]] <- as.factor(paste(df_tract$tractID,
                                            df_tract[[old_group.by]]))
    df_tract <- stats::na.omit(df_tract)
    tract_names <- unique(df_tract$tractID)
  }

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

    if (permute) {
      df_perm <- permutation_test(df_tract = this_df,
                                  n_permutations = n_permute,
                                  dwi_metric = dwi_metric,
                                  tract = this_tract,
                                  group.by = group.by,
                                  participant.id = participant.id,
                                  covariates = covariates,
                                  sample_uniform = TRUE,
                                  family = gam_fit$family,
                                  formula = gam_fit$formula,
                                  factor_a = this_comp_list[1],
                                  factor_b = this_comp_list[2])

      filename <- paste0("permutation_test_",
                         sub(" ", "_", this_tract),
                         ".csv")
      utils::write.csv(df_perm,
                       file.path(stats_dir, filename),
                       row.names = FALSE)

      if (group.by %in% covariates) {
        coef_name <- grep(paste0("^", group.by),
                          names(gam_fit$coefficients),
                          value = TRUE)
        observed_coef = gam_fit$coefficients[[coef_name]]
        group_p_value = sum(
          abs(df_perm$group_coefs) >= abs(observed_coef)
        ) / length(df_perm$group_coefs)

        print(paste0("Bootstrapped group p value = ", group_p_value))
      }
    }

    plot_gam_splines(gam_model = gam_fit,
                     tract = this_tract,
                     df_tract = this_df,
                     dwi_metric = dwi_metric,
                     covariates = covariates,
                     group.by = group.by,
                     participant.id = participant.id,
                     out_dir = plot_dir)

    df_diff <- spline_diff(gam_model = gam_fit,
                           tract = tract,
                           group.by = group.by,
                           factor_a = this_comp_list[1],
                           factor_b = this_comp_list[2],
                           out_dir = plot_dir)

    filename <- paste0("spline_diff_",
                       sub(" ", "_", this_tract),
                       ".csv")
    utils::write.csv(df_diff,
                     file.path(stats_dir, filename),
                     row.names = FALSE)
  }
}
