fit_gam <- function(df_tract, k=40, family="auto",
                           tract_name="", out_dir=".",
                           save_output=FALSE) {
  # Calculate GAM for tract node metrics (e.g. FA, MD)
  #
  # This function has a series of steps:
  #   1) Foo
  #   2) Bar
  #   3) Baz
  #
  # Arguments:
  #   df_tract = dataframe of node metric values for single tract
  #   k = dimension of the basis used to represent the node smoothing term.
  #       If k='auto', this function will attempt to find the best
  #   family = distribution to use for the gam. Must be either 'gamma', 'beta', or 'auto'
  #            if 'auto', this function will select the best fit between beta and gamma distributions.
  #   tract_name = name of the tract, used only for output file names
  #   out_dir = directory in which to save gam stats
  #   save_output = boolean flag to save gam stat files
  #
  # Writes:
  #   gam_stats_dir/Stats_GAM-*.txt
  #
  # Returns:
  #   fit = GAM object for each tract

  # Set link family
  if (tolower(family) == "gamma") {
    linkfamily <- Gamma(link = "logit")
  } else if (tolower(family) == "beta") {
    linkfamily <- betar(link = "logit")
  } else if (tolower(family) == "auto") {
    fit.beta <- fitdistrplus::fitdist(df_tract[[target]], "beta")
    fit.gamma <- fitdistrplus::fitdist(df_tract$dti_fa, "gamma")
    if (fit.beta$aic < fit.gamma$aic) {
      linkfamily <- mgcv::betar(link = "logit")
    } else {
      linkfamily <- stats::Gamma(link = "logit")
    }
  } else {
    errmessage = paste0(
      "The family argument must be 'gamma', 'beta', or 'auto'. Got ",
      family,
      " instead."
    )
    stop(errmessage)
  }

  if (k == "auto") {
    stop("Automatic setting of k is not yet implemented")
  }

  # Build formula
  vars <- paste0(covariates, collapse = "+")
  node_smooth <- paste0("s(nodeId, by = group, k=", k, ")")
  subject_random_effect <- "s(subjectID, bs = \"re\")"
  after_tilde <- paste0(list(vars, node_smooth, subject_random_effect), collapse = "+")
  dyn_string <- paste0(target, " ~ ", after_tilde)
  formula = as.formula(dyn_string)

  # Fit the gam
  gam_fit <- mgcv::bam(
    formula,
    data = df_tract,
    family = family,
    method = "REML"
  )

  if (save_output) {
    capture.output(
      mgcv::gam.check(gam_fit, rep = 500),
      file = paste0(
        out_dir, "k_check_gam_", family, "_", tract_name, ".txt"
      )
    )

    capture.output(
      mgcv::summary(gam_fit),
      file = paste0(
        out_dir, "fit_summary_gam_", family, "_", tract, ".txt"
      )
    )
  }

  return(gam_fit)
}
