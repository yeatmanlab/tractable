build_formula <- function(target, covariates, k) {
  # Build a GAM formula dynamically
  #
  # Arguments:
  #   target = diffusion metric to model
  #   covariates = list of strings of GAM covariates, not including the smoothing terms over nodes and the random effect due to subjectID
  #   k = dimension of the basis used to represent the node smoothing term.
  #
  # Returns:
  #   formula = GAM formula
  vars <- paste0(covariates, collapse = "+")
  node_smooth <- paste0("s(nodeID, by = group, k=", k, ")")
  subject_random_effect <- "s(subjectID, bs = \"re\")"
  after_tilde <- paste0(list(vars, node_smooth, subject_random_effect), collapse = "+")
  dyn_string <- paste0(target, " ~ ", after_tilde)
  formula = as.formula(dyn_string)
  return(formula)
}


fit_gam <- function(df_tract, target, covariates, k=40, family="auto",
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
  #   target = diffusion metric to model
  #   covariates = list of strings of GAM covariates, not including the smoothing terms over nodes and the random effect due to subjectID
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
    fit.gamma <- fitdistrplus::fitdist(df_tract[[target]], "gamma")
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
    # Initial k value. This will be multiplied by 2 in the while loop
    k.model <- 8

    # Initialize k check results to enter the while loop
    k.indices <- list(0.0, 0.0)
    k.pvals <- list(0.0, 0.0)

    while (any(k.indices < 1.0) & any(k.pvals <= 0.05)) {
      k.model <- k.model * 2
      formula <- build_formula(target = target, covariates = covariates, k = k.model)

      # Fit the gam
      gam_fit <- mgcv::bam(
        formula,
        data = df_tract,
        family = linkfamily,
        method = "REML"
      )

      k.output <- capture.output(mgcv::gam.check(gam_fit, rep = 500))
      empties <- which(k.output == "")
      table.start <- empties[length(empties)] + 1
      end.sep <- which(k.output == "---")
      table.end <- end.sep[length(end.sep)] - 1
      table.text = k.output[table.start:table.end]
      # Get rid of the significance codes
      table.text <- lapply(table.text, function(line) gsub("[*. ]+$", "", line))
      table.text <- unlist(table.text)
      k.check <- read.table(text = table.text)
      k.indices <- as.numeric(k.check[
        grep("s(nodeID):group", row.names(k.check), fixed = T), "k.index"
      ])
      k.pvals <- as.numeric(k.check[
        grep("s(nodeID):group", row.names(k.check), fixed = T), "p.value"
      ])
    }
  } else {
    k.model <- k
  }

  formula <- build_formula(target = target, covariates = covariates, k = k.model)

  # Fit the gam
  gam_fit <- mgcv::bam(
    formula,
    data = df_tract,
    family = linkfamily,
    method = "REML"
  )

  if (save_output) {
    capture.output(
      mgcv::gam.check(gam_fit, rep = 500),
      file = file.path(out_dir, paste0(
        "k_check_gam_", family, "_", tract_name, ".txt"
      ))
    )

    capture.output(
      mgcv::summary(gam_fit),
      file = file.path(out_dir, paste0(
        "fit_summary_gam_", family, "_", tract, ".txt"
      ))
    )
  }

  return(gam_fit)
}
