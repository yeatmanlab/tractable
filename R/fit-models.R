#' Build a GAM formula dynamically
#'
#' @param target Diffusion metric to model
#' @param covariates List of strings of GAM covariates, not including
#'     the smoothing terms over nodes and the random effect due to subjectID.
#'     This list can also include smoothing terms.
#' @param smooth_terms Smoothing terms, not including
#'     the smoothing terms over nodes and the random effect due to subjectID.
#' @param group.by The grouping variable used to group nodeID smoothing terms
#' @param participant.id The name of the column that encodes participant ID
#' @param k Dimension of the basis used to represent the node smoothing term
#'
#' @return A GAM formula
#' @export
#'
#' @examples
#' formula <- build_formula(target = "dti_fa",
#'                          covariates = c("group", "sex"),
#'                          k = 40)
#' formula <- build_formula(target = "dki_md",
#'                          covariates = c("group", "sex", "s(age, by=sex)"),
#'                          k = 32)
build_formula <- function(target, covariates, smooth_terms = NULL, group.by = "group", participant.id = "subjectID", k) {
  if (!is.null(covariates)) {
    vars <- paste0(covariates, collapse = "+")
  }
  if (!is.null(group.by)) {
    node_smooth <- paste0("s(nodeID, by = ", group.by, ", k=", k, ")")
  } else {
    node_smooth <- paste0("s(nodeID, k=", k, ")")
  }
  subject_random_effect <- paste0("s(", participant.id, ", bs = \"re\")")
  if (is.null(covariates)) {
    after_tilde <- paste0(list(node_smooth, subject_random_effect), collapse = "+")
  } else {
    after_tilde <- paste0(list(vars, node_smooth, subject_random_effect), collapse = "+")
  }
  if (!is.null(smooth_terms)) {
    after_tilde <- paste0(list(after_tilde, smooth_terms), collapse = "+")
  }
  dyn_string <- paste0(target, " ~ ", after_tilde)
  formula = stats::as.formula(dyn_string)
  return(formula)
}


#' Fit a GAM for tract node metrics (e.g. FA, MD)
#'
#' This function has a series of steps:
#' \itemize{
#'   \item If family == "auto", choose the distribution (either beta or gamma)
#'      that has the lowest AIC when fitting to the dMRI metric data
#'   \item If k == "auto", build an initial GAM model with k = 16 and continue to
#'      double the k value until gam.check shows that k is large enough
#'   \item Fit a GAM model such that: \cr \cr
#'      target ~ covariates + s(nodeID, by=group, k = k_value) + s(subjectID, bs = "re")
#'      \cr
#'   \item Optionally save the output of gam.check and summary to files.
#' }
#'
#' @param df_tract AFQ Dataframe of node metric values for single tract
#' @param target The diffusion metric to model (e.g. "FA", "MD")
#' @param covariates List of strings of GAM covariates, not including
#'     the smoothing terms over nodes and the random effect due to subjectID.
#' @param smooth_terms Smoothing terms, not including
#'     the smoothing terms over nodes and the random effect due to subjectID.
#' @param group.by The grouping variable used to group nodeID smoothing terms
#' @param participant.id The name of the column that encodes participant ID
#' @param formula Optional explicit formula to use for the GAM. If provided,
#'     this will override the dynamically generated formula build from the
#'     target and covariate inputs. Default = NULL.
#' @param k Dimension of the basis used to represent the node smoothing term,
#'     If k = 'auto', this function will attempt to find the best value
#' @param family Distribution to use for the gam. Must be either 'gamma',
#'     'beta', or 'auto'. If 'auto', this function will select the best fit
#'     between beta and gamma distributions.
#' @param tract_name Name of the tract, used only for output file names
#' @param out_dir Directory in which to save gam stats
#' @param save_output Boolean flag to save gam stat files
#'
#' @return Fitted GAM model
#' @export
#'
#' @examples
#' \dontrun{
#' df_afq <- read.csv("/path/to/afq/output.csv")
#' tract <- "CST_L"
#' df_tract <- df_afq[which(df_afq$tractID == tract), ]
#' gam_fit <- fit_gam(df_tract,
#'                    target = "dti_fa",
#'                    covariates = list("group", "sex"),
#'                    family = "gamma",
#'                    k = "auto")
#' }
fit_gam <- function(df_tract,
                    target,
                    covariates = NULL,
                    smooth_terms = NULL,
                    group.by = "group",
                    participant.id = "subjectID",
                    formula = NULL,
                    k = 40,
                    family = "auto",
                    tract_name = "",
                    out_dir = ".",
                    save_output = FALSE) {
  # Set link family
  if (is.character(family) | is.null(family)) {
    if (is.null(family) | tolower(family) == "auto") {
      fit.beta <- fitdistrplus::fitdist(df_tract[[target]], "beta")
      fit.gamma <- fitdistrplus::fitdist(df_tract[[target]], "gamma")
      if (fit.beta$aic < fit.gamma$aic) {
        linkfamily <- mgcv::betar(link = "logit")
      } else {
        linkfamily <- stats::Gamma(link = "logit")
      }
    } else if (tolower(family) == "gamma") {
      linkfamily <- stats::Gamma(link = "logit")
    } else if (tolower(family) == "beta") {
      linkfamily <- mgcv::betar(link = "logit")
    }
  } else {
    linkfamily <- family
  }

  if (is.null(formula)) {
    if (k == "auto") {
      # Initial k value. This will be multiplied by 2 in the while loop
      k.model <- 4

      # Initialize k check results to enter the while loop
      k.indices <- list(0.0, 0.0)
      k.pvals <- list(0.0, 0.0)

      while (any(k.indices < 1.0) | any(k.pvals <= 0.05)) {
        k.model <- k.model * 2
        formula <- build_formula(target = target,
                                 covariates = covariates,
                                 smooth_terms = smooth_terms,
                                 group.by = group.by,
                                 participant.id = participant.id,
                                 k = k.model)

        # Fit the gam
        gam_fit <- mgcv::bam(
          formula,
          data = df_tract,
          family = linkfamily,
          method = "REML"
        )

        k.check <- mgcv::k.check(gam_fit)
        k.indices <- na.omit(k.check[, "k-index"])
        k.pvals <- na.omit(k.check[, "p-value"])
       }
    } else {
      k.model <- k
    }
    formula <- build_formula(target = target,
                             covariates = covariates,
                             smooth_terms = smooth_terms,
                             group.by = group.by,
                             participant.id = participant.id,
                             k = k.model)
  }

  # Fit the gam
  gam_fit <- mgcv::bam(
    formula,
    data = df_tract,
    family = linkfamily,
    method = "REML"
  )

  if (save_output) {
    utils::capture.output(
      mgcv::gam.check(gam_fit, rep = 500),
      file = file.path(out_dir, paste0(
        "k_check_gam_", family, "_", sub(" ", "_", tract_name), ".txt"
      ))
    )
    gam_summary <- summary(gam_fit)
    utils::capture.output(
      gam_summary,
      file = file.path(out_dir, paste0(
        "fit_summary_gam_", family, "_", sub(" ", "_", tract_name), ".txt"
      ))
    )

    utils::write.csv(
      t(gam_summary$p.table[, "Pr(>|t|)"]),
      file.path(out_dir, paste0(
        "gam_pvals_", family, "_", sub(" ", "_", tract_name), ".csv"
        )))

  }

  return(gam_fit)
}
