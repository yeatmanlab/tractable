# Use colorblind-friendly palette from http://jfly.iam.u-tokyo.ac.jp/color/
# The palette with grey:
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# The palette with black:
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

#' Retrieve tract name from abbreviation
#'
#' @param tract AFQ tract abbreviation
#'
#' @return Formatted tract name
#' @export
#'
#' @examples
#' name <- tract_name("OR")
#' name <- tract_name("CST_L")
tract_name <- function(tract) {
  name <- switch(
    tract,
    "OR" = "Optic Radiation",
    "CST_R" = "Right Corticospinal",
    "CST_L" = "Left Corticospinal",
    "UNC_R" = "Right Uncinate",
    "UNC_L" = "Left Uncinate",
    "IFO_L" = "Left IFOF",
    "IFO_R" = "Right IFOF",
    "ARC_R" = "Right Arcuate",
    "ARC_L" = "Left Arcuate",
    "ATR_R" = "Right Thalamic Radiation",
    "ATR_L" = "Left Thalamic Radiation",
    "CGC_R" = "Right Cingulum Cingulate",
    "CGC_L" = "Left Cingulum Cingulate",
    "HCC_R" = "Right Cingulum Hippocampus",
    "HCC_L" = "Left Cingulum Hippocampus",
    "FP" = "Callosum Forceps Major",
    "FA" = "Callosum Forceps Minor",
    "ILF_R" = "RightILF",
    "ILF_L" = "LeftILF",
    "SLF_R" = "RightSLF",
    "SLF_L" = "LeftSLF",
    tract
  )

  return(name)
}


#' Plot GAM splines for each group
#'
#' @param gam_model GAM object, produced by gam/bam
#' @param tract AFQ tract name
#' @param df_tract A dataframe of AFQ nodes for certain tract
#' @param dwi_metric Diffusion MRI metric (e.g. FA, MD)
#' @param covariates List of strings of GAM covariates,
#'     not including the smoothing terms over nodes and the random effect due
#'     to subjectID.
#' @param participant.id The name of the column that encodes participant ID
#' @param group.by The grouping variable used to group nodeID smoothing terms
#' @param out_dir directory in which to save plots
#'
#' @export
#'
#' @examples
#' \dontrun{
#' df_afq <- read.csv("/path/to/afq/output.csv")
#' gam_fit <- fit_gam(df_afq,
#'                    target = "dti_fa",
#'                    covariates = list("group", "sex"),
#'                    family = "gamma",
#'                    k = 40)
#' plot_gam_splines(gam_model = gam_fit,
#'                  tract = "OR",
#'                  df_tract = df_afq,
#'                  dwi_metric = "dti_fa",
#'                  covariates = c("group", "sex"),
#'                  out_dir = ".")
#' }
plot_gam_splines <- function(gam_model, tract, df_tract, dwi_metric, covariates, group.by = "group", participant.id = "subjectID", out_dir) {
  # generate predictions
  values <- vector(mode = "list", length = length(covariates))
  names(values) <- covariates
  df_pred <- mgcv::predict.bam(
    gam_model,
    exclude_terms = c(covariates, "subjectID"),
    values = values,
    se.fit = T,
    type = "response"
  )

  # convert predictions to dataframe
  df_pred <- data.frame(
    nodeID = df_tract$nodeID,
    fit = df_pred$fit,
    se.fit = df_pred$se.fit
  )

  for (covar in covariates) {
    df_pred[[covar]] <- df_tract[[covar]]
  }

  df_pred[[group.by]] <- df_tract[[group.by]]
  df_pred[[participant.id]] <- df_tract[[participant.id]]

  # set up for plot
  h_tract <- tract_name(tract)
  h_title <- paste0("GAM fit of ", h_tract, " ", dwi_metric, " values")

  # draw plot
  options(warn=-1)
  p <- ggplot2::ggplot(data = df_pred) +
    ggplot2::geom_smooth(mapping = ggplot2::aes_string(
      x = "nodeID", y = "fit", color = group.by
    )) +
    ggplot2::ggtitle(h_title) +
    ggplot2::ylab(paste0("Fit ", dwi_metric)) +
    ggplot2::xlab("Tract Node") +
    ggplot2::theme(text = ggplot2::element_text(
      family = "Times New Roman", face = "bold", size = 14
    ))
  options(warn=0)

  # Use colorblind palette for fills and lines
  p + ggplot2::scale_color_manual(
    values = cbbPalette,
  )

  plot_filename <- file.path(out_dir, paste0("plot_gam_", tract, ".png"))

  ggplot2::ggsave(
    plot_filename,
    units = "in",
    width = 6,
    height = 6,
    device = "png"
  )
}
