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
#' plot_gam_splines(gam_fit, "OR", df_afq, "dti_fa", ".")
#' }
plot_gam_splines <- function(gam_model, tract, df_tract, dwi_metric, out_dir) {
  # generate predictions
  df_pred <- mgcv::predict.bam(
    gam_model,
    exclude_terms = c("age", "sex", "subjectID"),
    values = list(age = NULL, sex = NULL),
    se.fit = T,
    type = "response"
  )

  # convert predictions to dataframe
  df_pred <- data.frame(
    Group = df_tract$group,
    sex = df_tract$sex,
    subjectID = df_tract$subjectID,
    age = df_tract$age,
    nodeID = df_tract$nodeID,
    fit = df_pred$fit,
    se.fit = df_pred$se.fit
  )

  # set up for plot
  h_tract <- tract_name(tract)
  h_title <- paste0("GAM fit of ", h_tract, " ", dwi_metric, " values")

  # draw plot
  p <- ggplot2::ggplot(data = df_pred) +
    ggplot2::geom_smooth(mapping = ggplot2::aes_string(
      x = "nodeID", y = "fit", color = "Group"
    )) +
    ggplot2::ggtitle(h_title) +
    ggplot2::ylab(paste0("Fit ", dwi_metric)) +
    ggplot2::xlab("Tract Node") +
    ggplot2::theme(text = ggplot2::element_text(
      family = "Times New Roman", face = "bold", size = 14
    ))

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
