#' Calculate and plot difference between two splines
#'
#' This function both
#'   1) draws a spline-difference plot for 2 splines,
#'   2) and returns a dataframe of differences at each node
#'
#' @param gam_model GAM object, produced by gam/bam
#' @param tract AFQ tract name
#' @param factor_a First group factor, string
#' @param factor_b Second group factor, string
#' @param out_dir Directory in which to save plots
#'
#' @return A dataframe of spline differences at each node
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
#' df_diff = spline_diff(gam_model = gam_fit,
#'                       tract = "OR",
#'                       factor_a = "0",
#'                       factor_b = "1",
#'                       out_dir = ".")
#' }
spline_diff <- function(gam_model,
                        tract,
                        factor_a,
                        factor_b,
                        out_dir) {
  # determine bottom of plot
  df_pair <- itsadug::plot_diff(
    gam_model,
    view = "nodeID",
    comp = list(group = c(factor_a, factor_b)),
    rm.ranef = T,
    plot = F,
    print.summary = F
  )

  h_min <- min(df_pair$est)
  h_ci <- df_pair[which(df_pair$est == h_min), ]$CI
  min_val <- h_min - h_ci

  # add comparison column to df
  colnames(df_pair) <- c(colnames(df_pair[, 1:4]), "Comp")
  df_pair$Comp <- paste0(factor_a, factor_b)

  # set output
  grDevices::png(
    filename = file.path(out_dir, paste0(
      "plot_diff_", tract, "_pair.png"
    )),
    width = 600, height = 600
  )

  # draw plot
  graphics::par(mar = c(5, 5, 4, 2), family = "Times New Roman")
  p_summary <- utils::capture.output(itsadug::plot_diff(
    gam_model,
    view = "nodeID",
    comp = list(group = c(factor_a, factor_b)),
    rm.ranef = T,
    print.summary = T,
    main = paste0("Difference Scores, ", factor_a, "-", factor_b),
    ylab = "Est. difference",
    xlab = "Tract Node",
    cex.lab = 2,
    cex.axis = 2,
    cex.main = 2,
    cex.sub = 1.5,
    col.diff = "red"
  ))

  # determine sig nodes
  sig_regions <- p_summary[10:length(p_summary)]
  sig_regions <- gsub("\\t", "", sig_regions)

  # make list of start and end nodes, for shading
  sig_list <- as.list(strsplit(sig_regions, " - "))
  start_list <- as.numeric(sapply(sig_list, "[[", 1))
  end_list <- as.numeric(sapply(sig_list, "[[", 2))

  # shade significant regions
  for (h_ind in 1:length(start_list)) {
    graphics::polygon(
      x = c(rep(start_list[h_ind],2), rep(end_list[h_ind], 2)),
      y = c(0, min_val, min_val, 0),
      col = grDevices::rgb(1, 0, 0, 0.2),
      border = NA
    )
  }

  graphics::par(mar = c(5, 4, 4, 2))
  grDevices::dev.off()

  return(df_pair)
}
