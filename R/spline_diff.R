#' Calculate and plot difference between two splines
#'
#' This function both
#'   1) draws a spline-difference plot for 2 splines,
#'   2) and returns a dataframe of differences at each node
#'
#' @param gam_model GAM object, produced by gam/bam
#' @param tract AFQ tract name
#' @param group.by The grouping variable used to group nodeID smoothing terms
#' @param factor_a First group factor, string
#' @param factor_b Second group factor, string
#' @param save_output Boolean. If TRUE, save plot output
#' @param sim.ci Logical: Using simultaneous confidence intervals or not
#'   (default set to FALSE). The implementation of simultaneous CIs follows
#'   Gavin Simpson's blog of December 15, 2016:
#'   http://www.fromthebottomoftheheap.net/2016/12/15/simultaneous-interval-revisited/.
#'   This interval is calculated from simulations based. Please specify a seed
#'   (e.g., set.seed(123)) for reproducable results. Note: in contrast with
#'   Gavin Simpson's code, here the Bayesian posterior covariance matrix of
#'   the parameters is uncertainty corrected (unconditional=TRUE) to reflect
#'   the uncertainty on the estimation of smoothness parameters.
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
                        group.by = "group",
                        factor_a,
                        factor_b,
                        save_output = TRUE,
                        sim.ci = FALSE,
                        out_dir) {
  # determine bottom of plot
  comp <- list(c(factor_a, factor_b))
  names(comp) <- c(group.by)
  df_pair <- itsadug::plot_diff(
    gam_model,
    view = "nodeID",
    comp = comp,
    rm.ranef = TRUE,
    plot = FALSE,
    print.summary = FALSE,
    sim.ci = sim.ci,
    n.grid = 100,
  )

  h_min <- min(df_pair$est)
  h_ci <- df_pair[which(df_pair$est == h_min), ]$CI
  min_val <- h_min - h_ci

  # add comparison column to df
  df_pair$comp <- paste0(factor_a, "/", factor_b)

  if (save_output) {
    # set output
    grDevices::png(
      filename = file.path(out_dir, paste0(
        "plot_diff_", sub(" ", "_", tract), "_pair.png"
      )),
      width = 600, height = 600
    )

    # draw plot
    graphics::par(mar = c(5, 5, 4, 2), family = "Times New Roman")

    p_summary <- utils::capture.output(itsadug::plot_diff(
      gam_model,
      view = "nodeID",
      comp = comp,
      # sim.ci = TRUE,
      n.grid = 100,
      rm.ranef = TRUE,
      print.summary = TRUE,
      main = paste0(tract, "Difference Scores, ", factor_a, "-", factor_b),
      ylab = "Est. difference",
      xlab = "Tract Node",
      cex.lab = 2,
      cex.axis = 2,
      cex.main = 2,
      cex.sub = 1.5,
      col.diff = "red"
    ))

    # determine sig nodes
    if (p_summary[14] != "Difference is not significant.") {
      sig_regions <- p_summary[15:length(p_summary)]
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
    }

    graphics::par(mar = c(5, 4, 4, 2))
    grDevices::dev.off()
  }

  return(df_pair)
}
