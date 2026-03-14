#' Plot method for median_pcor objects
#'
#' Produces two types of plots for objects of class \code{"median_pcor"}:
#' \describe{
#'   \item{\code{"residuals"}}{scatterplot of the standardized residuals
#'         after removing conditional dependence, to check for independence.}
#'   \item{\code{"effects"}}{three-panel display showing the estimated
#'         median trends for each response and the varying partial correlation
#'         coefficient (rho) along the shared covariate.}
#' }
#'
#' @param x An object of class \code{"median_pcor"}.
#' @param type Character string; the type of plot. Can be abbreviated.
#'   Options: \code{"residuals"} or \code{"effects"}.
#' @param cex Numeric multiplier for text size. Passed to \code{par(cex = ...)}
#'   and used in point size for the residuals plot.
#' @param ... Additional graphical parameters.
#'
#' @details
#' For \code{type = "residuals"}, the function plots \code{x$y.c_res} against
#' \code{x$x.c} (the standardized residuals of the second response after
#' removing conditional dependence vs. the standardized residuals of the first
#' response). A reference line at zero is added, and a grid helps assess
#' independence.
#'
#' For \code{type = "effects"}, a layout with three panels is created:
#' \enumerate{
#'   \item{Median trend of the first response (with pointwise 1.96×SE bands).}
#'   \item{Median trend of the second response (with pointwise 1.96×SE bands).}
#'   \item{Estimated varying partial correlation coefficient
#'         \eqn{\rho(y_1, y_2 | \mathbf{z})} as a function of the covariate.}
#' }
#' The confidence bands are based on the standard errors returned by
#' \code{predict(..., se.fit = TRUE)} for the two median models. The correlation
#' panel does not display bands because uncertainty is not trivial to show. Future
#' package version will estimate such estimation error through a bootstrap resampling scheme.
#'
#' A professional colour palette is used, and graphical parameters are set to
#' produce clean, publication‑ready figures.
#'
#' @return The function is called for its side effects (plotting).
#'
#' @export
#' @import graphics
#' @importFrom scales alpha
#' @importFrom stats formula get_all_vars predict
#'
#' @seealso \code{\link{median_pcor}} for model fitting.
plot.median_pcor <- function(x, type = c("residuals", "effects"), cex = 1.1, ...) {
  type <- match.arg(type)

  # Professional color palette
  cols <- list(
    y1     = "#2D708E",
    y2     = "#440154",
    pcor   = "#1F9E89",
    ref    = "#D55E00",
    grid   = "#EEEEEE",
    text   = "#333333",
    accent = "#888888"
  )

  lwd_bold <- 2.5

  # Extract variable names for dynamic labeling
  y1_name <- all.vars(formula(x$mod.y1))[1]
  y2_name <- all.vars(formula(x$mod.y2))[1]

  if (type == "residuals") {
    # --- RESIDUALS PLOT ---
    old_par <- par(no.readonly = TRUE)
    on.exit(par(old_par))

    par(mar = c(5, 5, 4, 2), bg = "white")

    plot(x$x.c, x$y.c_res,
         pch = 21,
         bg = scales::alpha(cols$y1, 0.7), # Solid fill for clarity
         col = "black",                   # Black border for definition
         cex = cex * 1.1,
         xlab = paste(y1_name, "residuals"),
         ylab = paste(y2_name, "residuals"),
         main = "Bivariate Residual Independence",
         axes = FALSE, frame.plot = FALSE, ...)

    # Axes and Grid
    axis(1, col = cols$accent)
    axis(2, col = cols$accent, las = 1)
    grid(nx = NULL, ny = NULL, col = cols$grid, lty = 1)

    # Origin marked with a solid line for reference
    abline(h = 0, v = 0, col = "#444444", lty = 1, lwd = 1.2)

  } else if (type == "effects") {
    # --- EFFECTS PATH PLOT ---
    x_name  <- all.vars(formula(x$mod.pcor))[2]

    old_par <- par(no.readonly = TRUE)
    on.exit(par(old_par))

    # Layout with bottom panel slightly taller
    layout(matrix(c(1, 2, 3, 3), nrow = 2, byrow = TRUE), heights = c(1, 1.2))

    # Global styles: xaxs = "i" ensures the axis covers the full predictor range
    par(mar = c(4.5, 4.5, 3, 2), mgp = c(2.8, 0.7, 0), las = 1, tcl = -0.2, xaxs = "i")

    # --- INTERNAL FUNCTION: MEDIAN TRENDS ---
    draw_median_panel <- function(mod, col_main, title_str) {
      p_data <- mod$model
      all_v  <- all.vars(formula(mod))
      v_name <- all_v[length(all_v)]
      x_vals <- p_data[[v_name]]
      x_lim  <- range(x_vals)

      grid_x  <- seq(x_lim[1], x_lim[2], length.out = 200)
      grid_df <- data.frame(grid_x); names(grid_df) <- v_name

      preds <- predict(mod, newdata = grid_df, se.fit = TRUE)
      y_lim <- range(c(preds$fit + 2*preds$se.fit, preds$fit - 2*preds$se.fit))

      plot(grid_x, preds$fit, type = "n", xlab = v_name, ylab = "Estimate",
           xlim = x_lim, ylim = y_lim, axes = FALSE, frame.plot = FALSE)

      mtext(title_str, side = 3, line = 1, adj = 0, cex = cex, font = 2)
      grid(nx = NULL, ny = NULL, col = cols$grid, lty = 1)
      axis(1, at = pretty(x_lim), col = cols$accent)
      axis(2, col = cols$accent)

      # Confidence Band
      polygon(c(grid_x, rev(grid_x)),
              c(preds$fit + 2*preds$se.fit, rev(preds$fit - 2*preds$se.fit)),
              col = scales::alpha(col_main, 0.15), border = NA)
      lines(grid_x, preds$fit, col = col_main, lwd = lwd_bold)
      rug(x_vals, col = scales::alpha(cols$accent, 0.5))
    }

    # --- INTERNAL FUNCTION: CORRELATION (RHO) ---
    draw_pcor_panel <- function(mod, col_main, title_str, x_label) {
      p_data <- mod$model
      x_vals <- p_data[[x_label]]
      x_lim  <- range(x_vals)

      grid_x  <- seq(x_lim[1], x_lim[2], length.out = 200)
      grid_df <- data.frame(grid_x, x.c = 1); names(grid_df)[1] <- x_label

      fit_vals <- predict(mod, newdata = grid_df)

      # Dynamic zoom on the effect range
      y_range  <- range(fit_vals)
      y_buffer <- diff(y_range) * 0.15
      if(y_buffer == 0) y_buffer <- 0.1
      y_lim    <- c(y_range[1] - y_buffer, y_range[2] + y_buffer)

      plot(grid_x, fit_vals, type = "n", xlab = x_label, ylab = expression(Partial~rho),
           xlim = x_lim, ylim = y_lim, axes = FALSE, frame.plot = FALSE)

      mtext(title_str, side = 3, line = 1, adj = 0, cex = cex, font = 2)
      grid(nx = NULL, ny = NULL, col = cols$grid, lty = 1)
      axis(1, at = pretty(x_lim), col = cols$accent)
      axis(2, col = cols$accent)

      # Reference line at zero only if within visual range
      if(y_lim[1] < 0 && y_lim[2] > 0) {
        abline(h = 0, col = cols$ref, lty = 2, lwd = 1)
      }

      lines(grid_x, fit_vals, col = col_main, lwd = lwd_bold)
      rug(x_vals, col = scales::alpha(cols$accent, 0.5))
    }

    # Draw top panels
    draw_median_panel(x$mod.y1, cols$y1, paste("Median Trend:", y1_name))
    draw_median_panel(x$mod.y2, cols$y2, paste("Median Trend:", y2_name))

    # Draw bottom panel using Rho symbol
    pcor_label <- bquote(bold("Conditional Correlation:")~rho(.(y1_name)*","*.(y2_name))~"|"~.(x_name))
    draw_pcor_panel(x$mod.pcor, cols$pcor, pcor_label, x_name)
  }
}
