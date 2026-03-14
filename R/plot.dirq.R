# Avoid R CMD check notes for ggplot2 aesthetics
if (getRversion() >= "2.15.1")  utils::globalVariables(c("x", "y", "predictor"))

#' Plot method for dirq objects
#'
#' Produces a ggplot2 visualization of the directional quantile boundaries
#' fitted by \code{\link{dirq}}. The plot shows the standardized residuals
#' (after removing conditional dependence) as points, and overlays the
#' estimated boundaries for specified values of the predictor variable.
#' A reference dashed line through the origin is added for orientation.
#'
#' @param x An object of class \code{"dirq"} returned by \code{\link{dirq}}.
#' @param newdata Optional vector or one‑column data frame containing
#'   new predictor values for which to draw the boundaries. If \code{NULL}
#'   (the default), boundaries are shown at the quartiles of the original
#'   predictor.
#' @param ... Additional arguments (currently ignored).
#'
#' @return The function is called for producing a plot.
#'
#' @importFrom ggplot2 ggplot aes geom_point geom_path geom_hline geom_vline
#' @importFrom ggplot2 scale_color_manual labs theme_minimal theme element_text
#' @importFrom scales seq_gradient_pal
#' @importFrom stats binomial
#' @export
#'
#' @seealso \code{\link{dirq}} for model fitting.
plot.dirq <- function(x, newdata = NULL, ...) {

  pred_name <- x$predictor
  resp_names <- x$response_vars

  # 1. Predictor values for which to draw boundaries
  if (!is.null(newdata)) {
    pred_vals <- if(is.data.frame(newdata)) newdata[[1]] else newdata
  } else {
    pred_vals <- as.numeric(quantile(x$original_residuals$predictor, na.rm = TRUE))
  }

  # 2. Data: standardized residuals
  df_points <- data.frame(
    x = x$original_residuals$x,
    y = x$original_residuals$y
  )

  # 3. Create boundary circles for each predictor value
  boundary_list <- lapply(seq_along(pred_vals), function(i) {
    alpha_seq <- seq(-pi, pi, length.out = 500)
    nd_p <- data.frame(alpha = alpha_seq)
    nd_p[[pred_name]] <- pred_vals[i]
    r_pred <- exp(as.numeric(predict(x$qgam_fit, newdata = nd_p)))
    data.frame(
      x = r_pred * cos(alpha_seq),
      y = r_pred * sin(alpha_seq),
      predictor = pred_vals[i]
    )
  })
  df_boundaries <- do.call(rbind, boundary_list)

  # 4. Colour palette
  boundary_cols <- scales::seq_gradient_pal("#00ced1", "#ff4500", "Lab")(seq(0, 1, length.out = length(pred_vals)))

  # 5. Plot with ggplot2
  gg <- ggplot() +
    geom_point(data = df_points, aes(x = x, y = y), alpha = 0.3, color = "#34495e") +
    geom_path(data = df_boundaries, aes(x = x, y = y, group = predictor, color = as.factor(predictor)), linewidth = 1.5) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray60") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray60") +
    scale_color_manual(
      values = boundary_cols,
      name = pred_name,
      labels = round(pred_vals, 1)
    ) +
    labs(
      x = paste(resp_names[1], "(std. residuals)"),
      y = paste(resp_names[2], "(uncorrelated residuals)"),
      title = paste0("DIRQ Boundaries (tau=", x$qu, ")")
    ) +
    theme_minimal() +
    theme(text = element_text(size = 15), legend.position = "top")

  print(gg)
  invisible(gg)
}
