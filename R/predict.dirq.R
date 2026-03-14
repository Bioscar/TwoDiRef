#' Predict method for dirq objects
#'
#' Obtains the fitted directional quantile boundary for new predictor values.
#' The function reconstructs the boundary in the original scale of the response
#' variables by combining the quantile model for the log‑radius (fitted in
#' polar coordinates) with the conditional median models and scaling factors
#' stored in the \code{\link{dirq}} object.
#'
#' @param object An object of class \code{"dirq"} returned by \code{\link{dirq}}.
#' @param newdata A vector or one‑column data frame containing the new
#'   predictor value(s) for which the boundary is desired. If multiple values
#'   are supplied, only the first one is used (see Details).
#' @param ... Additional arguments (currently ignored, but kept for
#'   compatibility with the generic \code{\link[stats]{predict}}).
#'
#' @details
#' The prediction proceeds as follows:
#' \enumerate{
#'   \item For a fixed predictor value \eqn{z_0}, a sequence of angles
#'         \eqn{\alpha \in [-\pi, \pi]} is generated.
#'   \item The fitted quantile model for \eqn{r = \log(\sqrt{x^2 + y^2})}
#'         is used to obtain \eqn{\hat{r}(\alpha, z_0)}; the boundary radius
#'         is then \eqn{R(\alpha, z_0) = \exp(\hat{r})}.
#'   \item Standardised residuals are recovered as
#'         \eqn{x_{\text{std}} = R\cos(\alpha)} and
#'         \eqn{y_{\text{std,res}} = R\sin(\alpha)}.
#'   \item The conditional correlation effect is added back:
#'         \eqn{y_{\text{std}} = y_{\text{std,res}} + \hat{\beta}(z_0) x_{\text{std}}},
#'         where \eqn{\hat{\beta}(z_0)} comes from the varying‑coefficient model.
#'   \item Finally, the residuals are back‑transformed to the original scale
#'         using the median predictions and the MAD scaling factors:
#'         \eqn{x_{\text{orig}} = x_{\text{std}} \cdot \text{scale}_x + \hat{\mu}_1(z_0)},
#'         \eqn{y_{\text{orig}} = y_{\text{std}} \cdot \text{scale}_y + \hat{\mu}_2(z_0)}.
#' }
#' If \code{newdata} contains more than one value, only the first is employed,
#' because the function returns a single boundary curve (a closed loop) for a
#' given predictor. To obtain boundaries for multiple predictor values, call
#' \code{predict} repeatedly or use \code{\link{plot.dirq}}.
#'
#' @return A data frame with two columns named after the original response
#'   variables (as stored in \code{object$response_vars}). Each row corresponds
#'   to a point on the predicted boundary curve for the supplied predictor
#'   value.
#'
#' @importFrom stats predict
#' @export
#'
#' @seealso \code{\link{dirq}}, \code{\link{plot.dirq}}
#'
predict.dirq <- function(object, newdata, ...) {
  
  val <- if(is.data.frame(newdata)) newdata[[object$predictor]][1] else newdata[1]
  ga <- seq(-pi, pi, length.out = 250)
  
  # Predictions from the quantile model for the radius
  nd_p <- data.frame(alpha = ga)
  nd_p[[object$predictor]] <- val
  r_val <- exp(as.numeric(predict(object$qgam_fit, newdata = nd_p)))
  
  # Predictions from the conditional median models
  nd_val <- data.frame(temp = val)
  names(nd_val) <- object$predictor
  y1_c <- as.numeric(predict(object$mods$y1, newdata = nd_val))
  y2_c <- as.numeric(predict(object$mods$y2, newdata = nd_val))
  
  # Prediction from the varying‑coefficient (partial correlation) model
  nd_cor <- nd_val
  nd_cor$x.c <- 1
  B_val <- as.numeric(predict(object$mods$cor, newdata = nd_cor))
  
  # Reconstruction to standardized residuals
  x_std <- r_val * cos(ga)
  y_res_std <- r_val * sin(ga)
  y_std <- y_res_std + (B_val * x_std)
  
  # Back‑transform to the original scale using the stored scaling factors
  res_x <- (x_std * object$scales$x) + y1_c
  res_y <- (y_std * object$scales$y) + y2_c
  
  res_df <- data.frame(res_x, res_y)
  colnames(res_df) <- object$response_vars
  return(res_df)
}