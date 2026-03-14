#' Fit Directional Quantile Boundary on Residuals Dependent on Covariates
#'
#' After removing the conditional dependence between two responses via a
#' \code{\link{median_pcor}} model, this function fits a directional quantile
#' boundary to the bivariate residuals in polar coordinates. The boundary is
#' estimated as a quantile of the log‑radius \eqn{r = \log(\sqrt{x^2 + y^2})}
#' as a smooth function of the angle \eqn{\alpha = \mathrm{atan2}(y, x)} and the
#' original predictor (e.g., age). The result can be used to construct a
#' prediction region for the original variables.
#'
#' @param x An object of class \code{"median_pcor"} returned by
#'   \code{\link{median_pcor}}.
#' @param qu Numeric; the quantile to estimate for the boundary (e.g., 0.95 for a
#'   95\% region). Passed to \code{\link[qgam]{qgam}}.
#' @param err Numeric; the target error rate for the quantile fit in \code{qgam}.
#'   See \code{\link[qgam]{qgam}} for details.
#' @param ... Additional arguments passed to \code{\link[qgam]{qgam}}.
#'
#' @return An object of S3 class \code{"dirq"} with the following components:
#'   \item{qgam_fit}{The fitted \code{qgam} model for \code{r}.}
#'   \item{qu}{The quantile used.}
#'   \item{predictor}{Name of the predictor variable (the one shared by the two
#'     median models).}
#'   \item{response_vars}{Names of the two original response variables.}
#'   \item{original_residuals}{Data frame containing the standardized residuals
#'     (\code{x} and \code{y}) and the predictor values.}
#'   \item{mods}{List with the three sub‑models from \code{median_pcor}:
#'     \code{y1}, \code{y2}, and \code{cor}.}
#'   \item{scales}{List with the MAD scaling factors for the two responses
#'     (\code{x} and \code{y}).}
#'
#' @importFrom qgam qgam
#' @importFrom stats as.formula
#' @export
#'
#' @examples
#' \donttest{
#' if (requireNamespace("refreg", quietly = TRUE)) {
#'   # Load data and fit a median_pcor model
#'   data_no_dm <- subset(refreg::aegis, dm == "no")
#'   f <- list(fpg ~ s(age), hba1c ~ s(age))
#'   model <- median_pcor(f, data = data_no_dm, qu = 0.5)
#'
#'   # Fit a 95% directional quantile boundary
#'   region <- dirq(model, qu = 0.95, err = 0.01)
#'   print(region)
#'
#'   # Plotting the region on residuals scale
#'   plot(region)
#'   plot(region, newdata = data.frame(age = c(20, 40, 60)))
#'   # Predicting the region for a given covariate value
#'   predict(region, newdata = data.frame(age = 18))
#' } else {
#'   message("Package 'refreg' is needed to run this example.")
#' }
#' }
dirq <- function(x, qu = 0.95, err = 0.05, ...) {

  if (!requireNamespace("qgam", quietly = TRUE)) stop("Package 'qgam' is required.")

  # 1. Extraer nombres
  predictor_var <- all.vars(formula(x$mod.pcor))[2]
  response_vars <- c(names(x$mod.y1$model)[1], names(x$mod.y2$model)[1])

  # 2. Preparar Datos
  df <- x$data
  df$x_res <- x$x.c
  df$y_res <- x$y.c_res

  df$alpha <- atan2(df$y_res, df$x_res)
  df$r     <- log(sqrt(df$x_res^2 + df$y_res^2))

  # 3. Ajuste QGAM (sin especificar k, se usa el valor por defecto de mgcv)
  formula_str <- sprintf("r ~ s(alpha, bs = 'cc') + s(%s)", predictor_var)
  message("Fitting dirq boundary for tau = ", qu, " using predictor: ", predictor_var)
  fit <- qgam::qgam(as.formula(formula_str), data = df, qu = qu, err = err, ...)

  # 4. Construir Objeto S3
  obj <- list(
    qgam_fit = fit,
    qu = qu,
    predictor = predictor_var,
    response_vars = response_vars,
    original_residuals = data.frame(
      x = x$x.c,
      y = x$y.c_res,
      predictor = df[[predictor_var]]
    ),
    mods = list(y1 = x$mod.y1, y2 = x$mod.y2, cor = x$mod.pcor),
    scales = list(
      x = x$scale.x,
      y = x$scale.y
    )
  )

  class(obj) <- "dirq"
  return(obj)
}

#' Print method for \code{dirq} objects
#'
#' Displays basic information about the fitted directional quantile boundary.
#'
#' @param x An object of class \code{"dirq"}.
#' @param ... Additional arguments (ignored).
#' @export
print.dirq <- function(x, ...) {
  # Calculate effective degrees of freedom to show complexity
  edf_total <- sum(x$qgam_fit$edf)
  n_obs <- nrow(x$qgam_fit$model)

  cat("\n-- Directional Quantile Object (dirq) --\n")
  cat(sprintf("%-22s %s\n", "Class:", class(x)))
  cat(sprintf("%-22s %s\n", "Quantile (tau):", x$qu))
  cat(sprintf("%-22s %s\n", "Adjusted Predictor:", x$predictor))
  cat(sprintf("%-22s %d\n", "Observations:", n_obs))
  cat(sprintf("%-22s %s EDF\n", "Model Complexity:", round(edf_total, 2)))
  # Extraer la fórmula del modelo qgam
  form_text <- deparse(formula(x$qgam_fit))
  cat(sprintf("%-22s %s\n", "Formula:", form_text))
  cat("----------------------------------------\n")

  invisible(x)
}
