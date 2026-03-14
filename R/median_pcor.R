#' Conditional median and partial correlation model (via QGAM)
#'
#' Fits two quantile GAM (QGAM) models to estimate conditional medians,
#' standardizes the residuals using the Median Absolute Deviation (MAD),
#' and removes smooth conditional dependence between the two responses
#' by means of varying‑coefficient GAM terms.
#'
#' @param f A list of exactly two formulas of the form
#'   `list(y1 ~ s(z1) + s(z2), y2 ~ s(z1) + s(z2))`.
#'   Both formulas must share the same right‑hand side structure.
#' @param data A data frame containing all variables used in the formulas.
#' @param qu Quantile to estimate (default is 0.5 for the conditional median).
#'
#' @return An object of class `"median_pcor"` with components:
#'   \item{x.c}{Standardized residuals of the first response.}
#'   \item{y.c}{Centered residuals of the second response.}
#'   \item{y.c_res}{Residuals of the second response after removing conditional dependence.}
#'   \item{mod.y1}{QGAM model for the first response.}
#'   \item{mod.y2}{QGAM model for the second response.}
#'   \item{mod.pcor}{QGAM varying‑coefficient model for partial correlation.}
#'   \item{scale.x}{MAD of the first‑response residuals (used for scaling).}
#'   \item{scale.y}{MAD of the second‑response residuals.}
#'   \item{data}{Augmented data frame with the computed variables (`x.c`, `y.c`, `y.c_res`).}
#'   \item{qu}{Quantile levels used (default = 0.50).}
#'   \item{call}{Original call.}
#'
#' @importFrom qgam qgam
#' @importFrom stats mad predict
#' @importFrom mgcv s
#' @export
#'
#' @examples
#' \donttest{
#' if (requireNamespace("refreg", quietly = TRUE)) {
#'   # Subset of the aegis data where diabetes is absent
#'   data_no_dm <- subset(refreg::aegis, dm == "no")
#'
#'   # Formulas with the same right‑hand side
#'   f <- list(fpg ~ s(age), hba1c ~ s(age))
#'
#'   # Fit the model
#'   model <- median_pcor(f, data = data_no_dm, qu = 0.5)
#'   print(model)
#'
#'   # Plot methods would be called as:
#'   # plot(model, type = "residuals")
#'   # plot(model, type = "effects")
#' } else {
#'   message("Package 'refreg' is needed to run this example.")
#' }
#' }

median_pcor <- function(f, data, qu = 0.5) {
  
  if (length(f) != 2)
    stop("Argument 'f' must contain exactly two formulas.")
  
  # Check identical RHS structure
  if (!identical(f[[1]][[3]], f[[2]][[3]]))
    stop("Both formulas must share the same right-hand-side structure.")
  
  data_aux <- data
  
  # -------------------------
  # Conditional median models
  # -------------------------
  mod.y1 <- qgam::qgam(f[[1]], data = data_aux, qu = qu)
  mod.y2 <- qgam::qgam(f[[2]], data = data_aux, qu = qu)
  
  y1_name <- all.vars(f[[1]][[2]])
  y2_name <- all.vars(f[[2]][[2]])
  
  # -------------------------
  # Centering
  # -------------------------
  data_aux$x.c <- data_aux[[y1_name]] -
    predict(mod.y1, type = "response")
  
  data_aux$y.c <- data_aux[[y2_name]] -
    predict(mod.y2, type = "response")
  
  # -------------------------
  # MAD scaling
  # -------------------------
  scale.x <- mad(data_aux$x.c)
  scale.y <- mad(data_aux$y.c)
  
  data_aux$x.c <- data_aux$x.c / scale.x
  data_aux$y.c <- data_aux$y.c / scale.y
  
  # -------------------------
  # Build varying-coefficient formula
  # -------------------------
  build_by_formula <- function(formula_obj) {
    
    rhs <- formula_obj[[3]]
    
    # Separate additive terms
    if (inherits(rhs, "call") && rhs[[1]] == as.name("+")) {
      terms_list <- as.list(rhs)[-1]
    } else {
      terms_list <- list(rhs)
    }
    
    by_terms <- lapply(terms_list, function(term) {
      
      # Smooth term
      if (inherits(term, "call") && term[[1]] == as.name("s")) {
        
        args <- as.list(term)[-1]
        args$by <- as.name("x.c")
        
        return(as.call(c(as.name("s"), args)))
        
      } else {
        # Linear term
        return(call(":", as.name("x.c"), term))
      }
    })
    
    rhs_new <- Reduce(function(a, b) call("+", a, b), by_terms)
    
    as.formula(call("~", as.name("y.c"), rhs_new))
  }
  
  pcor_formula <- build_by_formula(f[[1]])
  
  # -------------------------
  # Partial correlation model
  # -------------------------
  mod.pcor <- qgam::qgam(
    pcor_formula,
    data = data_aux,
    qu = qu
  )
  
  # Evaluate beta(z)
  newdata_pcor <- data_aux
  newdata_pcor$x.c <- 1
  
  beta_hat <- predict(
    mod.pcor,
    newdata = newdata_pcor,
    type = "response"
  )
  
  # Remove conditional dependence
  data_aux$y.c_res <- data_aux$y.c -
    beta_hat * data_aux$x.c
  
  # -------------------------
  # Output object
  # -------------------------
  out <- list(
    x.c = data_aux$x.c,
    y.c = data_aux$y.c,
    y.c_res = data_aux$y.c_res,
    mod.y1 = mod.y1,
    mod.y2 = mod.y2,
    mod.pcor = mod.pcor,
    scale.x = scale.x,
    scale.y = scale.y,
    data = data_aux,
    qu = qu,
    call = match.call()
  )
  
  class(out) <- "median_pcor"
  
  return(out)
}

#' Print method for \code{median_pcor} objects
#'
#' @param x An object of class \code{median_pcor}.
#' @param ... Additional arguments (ignored).
#' @export
print.median_pcor <- function(x, ...) {
  
  # Retrieve response variable names
  y1_name <- if(!is.null(x$y1_name)) x$y1_name else all.vars(formula(x$mod.y1))[1]
  y2_name <- if(!is.null(x$y2_name)) x$y2_name else all.vars(formula(x$mod.y2))[1]
  
  # Extract covariates (predictors) from the formula
  covariates <- all.vars(formula(x$mod.y1))[-1]
  cov_text <- paste(covariates, collapse = ", ")
  
  n_obs <- nrow(x$data)
  
  cat("\nMedian-based Partial Correlation Model (via QGAM)")
  cat("\n--------------------------------------------------")
  cat(sprintf("\nQuantile: %s", x$qu))
  cat(sprintf("\nScale, MAD (%s): %.4f", y1_name, x$scale.x))
  cat(sprintf("\nScale, MAD (%s): %.4f", y2_name, x$scale.y))
  cat(sprintf("\nNumber of observations: %d", n_obs))
  cat(sprintf("\nResponse variables: %s and %s", y1_name, y2_name))
  cat(sprintf("\nCovariates used: %s", cov_text))
  
  cat("\n\nModels stored in components:")
  cat(sprintf("\n  $mod.y1   -> Median of %s conditioned on [%s]", y1_name, cov_text))
  cat(sprintf("\n  $mod.y2   -> Median of %s conditioned on [%s]", y2_name, cov_text))
  cat(sprintf("\n  $mod.pcor -> Correlation of %s ~ %s adjusted by [%s]", y2_name, y1_name, cov_text))
  
  cat("\n\nData stored in components ($data):")
  cat(sprintf("\n  $x.c      -> Residuals of %s (centered & scaled)", y1_name))
  cat(sprintf("\n  $y.c      -> Residuals of %s (centered & scaled)", y2_name))
  cat(sprintf("\n  $y.c_res  -> Residuals of %s uncorrelated with %s", y2_name, y1_name))
  cat("\n--------------------------------------------------\n")
  
  invisible(x)
}