#' Cross-validation evaluation of coverage for a conditional reference region (dirq)
#'
#' Performs a block cross-validation to evaluate the empirical coverage of a
#' conditional quantile region estimated by \code{dirq}. For each fold, a median regression
#' model is refitted and the corresponding region is recomputed, the coverage
#' indicator (whether the held‑out test points fall inside the region) is recorded.
#' Optionally fits a logistic GAM to explore how coverage varies with the predictor,
#' returns coverage by quartiles, and plots the GAM with confidence bands.
#'
#' @param object An object of class \code{dirq}. It must contain the components:
#'   \code{predictor} (character, name of the predictor variable),
#'   \code{response_vars} (character vector of length two, names of the two
#'   response variables), \code{qu} (numeric, the nominal coverage level, e.g. 0.95),
#'   and \code{mods} (a list of the two fitted models from which formulas can be
#'   extracted). If \code{mods} is not available, the formulas can be supplied
#'   directly via the \code{formulas} argument.
#' @param data Data frame containing the original data used to fit the
#'   \code{dirq} object. Must include the predictor and both response variables.
#' @param fold_size Integer, number of observations per cross‑validation group.
#'   The data are split into groups of (approximately) this size. Default 50.
#' @param err Numeric, tolerance parameter passed to \code{dirq} when
#'   recomputing the region. Default 0.01.
#' @param formulas Optional list of exactly two formulas specifying the models
#'   for the two responses. If \code{NULL} (default), the formulas are extracted
#'   from \code{object$mods}.
#' @param verbose Logical; if \code{TRUE}, print progress messages (in English).
#'   Default \code{FALSE}.
#' @param quiet Logical; if \code{TRUE}, suppress all output from
#'   \code{median_pcor} and \code{dirq} (e.g. the dots showing learning rate
#'   estimation). Default \code{TRUE}.
#' @param show_warnings Logical; if \code{FALSE} (default), warnings from the
#'   final GAM fit are suppressed. Set to \code{TRUE} for debugging.
#' @param plot_gam Logical; if \code{TRUE}, draw a plot of the GAM‑estimated
#'   coverage (with 95\% pointwise confidence intervals) on the percentage scale.
#'   Default \code{FALSE}.
#' @param ci Logical; if \code{TRUE} (default) and \code{plot_gam = TRUE},
#'   include the 95\% confidence band in the plot.
#' @param return_quartiles Logical; if \code{TRUE}, compute coverage by
#'   quartiles of the predictor variable and return the result. Default \code{TRUE}.
#' @param ... Additional arguments passed to \code{median_pcor} and \code{dirq}
#'   (e.g. \code{control} options).
#'
#' @return An invisible list with the following components:
#'   \describe{
#'     \item{overall_coverage}{Numeric, overall empirical coverage (in percent)
#'       rounded to two decimals.}
#'     \item{quartile_coverage}{A data frame with columns \code{quartile} (Q1–Q4),
#'       \code{coverage} (percentage) and \code{n} (number of observations in each
#'       quartile), or \code{NULL} if \code{return_quartiles = FALSE}.}
#'     \item{gam_fit}{A fitted \code{gam} object from \code{mgcv::gam} with
#'       binomial family, modelling \code{.coverage ~ s(predictor)}. \code{NULL}
#'       if the predictor is not found in the data.}
#'     \item{data}{The original data frame augmented with a column
#'       \code{.coverage} (0/1 indicator that the point fell inside its
#'       cross‑validated region).}
#'   }
#'
#' @details
#' The function implements the following steps for each fold:
#' \enumerate{
#'   \item Split the data into training (all but the current group) and test set.
#'   \item Refit a median regression model (quantile 0.5) using \code{median_pcor}
#'         with the two formulas.
#'   \item Recompute the quantile region using \code{dirq} with the original
#'         nominal level \code{object$qu} and tolerance \code{err}.
#'   \item For each observation in the test set, predict the polygon for its
#'         predictor value and test whether the point (response1, response2) lies
#'         inside it using \code{sp::point.in.polygon}.
#' }
#' After all folds, a GAM is fitted to explore how the probability of coverage
#' changes with the predictor. If \code{plot_gam = TRUE}, the estimated coverage
#' curve (with 95\% pointwise confidence intervals) is plotted on the percentage
#' scale. The y‑axis limits are set dynamically: if the nominal level
#' \eqn{\tau < 0.30} the minimum is 0\%, otherwise it is \eqn{\max(0, 100\tau - 15)\%};
#' the maximum is always 100\%.
#'
#' The function requires the packages \pkg{sp} and \pkg{mgcv} to be installed.
#' It also assumes that the functions \code{median_pcor} and \code{dirq} are
#' available in the calling environment (typically from the package that
#' produced the \code{dirq} object, e.g. \pkg{refreg}).
#'
#' @examples
#' \donttest{
#' # Load example data (from the refreg package)
#' library(refreg)
#' data("aegis", package = "refreg")
#' no_dm <- subset(aegis, dm == "no")
#'
#' # Fit median regression and compute 95% region
#' f <- list(fpg ~ s(age), hba1c ~ s(age))
#' modelo <- median_pcor(f = f, data = no_dm, qu = 0.50)
#' region <- dirq(modelo, qu = 0.95, err = 0.01)
#'
#' # Perform cross-validation with groups of 100 observations
#' cv_results <- cv_dirq_coverage(region, no_dm, fold_size = 100,
#'                                verbose = TRUE, quiet = TRUE,
#'                                plot_gam = TRUE, ci = TRUE,
#'                                return_quartiles = TRUE)
#'
#' # Overall coverage (percentage)
#' cv_results$overall_coverage
#'
#' # Coverage by age quartiles
#' cv_results$quartile_coverage
#'
#' # Data with .coverage column
#' head(cv_results$data)
#' }
#'
#' @importFrom sp point.in.polygon
#' @importFrom mgcv gam
#' @importFrom stats quantile aggregate plogis
#' @importFrom utils capture.output
#' @export
cv_dirq_coverage <- function(object, data, fold_size = 50, err = 0.01,
                             formulas = NULL, verbose = FALSE, quiet = TRUE,
                             show_warnings = FALSE, plot_gam = FALSE, ci = TRUE,
                             return_quartiles = TRUE, ...) {
  
  # Check required packages
  if (!requireNamespace("sp", quietly = TRUE)) {
    stop("Package 'sp' is required. Install with install.packages('sp').")
  }
  if (!requireNamespace("mgcv", quietly = TRUE)) {
    stop("Package 'mgcv' is required. Install with install.packages('mgcv').")
  }
  
  # Extract info from dirq object
  if (is.null(object$predictor) || is.null(object$response_vars) || is.null(object$qu)) {
    stop("The 'object' must contain fields: predictor, response_vars, qu")
  }
  predictor_var <- object$predictor
  response_vars <- object$response_vars
  region_qu <- object$qu
  
  # Get formulas
  if (is.null(formulas)) {
    if (is.null(object$mods) || length(object$mods) < 2) {
      stop("No two models found in object$mods. Provide 'formulas' explicitly.")
    }
    f1 <- tryCatch(formula(object$mods[[1]]), error = function(e) NULL)
    f2 <- tryCatch(formula(object$mods[[2]]), error = function(e) NULL)
    if (is.null(f1) || is.null(f2)) {
      stop("Could not extract formulas from object$mods. Provide 'formulas'.")
    }
    formulas <- list(f1, f2)
  } else {
    if (!is.list(formulas) || length(formulas) != 2) {
      stop("'formulas' must be a list of exactly two formulas.")
    }
  }
  
  # Prepare CV groups
  n <- nrow(data)
  indices <- seq_len(n)
  grupos <- split(indices, ceiling(indices / fold_size))
  
  if (verbose) cat("Performing CV with", length(grupos), "groups...\n")
  
  resultados <- integer(n)
  null_file <- if (.Platform$OS.type == "windows") "NUL" else "/dev/null"
  
  for (i in seq_along(grupos)) {
    idx_out <- grupos[[i]]
    train <- data[-idx_out, , drop = FALSE]
    test  <- data[ idx_out, , drop = FALSE]
    
    # Fit median model
    fit_temp <- if (quiet) {
      suppressWarnings(
        suppressMessages(
          capture.output(
            { tmp <- median_pcor(f = formulas, data = train, qu = 0.50, ...) },
            file = null_file
          )
        )
      )
      tmp
    } else {
      median_pcor(f = formulas, data = train, qu = 0.50, ...)
    }
    
    # Compute quantile region
    region_temp <- if (quiet) {
      suppressWarnings(
        suppressMessages(
          capture.output(
            { tmp2 <- dirq(fit_temp, qu = region_qu, err = err, ...) },
            file = null_file
          )
        )
      )
      tmp2
    } else {
      dirq(fit_temp, qu = region_qu, err = err, ...)
    }
    
    # Evaluate each test point
    for (j in seq_len(nrow(test))) {
      punto <- test[j, , drop = FALSE]
      newdata <- data.frame(punto[, predictor_var, drop = FALSE])
      names(newdata) <- predictor_var
      
      poligono <- predict(region_temp, newdata = newdata)
      
      dentro <- sp::point.in.polygon(
        punto[[response_vars[1]]],
        punto[[response_vars[2]]],
        poligono[, 1],
        poligono[, 2]
      )
      resultados[idx_out[j]] <- as.integer(dentro > 0)
    }
    
    if (verbose && i %% 5 == 0) {
      cat("  Processed group", i, "of", length(grupos), "\n")
    }
  }
  
  # Global coverage
  cobertura_total <- mean(resultados)
  
  if (verbose) {
    cat("\n--- CV RESULTS ---\n")
    cat("Observed coverage:", round(cobertura_total * 100, 2), "%\n")
    cat("Expected coverage:", region_qu * 100, "%\n")
  }
  
  # Augment data with coverage column
  data_con_cobertura <- data
  data_con_cobertura$.coverage <- resultados
  
  # Fit GAM (suppress warnings if requested)
  if (!(predictor_var %in% names(data_con_cobertura))) {
    warning("Predictor variable '", predictor_var, "' not found in data. GAM not fitted.")
    gam_fit <- NULL
  } else {
    form_gam <- as.formula(paste(".coverage ~ s(", predictor_var, ")"))
    if (!show_warnings) {
      gam_fit <- suppressWarnings(
        mgcv::gam(form_gam, family = binomial(), data = data_con_cobertura)
      )
    } else {
      gam_fit <- mgcv::gam(form_gam, family = binomial(), data = data_con_cobertura)
    }
  }
  
  # Coverage by quartiles of predictor
  quartile_coverage <- NULL
  if (return_quartiles && predictor_var %in% names(data)) {
    qs <- quantile(data[[predictor_var]], probs = c(0, 0.25, 0.5, 0.75, 1), na.rm = TRUE)
    data_con_cobertura$quartile <- cut(data_con_cobertura[[predictor_var]], 
                                       breaks = qs, 
                                       include.lowest = TRUE,
                                       labels = c("Covariate[Q0-Q25)", "Covariate[Q25-Q50)", 
                                                  "Covariate[Q50-Q75)", "Covariate[Q75-Q100)"))
    q_tab <- aggregate(.coverage ~ quartile, data = data_con_cobertura, FUN = mean)
    q_n   <- aggregate(.coverage ~ quartile, data = data_con_cobertura, FUN = length)
    quartile_coverage <- data.frame(
      quartile = q_tab$quartile,
      coverage = round(100 * q_tab$.coverage, 2),
      n = q_n$.coverage
    )
  }
  
  # Plot GAM if requested
  if (plot_gam && !is.null(gam_fit)) {
    # Sequence of predictor values
    x_seq <- seq(min(data[[predictor_var]], na.rm = TRUE),
                 max(data[[predictor_var]], na.rm = TRUE),
                 length.out = 100)
    newd <- data.frame(x = x_seq)
    names(newd) <- predictor_var
    
    # Predictions with SE on link scale
    pred_link <- predict(gam_fit, newdata = newd, se.fit = TRUE)
    fit_link <- pred_link$fit
    se_link  <- pred_link$se.fit
    
    # Confidence intervals on link scale (95%)
    z <- 1.96
    lower_link <- fit_link - z * se_link
    upper_link <- fit_link + z * se_link
    
    # Transform to response (probability) and then to percentage
    pred_prob <- plogis(fit_link)
    lower_prob <- plogis(lower_link)
    upper_prob <- plogis(upper_link)
    
    pred_pct <- 100 * pred_prob
    lower_pct <- 100 * lower_prob
    upper_pct <- 100 * upper_prob
    
    # Dynamic y-axis limits (in percentage)
    tau <- region_qu
    y_min <- if (tau < 0.30) 0 else max(0, 100 * tau - 15)
    y_max <- 100
    
    # Plot
    plot(x_seq, pred_pct, type = "n",  # start empty
         ylim = c(y_min, y_max),
         xlab = predictor_var, ylab = "Coverage (%)",
         main = "GAM-predicted coverage with 95% CI")
    
    # Add confidence band
    polygon(c(x_seq, rev(x_seq)), 
            c(lower_pct, rev(upper_pct)),
            col = "lightgray", border = NA)
    
    # Add estimated line
    lines(x_seq, pred_pct, col = "blue", lwd = 2)
    
    # Add nominal level
    abline(h = 100 * tau, lty = 2, col = "red")
    
    grid()
    legend("topright", 
           legend = c("Estimated coverage", "95% CI", "Nominal level"),
           lty = c(1, NA, 2), lwd = c(2, NA, 1),
           col = c("blue", NA, "red"),
           fill = c(NA, "lightgray", NA), border = NA,
           bty = "n")
  }
  
  # Return results
  invisible(list(
    overall_coverage = round(100 * cobertura_total, 2),
    quartile_coverage = quartile_coverage,
    gam_fit = gam_fit,
    data = data_con_cobertura
  ))
}