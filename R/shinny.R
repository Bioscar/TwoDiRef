#' Interactive Shiny App for Visualizing DIRQ Reference Regions
#'
#' Launches a Shiny application that displays how the bivariate reference
#' region estimated by \code{\link{dirq}} evolves along a predictor variable.
#' The app shows a sliding window of data points around a chosen predictor
#' value, overlays the corresponding region boundary, and computes the
#' empirical coverage of the region.
#'
#' @param modelo_region An object of class \code{"dirq"} returned by
#'   \code{\link{dirq}}. It contains the fitted directional quantile boundary
#'   and the necessary meta-information (predictor name, response names,
#'   scaling factors, and sub-models).
#' @param data A data frame containing the original variables used to fit the
#'   \code{median_pcor} model (and subsequently the \code{dirq} object). Must
#'   include the predictor and the two response variables with the same names
#'   as stored in \code{modelo_region}.
#' @param window_width Numeric; the half-width of the moving window used to
#'   select data points around the current predictor value. Defaults to 1.
#' @param ... Additional arguments passed to the \code{\link[shiny]{sliderInput}}
#'   for the predictor (e.g., \code{step}, \code{animate}) or to the plot
#'   (currently ignored).
#'
#' @return A Shiny app object. The function is called for its side effect of
#'   launching an interactive web application.
#'
#' @details
#' The app creates a slider for the predictor variable (with range taken from
#' \code{data}). For each selected predictor value \eqn{z_0}, it:
#' \enumerate{
#'   \item Computes the predicted boundary curve using \code{\link{predict.dirq}}.
#'   \item Subsets the data to those observations whose predictor lies within
#'         \eqn{[z_0 - w, z_0 + w]} where \eqn{w} is \code{window_width}.
#'   \item Plots the subset as points, fills the region polygon, and draws its
#'         outline.
#'   \item Calculates the proportion of points inside the polygon (using
#'         \code{\link[sp]{point.in.polygon}}) and displays it as "Empirical Coverage".
#' }
#' The slider can be animated (via the built-in play button) to visualise the
#' evolution smoothly.
#'
#' @importFrom sp point.in.polygon
#' @export
#'
#' @seealso \code{\link{dirq}}, \code{\link{predict.dirq}}, \code{\link{plot.dirq}}
#'
#' @examples
#' \donttest{
#' if (requireNamespace("refreg", quietly = TRUE) && interactive()) {
#'   # Example using the aegis data
#'   data_no_dm <- subset(refreg::aegis, dm == "no")
#'   f <- list(fpg ~ s(age), hba1c ~ s(age))
#'   model <- median_pcor(f, data = data_no_dm, qu = 0.5)
#'   region <- dirq(model, qu = 0.95, err = 0.01)
#'
#'   # Launch the interactive app
#'   shiny_region(region, data_no_dm)
#' }
#' }
shiny_region <- function(modelo_region, data, window_width = 1, ...) {
  # Check that required packages are available
  if (!requireNamespace("shiny", quietly = TRUE)) {
    stop("Package 'shiny' is required to run this function.")
  }
  if (!requireNamespace("sp", quietly = TRUE)) {
    stop("Package 'sp' is required to run this function.")
  }

  # Extract variable names from the dirq object
  pred_name <- modelo_region$predictor
  resp_names <- modelo_region$response_vars  # length 2 vector

  # Check that data contains the required columns
  if (!all(c(pred_name, resp_names) %in% names(data))) {
    stop("'data' must contain columns: ",
         paste(c(pred_name, resp_names), collapse = ", "))
  }

  # Compute range of the predictor for the slider
  pred_range <- range(data[[pred_name]], na.rm = TRUE)

  # Define UI using shiny:: prefix
  ui <- shiny::fluidPage(
    shiny::titlePanel("Temporal Evolution of the Reference Region"),

    shiny::sidebarLayout(
      shiny::sidebarPanel(
        # Slider for the predictor with animation
        shiny::sliderInput("pred_val",
                           label = paste("Select", pred_name, "value:"),
                           min = pred_range[1], max = pred_range[2],
                           value = pred_range[1], step = 1,
                           animate = shiny::animationOptions(interval = 800, loop = TRUE)),

        shiny::hr(),
        shiny::numericInput("window", "Window half-width:",
                            value = window_width, min = 0.1,
                            max = diff(pred_range)/2, step = 0.1),

        shiny::wellPanel(
          shiny::h4("Empirical Coverage"),
          shiny::textOutput("txt_coverage"),
          shiny::tags$small("Points inside the coloured region")
        )
      ),

      shiny::mainPanel(
        shiny::plotOutput("distPlot", height = "600px")
      )
    )
  )

  server <- function(input, output, session) {
    data_calculada <- shiny::reactive({
      newdata <- stats::setNames(data.frame(input$pred_val), pred_name)
      p <- predict(modelo_region, newdata = newdata)
      lower <- input$pred_val - input$window
      upper <- input$pred_val + input$window
      idx <- which(data[[pred_name]] >= lower & data[[pred_name]] <= upper)
      sub_data <- data[idx, , drop = FALSE]
      cov <- if (nrow(sub_data) > 0) {
        mean(sp::point.in.polygon(
          sub_data[[resp_names[1]]], sub_data[[resp_names[2]]],
          p[[resp_names[1]]], p[[resp_names[2]]]
        ) > 0)
      } else {
        NA_real_
      }
      list(p = p, sub_data = sub_data, cov = cov)
    })

    output$distPlot <- shiny::renderPlot({
      res <- data_calculada()
      x_range <- range(data[[resp_names[1]]], na.rm = TRUE)
      y_range <- range(data[[resp_names[2]]], na.rm = TRUE)
      x_pad <- diff(x_range) * 0.05
      y_pad <- diff(y_range) * 0.05
      xlim <- c(x_range[1] - x_pad, x_range[2] + x_pad)
      ylim <- c(y_range[1] - y_pad, y_range[2] + y_pad)
      graphics::plot(res$sub_data[[resp_names[1]]], res$sub_data[[resp_names[2]]],
                     col = "lightgrey", pch = 16, cex = 1.2,
                     xlab = resp_names[1], ylab = resp_names[2],
                     main = paste(pred_name, "=", round(input$pred_val, 2)),
                     xlim = xlim, ylim = ylim)
      graphics::polygon(res$p[[resp_names[1]]], res$p[[resp_names[2]]],
                        col = grDevices::rgb(0, 0, 1, 0.2), border = NA)
      graphics::lines(res$p[[resp_names[1]]], res$p[[resp_names[2]]],
                      col = "blue", lwd = 3)
      graphics::grid(col = "gray90")
    })

    output$txt_coverage <- shiny::renderText({
      res <- data_calculada()
      if (is.na(res$cov)) {
        "No points in window"
      } else {
        paste0(round(res$cov * 100, 1), "%")
      }
    })
  }

  shiny::shinyApp(ui = ui, server = server)
}
