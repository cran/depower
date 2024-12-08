#' Plot power objects
#'
#' An automatic plot method for objects returned by [depower::power()].
#'
#' If you are limited by the output from `plot.depower()`, keep in mind that the
#' object returned by [depower::power()] is a standard data frame. This allows
#' you to easily plot all results with standard plotting functions. In addition,
#' because `plot.depower()` uses ggplot2, you can modify the plot as you
#' normally would. For example:
#'
#' ```{r, plot_details, dev=c('svg','cairo_pdf'), fig.width=5, fig.height=4, fig.show='hide'}
#' set.seed(1234)
#' sim_log_lognormal(
#'   n1 = c(10, 15),
#'   n2 = c(10, 15),
#'   ratio = c(1.3, 1.5),
#'   cv1 = c(0.3),
#'   cv2 = c(0.3, 0.5),
#'   nsims = 1000
#' ) |>
#'   power(alpha = 0.05) |>
#'   plot(hline = 0.8, caption_width = 60) +
#'   ggplot2::theme_bw() +
#'   ggplot2::theme(plot.caption = ggplot2::element_text(hjust = 0)) +
#'   ggplot2::labs(title = "Power for the ratio of geometric means")
#' ```
#'
#' \if{html}{\out{<div style="display: flex; justify-content: center; padding-top: 10px; padding-bottom: 10px;">
#'   <img style="max-width: 100\%; height: auto;" src="figures/plot_details-1.svg" alt="Example of extending plot() with ggplot2 functions" />
#' </div>}}
#' \if{latex}{
#'   \out{\begin{center}}
#'   \figure{plot_details-1.pdf}{options: width=4in}
#'   \out{\end{center}}
#' }
#'
#' @param x (depower)\cr
#'        The data frame returned by [depower::power()].
#' @param x_axis (string: `NULL`; `names(x)`)\cr
#'        The name of the column to be used for the x-axis. Automatically chosen
#'        if `NULL`.
#' @param y_axis (string: `NULL`; `names(x)`)\cr
#'        The name of the column to be used for the y-axis. Automatically chosen
#'        if `NULL`. Generally, `"power"` (default) should be used for the y-axis.
#' @param color (string: `NULL`; `names(x)`)\cr
#'        The name of the column to be used for the [ggplot2::aes()] color
#'        aesthetic. Automatically chosen if `NULL`. Use `NA` to turn off.
#' @param facet_row (string: `NULL`; `names(x)`)\cr
#'        The name of the column to be used for the [ggplot2::facet_grid()] row.
#'        Automatically chosen if `NULL`. Use `NA` to turn off.
#' @param facet_col (string: `NULL`; `names(x)`)\cr
#'        The name of the column to be used for the [ggplot2::facet_grid()]
#'        column. Automatically chosen if `NULL`. Use `NA` to turn off.
#' @param hline (numeric: `NULL`; `(0, 1)`)\cr
#'        The y-intercept at which to draw a horizontal line.
#' @param caption (Scalar logical: `TRUE`)\cr
#'        If `TRUE` (default), a caption is added to the plot. The caption
#'        includes information on parameter values that were conditioned on to
#'        generate the plot. If `FALSE`, the caption is not included.
#' @param caption_width (Scalar integer: `70L`)\cr
#'        The target column number for wrapping the caption text.
#' @param ... Unused additional arguments.
#'
#' @return A [ggplot2::ggplot()] object.
#'
#' @seealso [depower::power()]
#'
#' @examples
#' #----------------------------------------------------------------------------
#' # plot() examples
#' #----------------------------------------------------------------------------
#' library(depower)
#'
#' # Power for independent two-sample t-test
#' set.seed(1234)
#' sim_log_lognormal(
#'   n1 = c(10, 15),
#'   n2 = c(10, 15),
#'   ratio = c(1.3, 1.5),
#'   cv1 = c(0.3),
#'   cv2 = c(0.3, 0.5),
#'   nsims = 500
#' ) |>
#'   power(alpha = 0.05) |>
#'   plot()
#'
#' # Power for dependent two-sample t-test
#' set.seed(1234)
#' sim_log_lognormal(
#'   n1 = c(10, 15),
#'   n2 = c(10, 15),
#'   ratio = c(1.3, 1.5),
#'   cv1 = c(0.3, 0.5),
#'   cv2 = c(0.3, 0.5),
#'   cor = c(0.3),
#'   nsims = 500
#' ) |>
#'   power(alpha = 0.01) |>
#'   plot()
#'
#' # Power for two-sample independent AND two-sample dependent t-test
#' set.seed(1234)
#' sim_log_lognormal(
#'   n1 = c(10, 15),
#'   n2 = c(10, 15),
#'   ratio = c(1.3, 1.5),
#'   cv1 = c(0.3),
#'   cv2 = c(0.3),
#'   cor = c(0, 0.3, 0.6),
#'   nsims = 500
#' ) |>
#'   power(alpha = c(0.05, 0.01)) |>
#'   plot(facet_row = "cor", color = "test")
#'
#' # Power for one-sample t-test
#' set.seed(1234)
#' sim_log_lognormal(
#'   n1 = c(10, 15),
#'   ratio = c(1.2, 1.4),
#'   cv1 = c(0.3, 0.5),
#'   nsims = 500
#' ) |>
#'   power(alpha = c(0.05, 0.01)) |>
#'   plot()
#'
#' \donttest{
#' # Power for independent two-sample NB test
#' set.seed(1234)
#' sim_nb(
#'   n1 = c(10, 15),
#'   mean1 = 10,
#'   ratio = c(1.8, 2),
#'   dispersion1 = 10,
#'   dispersion2 = 3,
#'   nsims = 100
#' ) |>
#'   power(alpha = 0.01) |>
#'   plot()
#'
#' # Power for BNB test
#' set.seed(1234)
#' sim_bnb(
#'   n = c(10, 12),
#'   mean1 = 10,
#'   ratio = c(1.3, 1.5),
#'   dispersion = 5,
#'   nsims = 100
#' ) |>
#'   power(alpha = 0.01) |>
#'   plot()
#' }
#'
#' @importFrom dplyr filter summarize across all_of any_of select left_join
#'                   mutate bind_cols
#' @importFrom ggplot2 ggplot aes facet_grid geom_point geom_line labs vars
#'                     labeller label_both scale_y_continuous scale_x_continuous
#'                     scale_color_discrete geom_hline
#' @importFrom scales percent
#'
#' @export
#' @rdname plot.depower
plot.depower <- function(
    x,
    x_axis = NULL,
    y_axis = NULL,
    color = NULL,
    facet_row = NULL,
    facet_col = NULL,
    hline = NULL,
    caption = TRUE,
    caption_width = 70L,
    ...
) {
  #-----------------------------------------------------------------------------
  # Check arguments
  #-----------------------------------------------------------------------------
  if(nrow(x) < 2L) {
    stop("Argument 'x' must have at least 2 rows.")
  }

  # axis should be null or length==1
  axis <- list(
    x_axis = x_axis, y_axis = y_axis, color = color,
    facet_row = facet_row, facet_col = facet_col
  )
  if(!all(lengths(axis) < 2L)) {
    idx_bad <- which(lengths(axis) > 1L)
    stop(paste0(
      "Argument(s) ",
      csw(names(axis)[idx_bad], and = TRUE),
      " should be a string of length 1 for the column to be used."
    ))
  }
  # axis should be variable in x
  axis <- unlist(axis)
  if(!is.null(axis)) {
    axis <- axis[!is.na(axis)]
    are_bad <- !axis %in% names(x)
    if(any(are_bad)) {
      stop(paste0(
        "Column(s) ",
        csw(axis[are_bad], and = TRUE),
        " are not found in object for argument 'x'."
      ))
    }
  }
  if(!is.null(hline)) {
    if(!is.numeric(hline) || length(hline) != 1L || hline <= 0 || hline >= 1) {
      stop("Argument 'hline' must be a scalar numeric from (0,1).")
    }
  }
  if(!(is.logical(caption) && length(caption) == 1L)) {
    stop("Argument 'caption' must be a scalar logical.")
  }
  if(!(is.numeric(caption_width) && length(caption_width) == 1L && caption_width > 0)) {
    stop("Argument 'caption_width' must be a positive scalar numeric.")
  }

  #-----------------------------------------------------------------------------
  # Jump to next method
  #-----------------------------------------------------------------------------
  NextMethod("plot")
}

#' @export
plot.log_lognormal_independent_two_sample <- function(
    x,
    x_axis = NULL,
    y_axis = NULL,
    color = NULL,
    facet_row = NULL,
    facet_col = NULL,
    hline = NULL,
    caption = TRUE,
    caption_width = 70L,
    ...
) {
  #-----------------------------------------------------------------------------
  # Get axis labels
  #-----------------------------------------------------------------------------
  axis <- axis_builder(
    data = x,
    sorted_axis = c("power", "n2", "n1", "ratio", "cv2", "cv1", "cor", "alpha", "test"),
    x_axis = x_axis,
    y_axis = y_axis,
    color = color,
    facet_row = facet_row,
    facet_col = facet_col
  )

  x_axis <- axis[["x_axis"]]
  y_axis <- axis[["y_axis"]]
  color <- axis[["color"]]
  facet_row <- axis[["facet_row"]]
  facet_col <- axis[["facet_col"]]

  #-----------------------------------------------------------------------------
  # Prepare plot data
  #-----------------------------------------------------------------------------
  plot_data <- data_builder(
    data = x,
    exclude = NULL,
    x_axis = x_axis,
    y_axis = y_axis,
    color = color,
    facet_row = facet_row,
    facet_col = facet_col
  )

  #-----------------------------------------------------------------------------
  # Create plot & return
  #-----------------------------------------------------------------------------
  ggplot(
    data = plot_data$df_plot,
    aes(x = .data[[x_axis]], y = .data[[y_axis]])
  ) +
    geom_builder(
      plot_data$condition_at, x_axis, y_axis, color,
      facet_row, facet_col, hline, caption, caption_width, plot_data$axis_labels
    )
}

#' @export
plot.log_lognormal_dependent_two_sample <- function(
    x,
    x_axis = NULL,
    y_axis = NULL,
    color = NULL,
    facet_row = NULL,
    facet_col = NULL,
    hline = NULL,
    caption = TRUE,
    caption_width = 70L,
    ...
) {
  #-----------------------------------------------------------------------------
  # Get axis labels
  #-----------------------------------------------------------------------------
  axis <- axis_builder(
    data = x,
    sorted_axis = c("power", "n2", "n1", "ratio", "cv2", "cv1", "cor", "alpha", "test"),
    x_axis = x_axis,
    y_axis = y_axis,
    color = color,
    facet_row = facet_row,
    facet_col = facet_col
  )

  x_axis <- axis[["x_axis"]]
  y_axis <- axis[["y_axis"]]
  color <- axis[["color"]]
  facet_row <- axis[["facet_row"]]
  facet_col <- axis[["facet_col"]]

  #-----------------------------------------------------------------------------
  # Prepare plot data
  #-----------------------------------------------------------------------------
  # two sample dependent t-test requires both groups have the same sample size.
  # So n1 or n2 should not be included in conditioning.
  if(any(c(x_axis, y_axis, color, facet_row, facet_col) %in% c("n1", "n2"))) {
    exclude <- c("n1", "n2")
  } else {
    exclude <- NULL
  }

  plot_data <- data_builder(
    data = x,
    exclude = exclude,
    x_axis = x_axis,
    y_axis = y_axis,
    color = color,
    facet_row = facet_row,
    facet_col = facet_col
  )

  # Add sample size descriptor to caption 'n1 = n2'
  # when n1 or n2 is used as axis/color/facet.
  if(any(c(x_axis, y_axis, color, facet_row, facet_col) %in% c("n1", "n2"))) {
    plot_data$condition_at <- c("n1" = "n2", plot_data$condition_at)
  }

  #-----------------------------------------------------------------------------
  # Create plot & return
  #-----------------------------------------------------------------------------
  ggplot(
    data = plot_data$df_plot,
    aes(x = .data[[x_axis]], y = .data[[y_axis]])
  ) +
    geom_builder(
      plot_data$condition_at, x_axis, y_axis, color,
      facet_row, facet_col, hline, caption, caption_width, plot_data$axis_labels
    )
}

#' @export
plot.log_lognormal_mixed_two_sample <- function(
    x,
    x_axis = NULL,
    y_axis = NULL,
    color = NULL,
    facet_row = NULL,
    facet_col = NULL,
    hline = NULL,
    caption = TRUE,
    caption_width = 70L,
    ...
) {
  #-----------------------------------------------------------------------------
  # Get axis labels
  #-----------------------------------------------------------------------------
  axis <- axis_builder(
    data = x,
    sorted_axis = c("power", "n2", "n1", "ratio", "cv2", "cv1", "cor", "alpha", "test"),
    x_axis = x_axis,
    y_axis = y_axis,
    color = color,
    facet_row = facet_row,
    facet_col = facet_col
  )

  x_axis <- axis[["x_axis"]]
  y_axis <- axis[["y_axis"]]
  color <- axis[["color"]]
  facet_row <- axis[["facet_row"]]
  facet_col <- axis[["facet_col"]]

  #-----------------------------------------------------------------------------
  # Prepare plot data
  #-----------------------------------------------------------------------------
  all_cor_0 <- all(x[["cor"]] == 0)
  any_cor_0 <- any(x[["cor"]] == 0)

  # two sample dependent t-test requires both groups have the same sample size.
  # So n1 or n2 should not be included in conditioning.
  if(
    !all_cor_0 &&
    (any(c(x_axis, y_axis, color, facet_row, facet_col) %in% c("n1", "n2")))
  ) {
    exclude <- c("n1", "n2")
  } else {
    exclude <- NULL
  }

  # If we have mixed correlation=0 and correlation>0 data, we need to remove
  # the mixed sample size rows (i.e. n1 != n2) where cor=0.
  if(any_cor_0 && !all_cor_0) {
    x <- x |>
      dplyr::filter(!(.data[["cor"]] == 0 & (.data[["n1"]] != .data[["n2"]])))
  }

  plot_data <- data_builder(
    data = x,
    exclude = exclude,
    x_axis = x_axis,
    y_axis = y_axis,
    color = color,
    facet_row = facet_row,
    facet_col = facet_col
  )

  # Add sample size descriptor to caption 'n1 = n2'
  # when n1 or n2 is used as axis/color/facet.
  if(
    !all_cor_0 &&
    (any(c(x_axis, y_axis, color, facet_row, facet_col) %in% c("n1", "n2")))
  ) {
    plot_data$condition_at <- c("n1" = "n2", plot_data$condition_at)
  }

  #-----------------------------------------------------------------------------
  # Create plot & return
  #-----------------------------------------------------------------------------
  ggplot(
    data = plot_data$df_plot,
    aes(x = .data[[x_axis]], y = .data[[y_axis]])
  ) +
    geom_builder(
      plot_data$condition_at, x_axis, y_axis, color,
      facet_row, facet_col, hline, caption, caption_width, plot_data$axis_labels
    )
}

#' @export
plot.log_lognormal_one_sample <- function(
    x,
    x_axis = NULL,
    y_axis = NULL,
    color = NULL,
    facet_row = NULL,
    facet_col = NULL,
    hline = NULL,
    caption = TRUE,
    caption_width = 70L,
    ...
) {
  #-----------------------------------------------------------------------------
  # Get axis labels
  #-----------------------------------------------------------------------------
  axis <- axis_builder(
    data = x,
    sorted_axis = c("power", "n1", "ratio", "cv1", "alpha", "test"),
    x_axis = x_axis,
    y_axis = y_axis,
    color = color,
    facet_row = facet_row,
    facet_col = facet_col
  )

  x_axis <- axis[["x_axis"]]
  y_axis <- axis[["y_axis"]]
  color <- axis[["color"]]
  facet_row <- axis[["facet_row"]]
  facet_col <- axis[["facet_col"]]

  #-----------------------------------------------------------------------------
  # Prepare plot data
  #-----------------------------------------------------------------------------
  plot_data <- data_builder(
    data = x,
    exclude = NULL,
    x_axis = x_axis,
    y_axis = y_axis,
    color = color,
    facet_row = facet_row,
    facet_col = facet_col
  )

  #-----------------------------------------------------------------------------
  # Create plot & return
  #-----------------------------------------------------------------------------
  ggplot(
    data = plot_data$df_plot,
    aes(x = .data[[x_axis]], y = .data[[y_axis]])
  ) +
    geom_builder(
      plot_data$condition_at, x_axis, y_axis, color,
      facet_row, facet_col, hline, caption, caption_width, plot_data$axis_labels
    )
}

#' @export
plot.nb <- function(
    x,
    x_axis = NULL,
    y_axis = NULL,
    color = NULL,
    facet_row = NULL,
    facet_col = NULL,
    hline = NULL,
    caption = TRUE,
    caption_width = 70L,
    ...
) {
  #-----------------------------------------------------------------------------
  # Get axis labels
  #-----------------------------------------------------------------------------
  axis <- axis_builder(
    data = x,
    sorted_axis = c("power", "n2", "n1", "ratio", "mean1", "dispersion2",
                    "dispersion1", "alpha", "test"),
    x_axis = x_axis,
    y_axis = y_axis,
    color = color,
    facet_row = facet_row,
    facet_col = facet_col
  )

  x_axis <- axis[["x_axis"]]
  y_axis <- axis[["y_axis"]]
  color <- axis[["color"]]
  facet_row <- axis[["facet_row"]]
  facet_col <- axis[["facet_col"]]

  #-----------------------------------------------------------------------------
  # Prepare plot data
  #-----------------------------------------------------------------------------
  plot_data <- data_builder(
    data = x,
    exclude = NULL,
    x_axis = x_axis,
    y_axis = y_axis,
    color = color,
    facet_row = facet_row,
    facet_col = facet_col
  )

  #-----------------------------------------------------------------------------
  # Create plot & return
  #-----------------------------------------------------------------------------
  ggplot(
    data = plot_data$df_plot,
    aes(x = .data[[x_axis]], y = .data[[y_axis]])
  ) +
    geom_builder(
      plot_data$condition_at, x_axis, y_axis, color,
      facet_row, facet_col, hline, caption, caption_width, plot_data$axis_labels
    )
}

#' @export
plot.bnb <- function(
    x,
    x_axis = NULL,
    y_axis = NULL,
    color = NULL,
    facet_row = NULL,
    facet_col = NULL,
    hline = NULL,
    caption = TRUE,
    caption_width = 70L,
    ...
) {
  #-----------------------------------------------------------------------------
  # Get axis labels
  #-----------------------------------------------------------------------------
  axis <- axis_builder(
    data = x,
    sorted_axis = c("power", "n1", "ratio", "mean1", "dispersion1", "alpha", "test"),
    x_axis = x_axis,
    y_axis = y_axis,
    color = color,
    facet_row = facet_row,
    facet_col = facet_col
  )

  x_axis <- axis[["x_axis"]]
  y_axis <- axis[["y_axis"]]
  color <- axis[["color"]]
  facet_row <- axis[["facet_row"]]
  facet_col <- axis[["facet_col"]]

  #-----------------------------------------------------------------------------
  # Prepare plot data
  #-----------------------------------------------------------------------------
  # BNB requires both groups have the same sample size.
  # So n1 or n2 should not be included in conditioning.
  if(any(c(x_axis, y_axis, color, facet_row, facet_col) %in% c("n1", "n2"))) {
    exclude <- c("n1", "n2")
  } else {
    exclude <- NULL
  }

  plot_data <- data_builder(
    data = x,
    exclude = exclude,
    x_axis = x_axis,
    y_axis = y_axis,
    color = color,
    facet_row = facet_row,
    facet_col = facet_col
  )

  # Add sample size descriptor to caption 'n1 = n2'
  # when n1 or n2 is used as axis/color/facet.
  if(any(c(x_axis, y_axis, color, facet_row, facet_col) %in% c("n1", "n2"))) {
    plot_data$condition_at <- c("n1" = "n2", plot_data$condition_at)
  }

  #-----------------------------------------------------------------------------
  # Create plot & return
  #-----------------------------------------------------------------------------
  ggplot(
    data = plot_data$df_plot,
    aes(x = .data[[x_axis]], y = .data[[y_axis]])
  ) +
    geom_builder(
      plot_data$condition_at, x_axis, y_axis, color,
      facet_row, facet_col, hline, caption, caption_width, plot_data$axis_labels
    )
}

#===============================================================================
# Plot helper functions
#===============================================================================
#-------------------------------------------------------------------------------
# Select plot axis
#-------------------------------------------------------------------------------
axis_builder <- function(
    data,
    sorted_axis,
    x_axis = NULL,
    y_axis = NULL,
    color = NULL,
    facet_row = NULL,
    facet_col = NULL
) {
  # Prepare data used for axis assignment
  ## Create a vector of axis names sorted in the order we want to assign them.
  ## Match the sorted vars to the actual data.
  names_data <- names(data)
  idx <- match(sorted_axis, names_data) |>
    remove_na()
  sorted_axis <- names_data[idx]

  ## Number of unique values for each variable of interest.
  n_uniques <- vapply(
    X = data,
    FUN = function(x) {length(unique(x))},
    FUN.VALUE = integer(1),
    USE.NAMES = TRUE
  )
  n_uniques <- n_uniques[sorted_axis]

  ## Only use variables which are not previously used and have 2 or more levels
  used <- c(y_axis, x_axis, color, facet_row, facet_col)
  ok_to_use <- vapply(
    X = sorted_axis,
    FUN = function(x) {!x %in% used},
    FUN.VALUE = logical(1),
    USE.NAMES = TRUE
  )

  axis <- n_uniques > 1L & ok_to_use
  axis <- axis[axis] |>
    remove_na() |>
    names()

  ## In a dependent two-sample t-Test, n1 shouldn't be used when n2 is already
  ## in use, because they must be the same value.
  if("cor" %in% names_data) {
    ok_to_use_n1 <- all(data[["cor"]] == 0)
    if(!ok_to_use_n1) {
      axis <- axis[axis != "n1"]
    }
  }
  ## In a BNB test, n2 shouldn't be used when n1 is already in use, because they
  ## must be the same value.
  if(inherits(data, "bnb")) {
    axis <- axis[axis != "n2"]
  }

  # Assign axis where needed
  if(is.null(y_axis)) {
    # handle pathological cases
    if(isTRUE(ok_to_use["power"]) && !"power" %in% axis) {
      y_axis <- "power"
    } else {
      y_axis <- axis[1L]
      axis <- remove_used_axis(axis)
    }
  }

  if(is.null(x_axis)) {
    # handle pathological cases
    if(isTRUE(ok_to_use["n2"]) && is.na(axis[1L])) {
      x_axis <- "n2"
    } else if(isTRUE(ok_to_use["n1"]) && is.na(axis[1L])) {
      x_axis <- "n1"
    } else if(isTRUE(ok_to_use["n"]) && is.na(axis[1L])) {
      x_axis <- "n"
    } else {
      x_axis <- axis[1L]
      axis <- remove_used_axis(axis)
    }
  }

  if(is.null(color)) {
    color <- axis[1L]
    axis <- remove_used_axis(axis)
  }

  if(is.null(facet_row)) {
    facet_row <- axis[1L]
    axis <- remove_used_axis(axis)
  }

  if(is.null(facet_col)) {
    facet_col <- axis[1L]
    axis <- remove_used_axis(axis)
  }

  # Return
  list(
    y_axis = y_axis,
    x_axis = x_axis,
    color = color,
    facet_row = facet_row,
    facet_col = facet_col
  )
}

#-------------------------------------------------------------------------------
# Create list of data items for creating the plot
#-------------------------------------------------------------------------------
data_builder <- function(
    data,
    exclude = NULL,
    x_axis,
    y_axis,
    color,
    facet_row,
    facet_col
) {
  # Combine factors and get labels
  axis <- c(x_axis, y_axis, color, facet_row, facet_col)
  axis <- axis[!is.na(axis)]

  axis_labels <- vapply(
    X = data[axis],
    FUN = function(x) {attr(x, "label", exact = TRUE)},
    FUN.VALUE = character(1L),
    USE.NAMES = TRUE
  )

  #-----------------------------------------------------------------------------
  # Determine what we are conditioning on and get values to hold constant for
  # the left_join below. We don't want to plot data that is being conditioned on.
  # Because nsims may unexpectedly be reduced, we will need to handle this separately.
  #-----------------------------------------------------------------------------
  # Exclude unnecessary variables.
  unused <- setdiff(names(data), axis)
  unused <- unused[!unused %in% c("nsims", "data", "result")]

  # Don't lose labels during mutate (middle.labelled)
  condition_at <- data |>
    summarize(across(.cols = all_of(unused), .fns = middle.labelled)) |>
    select(
      any_of(
        c("n1", "n2", "n", "ratio", "mean1", "cv", "cv1", "cv2", "dispersion1",
          "dispersion2", "dispersion", "cor", "alpha", "alternative", "mean_null",
          "ratio_null", "test")
      )
    ) |>
    select(
      -any_of({{exclude}})
    )

  # Don't want these variables to be numeric, so factorize them
  vars_to_factorize <- c(color)
  vars_to_factorize <- vars_to_factorize[!is.na(vars_to_factorize)]
  if(length(vars_to_factorize) == 0L) {vars_to_factorize <- NULL}

  # Create data frame for plotting
  df_plot <- condition_at |>
    left_join(data, by = names(condition_at)) |>
    mutate(across(.cols = any_of({{vars_to_factorize}}), .fns = as.factor))

  # now add nsims back to condition_at
  if("nsims" %in% names(data)) {
    df_nsims <- data |>
      select(nsims) |>
      filter(!duplicated(nsims))
    if(nrow(df_nsims) > 1L) {
      df_nsims <- df_nsims |>
        summarize(nsims = paste0("Varying: ", "min=", min(nsims), ",max=", max(nsims)))
      attr(df_nsims[["nsims"]], "label") <- "N Simulations"
    }
    condition_at <- bind_cols(condition_at, df_nsims)
  }

  # Replace names with labels and convert to list
  names(condition_at) <- vapply(
    X = condition_at,
    FUN = function(x) {attr(x, "label", exact = TRUE)},
    FUN.VALUE = character(1L),
    USE.NAMES = FALSE
  )
  condition_at <- as.list(condition_at)

  # Return
  list(df_plot = df_plot, axis_labels = axis_labels, condition_at = condition_at)
}

#-------------------------------------------------------------------------------
# Build geoms using a list
#-------------------------------------------------------------------------------
geom_builder <- function(
    condition_at,
    x_axis,
    y_axis,
    color,
    facet_row,
    facet_col,
    hline,
    caption,
    caption_width,
    axis_labels
) {
  list(
    if(is.na(color)) {
      geom_point()
    } else {
      geom_point(aes(color = .data[[color]], group = .data[[color]]))
    },
    if(is.na(color)) {
      geom_line()
    } else {
      geom_line(aes(color = .data[[color]], group = .data[[color]]))
    },
    if(x_axis %in% c("n1", "n2", "n")) {
      scale_x_continuous(breaks = ~round(unique(pretty(.)))) # integer breaks
    },
    if(x_axis %in% c("power")) {
      scale_x_continuous(limits = c(0, 1), labels = percent)
    },
    if(y_axis %in% c("n1", "n2", "n")) {
      scale_y_continuous(breaks = ~round(unique(pretty(.)))) # integer breaks
    },
    if(y_axis %in% c("power")) {
      scale_y_continuous(limits = c(0, 1), labels = percent)
    },
    if(!is.na(facet_row) && is.na(facet_col)) {
      facet_grid(
        rows = vars(.data[[facet_row]]),
        labeller = labeller(.rows = label_both2)
      )
    },
    if(is.na(facet_row) && !is.na(facet_col)) {
      facet_grid(
        cols = vars(.data[[facet_col]]),
        labeller = labeller(.cols = label_both2)
      )
    },
    if(!is.na(facet_row) && !is.na(facet_col)) {
      facet_grid(
        rows = vars(.data[[facet_row]]),
        cols = vars(.data[[facet_col]]),
        labeller = labeller(.rows = label_both2, .cols = label_both2)
      )
    },
    scale_color_discrete(labels = round2),
    labs(x = axis_labels[match(x_axis, names(axis_labels))]),
    labs(y = axis_labels[match(y_axis, names(axis_labels))]),
    if(!is.na(color)) {
      labs(color = axis_labels[match(color, names(axis_labels))])
    },
    if(!is.null(hline)) {geom_hline(yintercept = hline, alpha = 0.5)},
    if(caption) {labs(caption = caption_builder(condition_at, caption_width))}
  )
}

# Build caption for plot
caption_builder <- function(condition_at, caption_width) {
  if(length(condition_at) == 0L) {
    return(NULL)
  }
  paste0(
    "Conditioned on: ",
    csw(
      paste0(names(condition_at), "=", round2(condition_at)),
      quote = FALSE,
      and = TRUE,
      period = TRUE,
      new_line = Inf
    )
  ) |>
    str_wrap(width = caption_width, collapse = "\n")
}

# Add rounding to ggplot2::label_both() because we don't want long decimals
label_both2 <- function(labels, sep = ": ") {
  value <- lapply(labels, round2)
  variable <- lapply(
    as.list(names(labels)),
    rep,
    if(is.null(nrow(labels))) length(labels[[1]]) else nrow(labels)
  )
  out <- vector("list", length(value))
  for(i in seq_along(out)) {
    out[[i]] <- paste(variable[[i]], value[[i]], sep = sep)
  }
  out
}
attr(label_both2, "class") <- c("function", "labeller")

remove_used_axis <- function(x) {
 if(length(x) > 1L) {
   x <- x[-1L]
 } else {
   x <- NA
 }
  x
}

utils::globalVariables(c(".data", "nsims"))
