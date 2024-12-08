#' Capture unevaluated dots
#'
#' `dots()` returns list of unevaluated expressions.
#'
#' @param ... The `...` arguments passed down from the original call.
#' @param .names (scalar logical: `TRUE`)\cr
#'        Whether or not to fill in missing names slot using a character version
#'        of the expression.
#' @param .duplicate_names (string: `"warn"`)\cr
#'        What should be done if duplicate names are found? Choices are
#'        `"ignore"`, `"warn"`, or `"stop"`. Only used if `.names = TRUE`.
#'
#' @return
#' - `dots()`: list
#'
#' @keywords Internal
#' @noRd
dots <- function(..., .names = TRUE, .duplicate_names = "stop") {
  # Capture the dots
  dots <- as.list(substitute(list(...)))[-1L]

  # Error if duplicated dots
  if(any(duplicated.default(dots))) {
    stop("Argument '...' must not contain duplicated calls.")
  }

  # Process names. Check for missings or duplications.
  if(.names) {
    dots_names <- names(dots)
    if(any(dots_names == "") || is.null(dots_names)) {
      are_missing <- if(is.null(dots_names)) {
        rep(TRUE, length(dots))
      } else {
        dots_names == ""
      }
      dots_char <- vapply(
        X = dots,
        FUN = deparse1,
        FUN.VALUE = NA_character_,
        USE.NAMES = FALSE
      )
      dots_names[are_missing] <- dots_char[are_missing]
      names(dots) <- dots_names
    }

    # Duplicate names
    if(any(duplicated.default(dots_names))) {
      switch(.duplicate_names,
        "warn" = warning("Argument '...' must not contain duplicated names of calls."),
        "stop" = stop("Argument '...' must not contain duplicated names of calls.")
      )
    }
  }

  # Return
  dots
}

#-------------------------------------------------------------------------------
# Don't drop class attribute
#-------------------------------------------------------------------------------
#' @export
`[.depower` <- function(x, ...) {
  res <- NextMethod("[")
  class(res) <- class(x)
  res
}

#-------------------------------------------------------------------------------
# Simple function to negate a logical vector
#-------------------------------------------------------------------------------
negate <- `!`

#-------------------------------------------------------------------------------
# Simple function to remove NAs from an atomic vector
#-------------------------------------------------------------------------------
remove_na <- function(x) x[!is.na(x)]

#-------------------------------------------------------------------------------
# Round numbers in a mixed-type vector. Returns a character vector.
#-------------------------------------------------------------------------------
round2 <- function(x, digits = 2) {
  # Don't want logicals converted to numeric
  x <- as.character(x)

  # Convert numbers back to numeric
  x_num <- suppressWarnings(as.numeric(x))
  is_x_num <- !is.na(x_num)

  # Round numbers
  x[is_x_num] <- vapply(
    X = x_num[is_x_num],
    FUN = function(x) {
      if(x >= 1/(10^digits)) {
        x <- round(x, digits)
      } else if (x > 0 && x < 1/(10^digits)) {
        x <- sprintf("%.1e", x)
      }
      as.character(x)
    },
    FUN.VALUE = character(1L)
  )

  # Return
  x
}

#-------------------------------------------------------------------------------
# Simple function for comma separated words.
#-------------------------------------------------------------------------------
csw <- function(
    x,
    n = Inf,
    quote = TRUE,
    ellipsis = TRUE,
    and = FALSE,
    period = FALSE,
    new_line = 70L
) {
  if(quote) {
    x <- paste0("'", x, "'", recycle0 = TRUE)
  }
  if(length(x) > n && ellipsis) {
    x <- c(x[seq_len(n)], "...")
    and <- FALSE
    period <- FALSE
  }
  if(length(x) == 2) {
    if(and) {
      x <- paste0(x, collapse = " and ")
    } else {
      x <- paste0(x, collapse = ", ")
    }
  }
  if(length(x) > 2) {
    x <- paste0(x, collapse = ", ")
    if(and) {
      x <- sub("(.*), ", "\\1, and ", x)
    }
  }
  if(is.finite(new_line)) {
    x <- str_wrap(x, width = new_line, collapse = "\n")
  }
  if(period) {
    x <- paste0(x, ".")
  }

  # Return
  x
}

#-------------------------------------------------------------------------------
# Wrap character vectors for paragraphs
#-------------------------------------------------------------------------------
str_wrap <- function(x, width = 70, collapse = "\n\n") {
  paste0(strwrap(x, width = width), collapse = collapse)
}

#-------------------------------------------------------------------------------
# Fast conversion from simulated list to tall data frame.
#-------------------------------------------------------------------------------
list_to_df <- function(x) {
  if(length(x) < 2L) { # One sample case
    value <- x[[1L]]
    n <- length(value)
    # data.frame is slow.
    res <- list(
      item = as.character(seq_len(n)),
      value = value
    )
    attr(res, "names") <- c("item", "value")
    attr(res, "row.names") <- .set_row_names(n)
    attr(res, "class") <- "data.frame"
  } else { # Two sample case
    value1 <- x[[1L]]
    value2 <- x[[2L]]
    n1 <- length(value1)
    n2 <- length(value2)
    # factor is slow. Just need to make sure input is integer.
    condition <- c(rep(1L, times = n1), rep(2L, times = n2))
    attr(condition, "levels") <- c("1", "2")
    attr(condition, "class") <- "factor"
    # data.frame is slow.
    res <- list(
      item = as.character(c(seq_len(n1), seq_len(n2))),
      condition = condition,
      value = c(value1, value2)
    )
    attr(res, "names") <- c("item", "condition", "value")
    attr(res, "row.names") <- .set_row_names(n1 + n2)
    attr(res, "class") <- "data.frame"
  }

  # Return
  res
}

#-------------------------------------------------------------------------------
# Fast conversion from two-sample simulated data frame to list.
#-------------------------------------------------------------------------------
# Two-sample case only, assumes sorted by item/subject index
df_to_list <- function(x) {
  idx_value1 <- x[["condition"]] == "1"
  idx_value2 <- !idx_value1
  list(
    x[["value"]][idx_value1],
    x[["value"]][idx_value2]
  )
}

#-------------------------------------------------------------------------------
# Fast summary statistics.
# Manual calculation should be fine because we are working with clean simulated
# data. values are not extremely small or large, the sample size is small, and
# the standard deviation is not extremely small relative to the mean.
# These are faster than base and reduce external dependencies.
#-------------------------------------------------------------------------------
fmean <- function(x, n = length(x)) {
  sum(x) / n
}

fvar <- function(x, n = length(x), mean = fmean(x)) {
  sum((x - mean)^2) / (n - 1L)
}

# Middle value of a sorted vector
middle <- function(x) {
  sort(x)[ceiling(length(x) / 2)]
}
# keeps data frame labels during mutate
middle.labelled <- function(x) {
  label <- attr(x, "label", exact = TRUE)
  res <- sort(x)[ceiling(length(x) / 2)]
  attr(res, "label") <- label
  res
}

#' A unified interface to four optimization methods.
#'
#' @param nll (function)\cr
#'        A function which returns the negative log likelihood.
#' @param parameters (numeric)\cr
#'        A numeric vector of initial values for each parameter to be optimized
#'        in `nll`.
#' @param method (string: `"optim"`)\cr
#'        A string for the optimization method. Must be from `"nlm"`,
#'        `"nlm_constrained"`, `"optim"`, or `"optim_constrained"`.
#' @param lower (Scalar numeric: `-Inf`)\cr
#'        A scalar numeric for the lower bound.
#' @param upper (Scalar numeric: `Inf`)\cr
#'        A scalar numeric for the upper bound.
#' @param ... Additional arguments passed to the corresponding `method`.
#'
#' @return A list with 5 elements
#'
#' 1. estimate
#' 2. minimum
#' 3. iterations
#' 4. code
#' 5. message
#'
#' @keywords Internal
#' @noRd
mle <- function(
    method,
    nll,
    parameters,
    lower = 1e-03,
    upper = 1e06,
    warnings = FALSE,
    ...
  ) {
  if(!is.logical(warnings) || length(warnings) != 1L) {
    stop("Argument 'warnings' must be a scalar logical.")
  }
  if(!(length(lower) == 1L || length(lower) == length(parameters))) {
    stop("Argument 'lower' must be scalar or length(parameters).")
  }
  if(!(length(upper) == 1L || length(upper) == length(parameters))) {
    stop("Argument 'upper' must be scalar or length(parameters).")
  }

  #-----------------------------------------------------------------------------
  # Run optimization
  #-----------------------------------------------------------------------------
  res <- switch(
    method,
    nlm = nlm(
      f = nll,
      p = parameters,
      ...
    ),
    nlm_constrained = nlminb(
      start = parameters,
      objective = nll,
      ...,
      lower = lower,
      upper = upper
    ),
    optim = optim(
      par = parameters,
      fn = nll,
      ...
    ),
    optim_constrained = optim(
      par = parameters,
      fn = nll,
      ...,
      method = "L-BFGS-B",
      lower = lower,
      upper = upper
    ),
    stop("Argument 'method' must be one of 'nlm', 'nlm_constrained', 'optim', or 'optim_constrained'.")
  )

  #-----------------------------------------------------------------------------
  # Unify return object
  #-----------------------------------------------------------------------------
  out <- vector(mode = "list", length = 5L)
  names(out) <- c("estimate", "minimum", "iterations", "code", "message")

  if(method == "nlm") {
    out[["estimate"]] <- res[["estimate"]]
    out[["minimum"]] <- res[["minimum"]]
    out[["iterations"]] <- res[["iterations"]]
    code <- res[["code"]]
    out[["code"]] <- code
    out[["message"]] <- switch(
      as.character(code),
      "1" = "Relative gradient is close to zero, found probable solution.",
      "2" = "Successive iterates within tolerance, found probable solution.",
      "3" = "Last global step failed to locate a point lower than estimate. Either estimate is an approximate local minimum of the function or steptol is too small.",
      "4" = "Iteration limit exceeded.",
      "5" = "Maximum step size stepmax exceeded five consecutive times. Either the function is unbounded below, becomes asymptotic to a finite value from above in some direction or stepmax is too small.",
      "Wasn't able to determine optimization message"
    )
    if(code > 2L && warnings) {
      warning("The MLE algorithm may not have found a reliable solution.")
    }
  }

  if(method == "nlm_constrained") {
    out[["estimate"]] <- res[["par"]]
    out[["minimum"]] <- res[["objective"]]
    out[["iterations"]] <- res[["iterations"]]
    code <- res[["convergence"]]
    out[["code"]] <- code
    out[["message"]] <- res[["message"]]
    if(code > 0L && warnings) {
      warning("The MLE algorithm may not have found a reliable solution.")
    }
  }

  if(method %in% c("optim", "optim_constrained")) {
    out[["estimate"]] <- res[["par"]]
    out[["minimum"]] <- res[["value"]]
    out[["iterations"]] <- res[["counts"]][[1L]]
    code <- res[["convergence"]]
    out[["code"]] <- code
    if(code > 0L && warnings) {
      warning("The MLE algorithm may not have found a reliable solution.")
    }
    message <- res[["message"]]
    out[["message"]] <- if(is.null(message)) {
      switch(
        as.character(code),
        "0" = "Found probable solution.",
        "1" = "Iteration limit 'maxit' had been reached.",
        "10" = "Degeneracy of the Nelder-Mead simplex",
        "Wasn't able to determine optimization message"
      )
    } else {
      message
    }
  }

  #-----------------------------------------------------------------------------
  # Return
  #-----------------------------------------------------------------------------
  out
}
