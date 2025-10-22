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
  if (any(duplicated.default(dots))) {
    stop("Argument '...' must not contain duplicated calls.")
  }

  # Process names. Check for missings or duplications.
  if (.names) {
    dots_names <- names(dots)
    if (any(dots_names == "") || is.null(dots_names)) {
      are_missing <- if (is.null(dots_names)) {
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
    if (any(duplicated.default(dots_names))) {
      switch(
        .duplicate_names,
        "warn" = warning(
          "Argument '...' must not contain duplicated names of calls."
        ),
        "stop" = stop(
          "Argument '...' must not contain duplicated names of calls."
        )
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
# Convenience function to extract elements from a list-column.
#-------------------------------------------------------------------------------
## data = data frame
## column = list-column name
## element = list-column list element name
## value = template (type) of the return value
extract <- function(data, column, element, value) {
  vapply(
    X = data[[column]],
    FUN = function(x) x[[element]],
    FUN.VALUE = value
  )
}

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
      if (x >= 1 / (10^digits)) {
        x <- round(x, digits)
      } else if (x > 0 && x < 1 / (10^digits)) {
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
  if (quote) {
    x <- paste0("'", x, "'", recycle0 = TRUE)
  }
  if (length(x) > n && ellipsis) {
    x <- c(x[seq_len(n)], "...")
    and <- FALSE
    period <- FALSE
  }
  if (length(x) == 2) {
    if (and) {
      x <- paste0(x, collapse = " and ")
    } else {
      x <- paste0(x, collapse = ", ")
    }
  }
  if (length(x) > 2) {
    x <- paste0(x, collapse = ", ")
    if (and) {
      x <- sub("(.*), ", "\\1, and ", x)
    }
  }
  if (is.finite(new_line)) {
    x <- str_wrap(x, width = new_line, collapse = "\n")
  }
  if (period) {
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
# Convert string to title case (first character)
#-------------------------------------------------------------------------------
str_to_title <- function(x) {
  substr(x, 1L, 1L) <- toupper(substr(x, 1L, 1L))
  x
}

#-------------------------------------------------------------------------------
# Fast conversion from simulated list to tall data frame.
#-------------------------------------------------------------------------------
list_to_df <- function(x) {
  if (length(x) < 2L) {
    # One sample case
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
  } else {
    # Two sample case
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

#-------------------------------------------------------------------------------
# Functions for randomization tests and simulation of test statistic null
# distributions
#-------------------------------------------------------------------------------
pchisq2 <- function(
  q,
  df = 1L,
  lower.tail = FALSE,
  q_null = NULL
) {
  # Check arguments
  if (!is.null(q_null)) {
    if (!(is.numeric(q_null) && length(q_null) > 1L)) {
      stop("Argument 'q_null' must be a numeric vector.")
    }
    if (anyNA(q_null)) {
      q_null <- remove_na(q_null)
      warning(
        "The simulated null chi-square test statistic distribution had missing values."
      )
    }
  }

  # Calculate p-value
  # \insertRef{phipson_2010}{depower}
  p <- if (is.null(q_null)) {
    pchisq(q = q, df = df, lower.tail = lower.tail)
  } else {
    if (lower.tail) {
      (sum(q_null <= q) + 1L) / (length(q_null) + 1L)
    } else {
      (sum(q_null >= q) + 1L) / (length(q_null) + 1L)
    }
  }

  # Return
  p
}

# Perform the randomization tests
randomize_tests <- function(data, distribution, paired, call) {
  # Check arguments
  if (!(is.list(data) && length(data) == 2L)) {
    stop("Argument 'data' must be a list with 2 elements.")
  }

  check_simulated(distribution)
  ncores <- distribution$ncores
  if (ncores > 1L) {
    if (isTRUE(ncores > parallel::detectCores())) {
      max <- parallel::detectCores()
      warning(
        "Argument 'ncores' in 'depower::simulated()' should not be greater than ",
        max
      )
    }
  }

  # Simulate randomized data
  data_rand <- if (paired) {
    switch(
      distribution$method,
      "approximate" = randomize_dependent_two_sample(
        data = data,
        distribution = distribution
      ),
      "exact" = permute_dependent_two_sample(
        data = data,
        distribution = distribution
      ),
      stop("Argument 'distribution' had unknown value for '$method'.")
    )
  } else {
    switch(
      distribution$method,
      "approximate" = randomize_independent_two_sample(
        data = data,
        distribution = distribution
      ),
      "exact" = permute_independent_two_sample(
        data = data,
        distribution = distribution
      ),
      stop("Argument 'distribution' had unknown value for '$method'.")
    )
  }

  # Create data for rowwise operations
  res <- cbind(
    simulation = seq_len(length(data_rand)),
    list2DF(list(data = data_rand))
  )

  # Ensure call has correct arguments
  call$data <- quote(.data)
  call$distribution <- asymptotic()
  call$ci_level <- NULL

  # Run tests
  if (ncores > 1L) {
    cluster <- multidplyr::new_cluster(ncores)
    multidplyr::cluster_library(cluster, 'depower')
  }
  res <- res |>
    dplyr::rowwise() |>
    {
      \(.) {
        if (ncores > 1L) {
          multidplyr::partition(data = ., cluster = cluster)
        } else {
          .
        }
      }
    }() |>
    dplyr::mutate(
      result = list(eval(
        expr = call,
        envir = list(.data = data),
        enclos = parent.frame()
      ))
    ) |>
    {
      \(.) {
        if (ncores > 1L) {
          dplyr::collect(x = .)
        } else {
          .
        }
      }
    }() |>
    dplyr::ungroup()

  # Return
  res
}

# Independent data randomization
randomize_independent_two_sample <- function(data, distribution) {
  # Check arguments
  n1 <- length(data[[1L]])
  n2 <- length(data[[2L]])

  # Randomly assign groups to combined data.
  data_combined <- c(data[[1L]], data[[2L]])
  out <- lapply(
    X = seq_len(distribution$nsims),
    FUN = function(x) {
      rand <- sample.int(n1 + n2, replace = FALSE)
      list(
        value1 = data_combined[rand[seq_len(n1)]],
        value2 = data_combined[rand[-seq_len(n1)]]
      )
    }
  )

  # Return
  c(list(data), out)
}

# Dependent data randomization
randomize_dependent_two_sample <- function(data, distribution) {
  # Check arguments
  value1 <- data[[1L]]
  value2 <- data[[2L]]

  n1 <- length(value1)
  n2 <- length(value2)

  if (n1 != n2) {
    stop("Argument 'data' must have the same sample size for both samples.")
  }

  if (anyNA(value1) || anyNA(value2)) {
    not_na <- complete.cases(value1, value2)
    value1 <- value1[not_na]
    value2 <- value2[not_na]

    n1 <- length(value1)
    n2 <- length(value2)
  }

  # Randomize data pairs. i.e. randomly switch the values for each pair.
  # Below doesn't keep original order.
  data_combined <- c(value1, value2)
  out <- lapply(
    X = seq_len(distribution$nsims),
    FUN = function(x) {
      rand <- sample(x = c(TRUE, FALSE), size = n1, replace = TRUE)
      rand_value1 <- c(which(rand), n1 + which(!rand))
      rand_value2 <- c(n1 + which(rand), which(!rand))
      out <- list(
        value1 = data_combined[rand_value1],
        value2 = data_combined[rand_value2]
      )
    }
  )

  # Return
  c(list(data), out)
}

# Independent data permutation
#' @importFrom utils combn
permute_independent_two_sample <- function(data, distribution) {
  # Check arguments
  n1 <- length(data[[1L]])
  n2 <- length(data[[2L]])
  ncombs <- choose(n1 + n2, n1)

  if (ncombs > 1e6) {
    warning(
      paste0(
        "Using approximate randomization with ",
        distribution$nsims,
        " resamples. ",
        "The exact randomization test is only used when choose(",
        n1,
        "+",
        n2,
        ",",
        n1,
        ")=",
        formatC(ncombs, format = "e", digits = 2),
        " < 1e6."
      )
    )
    return(randomize_independent_two_sample(
      data = data,
      distribution = distribution
    ))
  }

  # Assign groups for all combn(n1+n2, n1) combinations.
  data_combined <- c(data[[1L]], data[[2L]])
  combs <- combn(x = n1 + n2, m = n1, FUN = NULL, simplify = FALSE)
  out <- lapply(combs, function(x) {
    list(
      value1 = data_combined[x],
      value2 = data_combined[-x]
    )
  })

  # Return
  out
}

# Dependent data permutation
# Can be thought of as all combn(x,m) for each value of m.
# x is seq_len(n_pairs)
powerset <- function(x) {
  x_len <- length(x)
  out <- vector(mode = "list", length = 2^x_len)
  out[[1]] <- integer()
  k <- 1L
  for (i in seq_len(x_len)) {
    for (j in seq_len(k)) {
      k <- k + 1L
      out[[k]] <- c(out[[j]], x[i])
    }
  }
  out
}

permute_dependent_two_sample <- function(data, distribution) {
  # Check arguments
  value1 <- data[[1L]]
  value2 <- data[[2L]]

  n1 <- length(value1)
  n2 <- length(value2)

  if (n1 != n2) {
    stop("Argument 'data' must have the same sample size for both samples.")
  }

  if (anyNA(value1) || anyNA(value2)) {
    not_na <- complete.cases(value1, value2)
    value1 <- value1[not_na]
    value2 <- value2[not_na]
    data <- list(value1 = value1, value2 = value2)

    n1 <- length(value1)
    n2 <- length(value2)
  }

  if (n1 > 20L) {
    warning(
      paste0(
        "Using approximate randomization with ",
        distribution$nsims,
        " resamples. ",
        "The exact randomization test is only used for sample size per group 20 or fewer (max 2^20=1048576 resamples)."
      )
    )
    return(randomize_dependent_two_sample(
      data = data,
      distribution = distribution
    ))
  }

  # Use the powerset (minus null set) to randomize all 2^n data pairs.
  # i.e. randomly switch the values for each pair.
  # Below doesn't keep original order.
  idx <- powerset(seq_len(n1))[-1]
  out <- lapply(idx, function(x) {
    list(
      value1 = c(value2[x], value1[-x]),
      value2 = c(value1[x], value2[-x])
    )
  })

  # Return
  c(list(data), out)
}

#-------------------------------------------------------------------------------
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
  if (!is.logical(warnings) || length(warnings) != 1L) {
    stop("Argument 'warnings' must be a scalar logical.")
  }
  if (!(length(lower) == 1L || length(lower) == length(parameters))) {
    stop("Argument 'lower' must be scalar or length(parameters).")
  }
  if (!(length(upper) == 1L || length(upper) == length(parameters))) {
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
    stop(
      "Argument 'method' must be one of 'nlm', 'nlm_constrained', 'optim', or 'optim_constrained'."
    )
  )

  #-----------------------------------------------------------------------------
  # Unify return object
  #-----------------------------------------------------------------------------
  out <- vector(mode = "list", length = 5L)
  names(out) <- c("estimate", "minimum", "iterations", "code", "message")

  if (method == "nlm") {
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
    if (code > 2L && warnings) {
      warning("The MLE algorithm may not have found a reliable solution.")
    }
  }

  if (method == "nlm_constrained") {
    out[["estimate"]] <- res[["par"]]
    out[["minimum"]] <- res[["objective"]]
    out[["iterations"]] <- res[["iterations"]]
    code <- res[["convergence"]]
    out[["code"]] <- code
    out[["message"]] <- res[["message"]]
    if (code > 0L && warnings) {
      warning("The MLE algorithm may not have found a reliable solution.")
    }
  }

  if (method %in% c("optim", "optim_constrained")) {
    out[["estimate"]] <- res[["par"]]
    out[["minimum"]] <- res[["value"]]
    out[["iterations"]] <- res[["counts"]][[1L]]
    code <- res[["convergence"]]
    out[["code"]] <- code
    if (code > 0L && warnings) {
      warning("The MLE algorithm may not have found a reliable solution.")
    }
    message <- res[["message"]]
    out[["message"]] <- if (is.null(message)) {
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
