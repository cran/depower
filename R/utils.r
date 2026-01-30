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

#' Beta-Binomial quantile function
#'
#' @description
#' Computes quantiles of the Beta-Binomial distribution.
#'
#' @details
#' Returns the smallest integer \eqn{k} such that \eqn{P(Y \leq k) \geq p},
#' where \eqn{Y \sim \text{BetaBinomial}(n, \alpha, \beta)}.
#'
#' The Beta-Binomial PMF is:
#' \deqn{
#'   P(Y = k) = \binom{n}{k} \frac{B(k + \alpha, n - k + \beta)}{B(\alpha, \beta)}
#' }
#'
#' This implementation uses a recurrence relation for efficiency:
#'
#' \deqn{
#'   \frac{P(k)}{P(k-1)} = \frac{n-k+1}{k} \cdot \frac{k - 1 + \alpha}{(n - k) + \beta}
#' }
#'
#' or it's equivalent shifted form
#'
#' \deqn{
#'   \frac{P(k+1)}{P(k)} = \frac{n-k}{k+1} \cdot \frac{k + \alpha}{(n - k - 1) + \beta}
#' }
#'
#' @param p
#' (numeric: `[0, 1]`)\cr
#' Probability or vector of probabilities.
#'
#' @param size
#' (integer: `[0, Inf)`)\cr
#' Number of trials (n) or vector of trial counts.
#'
#' @param shape1
#' (numeric: `(0, Inf)`)\cr
#' First shape parameter \eqn{\alpha} of the Beta distribution, or vector.
#'
#' @param shape2
#' (numeric: `(0, Inf)`)\cr
#' Second shape parameter \eqn{\beta} of the Beta distribution, or vector.
#'
#' @returns
#' Integer vector of quantiles, same length as the longest input.
#'
#' @examples
#' library(depower)
#'
#' # Single quantile
#' depower:::qbetabinom(0.5, size = 100, shape1 = 2, shape2 = 2)
#'
#' # Multiple quantiles
#' depower:::qbetabinom(c(0.025, 0.975), size = 100, shape1 = 10, shape2 = 5)
#'
#' # Vectorized over all arguments
#' depower:::qbetabinom(
#'   p = c(0.025, 0.975),
#'   size = c(100, 200),
#'   shape1 = c(10, 20),
#'   shape2 = c(5, 10)
#' )
#'
#' @noRd
qbetabinom <- function(p, size, shape1, shape2) {
  #-----------------------------------------------------------------------------
  # Check arguments
  #-----------------------------------------------------------------------------
  if (anyNA(p) || anyNA(size) || anyNA(shape1) || anyNA(shape2)) {
    stop("NA values not allowed in arguments.")
  }
  if (any(size < 0)) {
    stop("Argument 'size' must be non-negative.")
  }
  if (any(shape1 <= 0) || any(shape2 <= 0)) {
    stop("Arguments 'shape1' and 'shape2' must be positive.")
  }

  #-----------------------------------------------------------------------------
  # Vectorize inputs
  #-----------------------------------------------------------------------------
  # vector names are used in the error messages
  lens <- c(
    p = length(p),
    size = length(size),
    shape1 = length(shape1),
    shape2 = length(shape2)
  )
  max_len <- max(lens)

  if (max_len == 1L) {
    return(qbetabinom_scalar(p, size, shape1, shape2))
  }

  # Check for lengths that don't evenly divide max_len (would cause incomplete
  # recycling). In contrast to e.g. qbinom() which silently recycles...
  bad <- lens > 1L & max_len %% lens != 0L
  if (any(bad)) {
    stop(
      "Lengths of arguments are not compatible because\nthey do not all evenly divide with the maximum length:\n",
      paste0(names(lens), ": ", lens, collapse = "\n")
    )
  }

  p <- rep_len(p, max_len)
  size <- rep_len(size, max_len)
  shape1 <- rep_len(shape1, max_len)
  shape2 <- rep_len(shape2, max_len)

  vapply(
    X = seq_len(max_len),
    FUN = function(i) qbetabinom_scalar(p[i], size[i], shape1[i], shape2[i]),
    FUN.VALUE = integer(1L)
  )
}

#' Scalar Beta-Binomial quantile function
#'
#' The non-vectorized backend of `qbetabinom()`.
#'
#' @inheritParams qbetabinom
#'
#' @returns
#' Single integer quantile.
#'
#' @noRd
qbetabinom_scalar <- function(p, size, shape1, shape2) {
  n <- size
  alpha <- shape1
  beta <- shape2

  #-----------------------------------------------------------------------------
  # Edge cases
  #-----------------------------------------------------------------------------
  if (p <= 0) {
    return(0L)
  }
  if (p >= 1) {
    return(as.integer(n))
  }

  #-----------------------------------------------------------------------------
  # Normal approximation
  #-----------------------------------------------------------------------------
  # Normal approximation only for roughly symmetric cases
  # ok <- n > 100 &&
  #   min(alpha, beta) > 5
  # if (ok) {
  #   ab_sum <- alpha + beta
  #   mu <- n * alpha / ab_sum
  #   sigma <- sqrt(n * alpha * beta * (ab_sum + n) / (ab_sum^2 * (ab_sum + 1)))
  #   k_approx <- stats::qnorm(p, mu, sigma)
  #   return(as.integer(max(0L, min(n, round(k_approx)))))
  # }

  #-----------------------------------------------------------------------------
  # If P(Y=0) underflows to 0 (common when size is large and/or shape1 is small)
  # the recurrence stays at 0 and the CDF never increases. In this case, the
  # loop finishes and returns size, which can be incorrect for p < 1.
  # Instead, use slow but more robust calculation here.
  #-----------------------------------------------------------------------------
  log_pmf_0 <- lbeta(alpha, n + beta) - lbeta(alpha, beta)
  if (log_pmf_0 < log(.Machine$double.xmin)) {
    # P(Y=0) underflows, so start CDF at k=1. quantile cannot be 0 here.
    k <- seq_len(n)
    log_pmf <- lchoose(n, k) +
      lbeta(k + alpha, n - k + beta) -
      lbeta(alpha, beta)

    # m is the largest log-probability. Subtracting m shifts all log-probs so
    # the largest becomes 0. This guarantees exp(log_pmf - m) <= 1, preventing
    # overflow, and makes very small terms less likely to underflow to 0.
    m <- max(log_pmf)
    # These are scaled weights, proportional to the original probabilities:
    # w_i = exp(log_pmf_i - m) = exp(log_pmf_i) / exp(m).
    # They have the same ratios as the original pmf values, just rescaled.
    w <- exp(log_pmf - m)
    # Normalization constant on the scaled scale:
    # If we multiplied back by exp(m), we'd recover the original sum of pmf's.
    Z <- sum(w)
    # The normalized cumulative distribution, computed stably.
    cdf <- cumsum(w) / Z
    return(as.integer(min(k[cdf >= p])))
  }

  #-----------------------------------------------------------------------------
  # Sequential computation with recurrence relation
  #
  # Recurrence: P(k+1) / P(k) = [(n-k) / (k+1)] * [(k+a) / (n-k-1+b)]
  # Equivalently: P(k) / P(k-1) = [(n-k+1) / k] * [(k-1+a) / (n-k+b)]
  #
  # This avoids computing lbeta() for every k, which is the main cost.
  # Early stopping when CDF reaches target probability.
  #-----------------------------------------------------------------------------
  # Compute P(Y = 0) using log-scale for numerical stability
  pmf <- exp(log_pmf_0)
  cdf <- pmf

  if (cdf >= p) {
    return(0L)
  }

  # Sequential search with early stopping
  for (k in seq_len(n)) {
    ratio <- ((n - k + 1) / k) * ((k - 1 + alpha) / (n - k + beta))
    pmf <- pmf * ratio
    cdf <- cdf + pmf

    if (cdf >= p) {
      return(as.integer(k))
    }
  }

  as.integer(n)
}

#' @title
#' Clopper-Pearson exact confidence interval for a binomial proportion
#'
#' @description
#' Calculates the Clopper-Pearson exact confidence interval for a binomial proportion.
#'
#' @details
#' The Clopper-Pearson exact interval inverts the binomial test via Beta quantiles.
#' The bounds \eqn{(\pi_L, \pi_U)} satisfy:
#'
#' \deqn{P(X \geq x \mid \pi = \pi_L) = \alpha/2}
#' \deqn{P(X \leq x \mid \pi = \pi_U) = \alpha/2}
#'
#' With \eqn{x} successes in \eqn{n} trials,
#'
#' \deqn{\pi_L = B^{-1}\left(\frac{\alpha}{2}; x, n-x+1\right)}
#' \deqn{\pi_U = B^{-1}\left(1-\frac{\alpha}{2}; x+1, n-x\right)}
#'
#' where \eqn{B^{-1}(q; a, b)} is the \eqn{q}-th quantile of
#' \eqn{\text{Beta}(a, b)}.
#'
#' This method guarantees at least nominal coverage but is conservative
#' (intervals are wider than necessary).
#'
#' @references
#' \insertRef{newcombe_1998}{depower},
#'
#' \insertRef{wilson_1927}{depower},
#'
#' \insertRef{clopper_1934}{depower}
#'
#' @param x
#' (integer: `[0, Inf)`)\cr
#' Number of successes.
#' Must be non-negative integers with `x <= n`.
#'
#' @param n
#' (integer: `[1, Inf)`)\cr
#' Number of trials.
#' Must be positive integers.
#'
#' @param conf_level
#' (Scalar numeric: `0.95`; `(0,1)`)\cr
#' The confidence level.
#'
#' @returns
#' A list with elements:
#' \tabular{ll}{
#'   Name \tab Description \cr
#'   `lower` \tab Lower bound of the confidence interval. \cr
#'   `upper` \tab Upper bound of the confidence interval.
#' }
#'
#' @seealso
#' [depower::eval_power_ci()],
#' [depower::add_power_ci()],
#' [depower::binom_ci_wilson()],
#' [depower::binom_pi_bayes()]
#'
#' @examples
#' #----------------------------------------------------------------------------
#' # binom_ci_clopper_pearson() examples
#' #----------------------------------------------------------------------------
#' library(depower)
#'
#' # Single proportion: 80 successes out of 100 trials
#' binom_ci_clopper_pearson(x = 80, n = 100)
#'
#' # Vectorized: multiple proportions
#' binom_ci_clopper_pearson(x = c(8, 80, 800), n = c(10, 100, 1000))
#'
#' # 99% confidence interval
#' binom_ci_clopper_pearson(x = 80, n = 100, conf_level = 0.99)
#'
#' @noRd
binom_ci_clopper_pearson <- function(x, n, conf_level = 0.95) {
  #-----------------------------------------------------------------------------
  # Check arguments
  #-----------------------------------------------------------------------------
  if (!is.numeric(x) || any(x < 0) || any(x != round(x))) {
    stop("Argument 'x' must be non-negative integers.")
  }

  if (!is.numeric(n) || any(n <= 0) || any(n != round(n))) {
    stop("Argument 'n' must be positive integers.")
  }

  if (length(x) != length(n)) {
    stop("Arguments 'x' and 'n' must be the same length.")
  }

  if (any(x > n)) {
    stop("Argument 'x' must be less than or equal to 'n'.")
  }

  if (
    !is.numeric(conf_level) ||
      length(conf_level) != 1L ||
      conf_level <= 0 ||
      conf_level >= 1
  ) {
    stop("Argument 'conf_level' must be a scalar numeric in (0, 1).")
  }

  #-----------------------------------------------------------------------------
  # Compute Clopper-Pearson exact interval
  #-----------------------------------------------------------------------------
  alpha <- 1 - conf_level

  lower <- ifelse(
    x == 0,
    0,
    stats::qbeta(alpha / 2, x, n - x + 1)
  )

  upper <- ifelse(
    x == n,
    1,
    stats::qbeta(1 - alpha / 2, x + 1, n - x)
  )

  #-----------------------------------------------------------------------------
  # Return
  #-----------------------------------------------------------------------------
  list(
    lower = lower,
    upper = upper
  )
}

#' @title
#' Wilson score confidence interval for a binomial proportion
#'
#' @description
#' Calculates the Wilson score confidence interval for a binomial proportion.
#'
#' @details
#' The Wilson score interval is derived from inverting the score test.
#' Starting with the inequality
#'
#' \deqn{
#'   \left| \frac{\hat{\pi}-\pi}{\sqrt{\pi(1-\pi)/n}} \right| \le z_{1-\alpha/2},
#' }
#'
#' and solving the resulting quadratic for \eqn{\pi} yields
#'
#' \deqn{
#'   \frac{\hat{\pi}+\frac{z^2}{2n} \pm z \sqrt{\frac{\hat{\pi}(1-\hat{\pi})}{n}+\frac{z^2}{4n^2}}}{1+\frac{z^2}{n}},
#' }
#'
#' with \eqn{z = z_{1-\alpha/2}} and \eqn{\hat{\pi} = x/n}.
#'
#' @references
#' \insertRef{newcombe_1998}{depower},
#'
#' \insertRef{wilson_1927}{depower},
#'
#' \insertRef{clopper_1934}{depower}
#'
#' @param x
#' (integer: `[0, Inf)`)\cr
#' Number of successes.
#' Must be non-negative integers with `x <= n`.
#'
#' @param n
#' (integer: `[1, Inf)`)\cr
#' Number of trials.
#' Must be positive integers.
#'
#' @param conf_level
#' (Scalar numeric: `0.95`; `(0,1)`)\cr
#' The confidence level.
#'
#' @returns
#' A list with elements:
#' \tabular{ll}{
#'   Name \tab Description \cr
#'   `lower` \tab Lower bound of the confidence interval. \cr
#'   `upper` \tab Upper bound of the confidence interval.
#' }
#'
#' @seealso
#' [depower::eval_power_ci()],
#' [depower::add_power_ci()],
#' [depower::binom_ci_clopper_pearson()],
#' [depower::binom_pi_bayes()]
#'
#' @examples
#' #----------------------------------------------------------------------------
#' # binom_ci_wilson() examples
#' #----------------------------------------------------------------------------
#' library(depower)
#'
#' # Single proportion: 80 successes out of 100 trials
#' binom_ci_wilson(x = 80, n = 100)
#'
#' # Vectorized: multiple proportions
#' binom_ci_wilson(x = c(8, 80, 800), n = c(10, 100, 1000))
#'
#' # 99% confidence interval
#' binom_ci_wilson(x = 80, n = 100, conf_level = 0.99)
#'
#' @noRd
binom_ci_wilson <- function(x, n, conf_level = 0.95) {
  #-----------------------------------------------------------------------------
  # Check arguments
  #-----------------------------------------------------------------------------
  if (!is.numeric(x) || any(x < 0) || any(x != round(x))) {
    stop("Argument 'x' must be non-negative integers.")
  }

  if (!is.numeric(n) || any(n <= 0) || any(n != round(n))) {
    stop("Argument 'n' must be positive integers.")
  }

  if (length(x) != length(n)) {
    stop("Arguments 'x' and 'n' must be the same length.")
  }

  if (any(x > n)) {
    stop("Argument 'x' must be less than or equal to 'n'.")
  }

  if (
    !is.numeric(conf_level) ||
      length(conf_level) != 1L ||
      conf_level <= 0 ||
      conf_level >= 1
  ) {
    stop("Argument 'conf_level' must be a scalar numeric in (0, 1).")
  }

  #-----------------------------------------------------------------------------
  # Compute Wilson score interval
  #-----------------------------------------------------------------------------
  alpha <- 1 - conf_level
  p_hat <- x / n
  z <- stats::qnorm(1 - alpha / 2)

  denom <- 1 + z^2 / n
  center <- (p_hat + z^2 / (2 * n)) / denom
  margin <- z * sqrt(p_hat * (1 - p_hat) / n + z^2 / (4 * n^2)) / denom

  lower <- center - margin
  upper <- center + margin

  #-----------------------------------------------------------------------------
  # Return
  #-----------------------------------------------------------------------------
  list(
    lower = lower,
    upper = upper
  )
}

#' @title
#' Bayesian posterior predictive interval for a binomial proportion
#'
#' @description
#' Calculates the Bayesian posterior predictive interval for a binomial proportion.
#' The interval quantifies the expected range of proportion estimates from a future study.
#'
#' @details
#' With a \eqn{\text{Beta}(\alpha, \beta)} prior on the true proportion \eqn{\pi},
#' the posterior after observing \eqn{x} successes in \eqn{n} trials is:
#' \deqn{
#'   \pi \mid X = x \sim \text{Beta}(\alpha + x, \beta + n - x)
#' }
#'
#' The posterior predictive distribution for \eqn{Y}, the number of successes in a future
#' study with \eqn{m} trials, is Beta-Binomial:
#' \deqn{
#'   Y \mid X = x \sim \text{BetaBinomial}(m, \alpha + x, \beta + n - x)
#' }
#'
#' The posterior predictive interval is constructed from quantiles of this
#' distribution, expressed as proportions \eqn{Y/m}.
#'
#' The posterior predictive mean and variance of \eqn{\hat{\pi}_{\text{new}} = Y/m} are:
#' \deqn{
#' \begin{aligned}
#'   E[\hat{\pi}_{\text{new}} \mid X = x] &= \frac{\alpha + x}{\alpha + \beta + n} \\
#'   \text{Var}[\hat{\pi}_{\text{new}} \mid X = x]
#'   &= \frac
#'       {(\alpha + x)(\beta + n - x)(\alpha + \beta + n + m)}
#'       {m (\alpha + \beta + n)^{2} (\alpha + \beta + n + 1)}.
#' \end{aligned}
#' }
#'
#' @references
#' \insertRef{gelman_2013}{depower}
#'
#' @param x
#' (integer: `[0, Inf)`)\cr
#' Number of successes.
#' Must be non-negative integers with `x <= n`.
#'
#' @param n
#' (integer: `[1, Inf)`)\cr
#' Number of trials.
#'
#' @param future_n
#' (integer or `NULL`: `NULL`; `[1, Inf)`)\cr
#' Number of trials in the future study.
#' If `NULL` (default), uses the same number as the observed study (`n`).
#'
#' @param pred_level
#' (Scalar numeric: `0.95`; `(0,1)`)\cr
#' The posterior predictive interval level.
#'
#' @param prior
#' (Numeric vector of length 2: `c(1, 1)`; each `(0, Inf)`)\cr
#' Parameters \eqn{(\alpha, \beta)} for the Beta prior on the true proportion.
#' Default `c(1, 1)` is the uniform prior.
#' Use `c(0.5, 0.5)` for the Jeffreys prior.
#'
#' @returns
#' A list with elements:
#' \tabular{ll}{
#'   Name \tab Description \cr
#'   `mean` \tab Predictive mean of future proportion estimate. \cr
#'   `lower` \tab Lower bound of posterior predictive interval. \cr
#'   `upper` \tab Upper bound of posterior predictive interval.
#' }
#'
#' @seealso
#' [depower::eval_power_pi()],
#' [depower::add_power_pi()],
#' [depower::binom_ci_wilson()],
#' [depower::binom_ci_clopper_pearson()]
#'
#' @examples
#' #----------------------------------------------------------------------------
#' # binom_pi_bayes() examples
#' #----------------------------------------------------------------------------
#' library(depower)
#'
#' # Single proportion: 80 successes out of 100 trials
#' binom_pi_bayes(x = 80, n = 100)
#'
#' # Predict for a larger future study (narrower interval)
#' binom_pi_bayes(x = 80, n = 100, future_n = 1000)
#'
#' # Predict for a smaller future study (wider interval)
#' binom_pi_bayes(x = 80, n = 100, future_n = 50)
#'
#' # Use Jeffreys prior instead of uniform
#' binom_pi_bayes(x = 80, n = 100, prior = c(0.5, 0.5))
#'
#' # Vectorized: multiple proportions
#' binom_pi_bayes(x = c(8, 80, 800), n = c(10, 100, 1000))
#'
#' # 99% predictive interval
#' binom_pi_bayes(x = 80, n = 100, pred_level = 0.99)
#'
#' @noRd
binom_pi_bayes <- function(
  x,
  n,
  future_n = NULL,
  pred_level = 0.95,
  prior = c(1, 1)
) {
  #-----------------------------------------------------------------------------
  # Check arguments
  #-----------------------------------------------------------------------------
  if (!is.numeric(x) || any(x < 0) || any(x != round(x))) {
    stop("Argument 'x' must be non-negative integers.")
  }

  if (!is.numeric(n) || any(n <= 0) || any(n != round(n))) {
    stop("Argument 'n' must be positive integers.")
  }

  if (length(x) != length(n)) {
    stop("Arguments 'x' and 'n' must be the same length.")
  }

  if (any(x > n)) {
    stop("Argument 'x' must be less than or equal to 'n'.")
  }

  if (!is.null(future_n)) {
    if (
      !is.numeric(future_n) ||
        any(future_n <= 0) ||
        any(future_n != round(future_n))
    ) {
      stop("Argument 'future_n' must be positive integers or NULL.")
    }
    if (length(future_n) != 1L && length(future_n) != length(x)) {
      stop("Argument 'future_n' must be length 1 or the same length as 'x'.")
    }
  }

  if (
    !is.numeric(pred_level) ||
      length(pred_level) != 1L ||
      pred_level <= 0 ||
      pred_level >= 1
  ) {
    stop("Argument 'pred_level' must be a scalar numeric in (0, 1).")
  }

  if (!is.numeric(prior) || length(prior) != 2L || any(prior <= 0)) {
    stop("Argument 'prior' must be a positive numeric vector of length 2.")
  }

  #-----------------------------------------------------------------------------
  # Compute posterior predictive interval
  #-----------------------------------------------------------------------------
  alpha_pred <- 1 - pred_level

  # Future study size: use future_n if provided, otherwise match observed
  m <- if (is.null(future_n)) n else rep_len(future_n, length(x))

  # Posterior parameters: Beta(alpha_post, beta_post)
  alpha_post <- prior[1] + x
  beta_post <- prior[2] + n - x

  # Predictive mean: E[Y/m | X = x]
  pred_mean <- alpha_post / (alpha_post + beta_post)

  # Posterior predictive interval from Beta-Binomial quantiles
  pred_lower <- qbetabinom(alpha_pred / 2, m, alpha_post, beta_post) / m
  pred_upper <- qbetabinom(1 - alpha_pred / 2, m, alpha_post, beta_post) / m

  #-----------------------------------------------------------------------------
  # Return
  #-----------------------------------------------------------------------------
  list(
    mean = pred_mean,
    lower = pred_lower,
    upper = pred_upper,
    future_n = m
  )
}
