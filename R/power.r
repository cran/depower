#' Simulated power
#'
#' A method to calculate power for objects returned by [depower::sim_log_lognormal()],
#' [depower::sim_nb()], and [depower::sim_bnb()].
#'
#' @details
#' Power is calculated as the proportion of hypothesis tests which result in a
#' p-value less than or equal to `alpha`. e.g.
#'
#' ```{r, eval = FALSE}
#' sum(p <= alpha) / nsims
#' ```
#'
#' Power is defined as the expected probability of rejecting the null hypothesis
#' for a chosen value of the unknown effect. In a multiple comparisons scenario,
#' power is defined as the marginal power, which is the expected power of the
#' test for each individual null hypothesis assumed to be false.
#'
#' Other forms of power under the multiple comparisons scenario include
#' disjunctive or conjunctive power. Disjunctive power is defined as the
#' expected probability of correctly rejecting one or more null hypotheses.
#' Conjunctive power is defined as the expected probability of correctly
#' rejecting all null hypotheses. In the simplest case, and where all hypotheses
#' are independent, if the marginal power is defined as \eqn{\pi} and \eqn{m} is
#' the number of null hypotheses assumed to be false, then disjunctive power may
#' be calculated as \eqn{1 - (1 - \pi)^m} and conjunctive power may be calculated
#' as \eqn{\pi^m}. Disjunctive power tends to decrease with increasingly
#' correlated hypotheses and conjunctive power tends to increase with
#' increasingly correlated hypotheses.
#'
#' ## Argument `...`
#' `...` are the name-value pairs for the functions used to perform the tests.
#' If not named, the functions coerced to character will be used for the
#' name-value pairs. Typical in non-standard evaluation, `...` accepts bare
#' functions and converts them to a list of expressions. Each element in this
#' list will be validated as a `call` and then evaluated on the simulated data.
#' A [base::call()] is simply an unevaluated function. Below are some examples
#' of specifying `...` in `power()`.
#'
#' ```{r, eval=FALSE}
#' # Examples of specifying ... in power()
#' data <- sim_nb(
#'   n1 = 10,
#'   mean1 = 10,
#'   ratio = c(1.6, 2),
#'   dispersion1 = 2,
#'   dispersion2 = 2,
#'   nsims = 200
#' )
#'
#' # ... is empty, so an appropriate default function will be provided
#' power(data)
#'
#' # This is equivalent to leaving ... empty
#' power(data, "NB Wald test" = wald_test_nb())
#'
#' # If not named, "wald_test_nb()" will be used to label the function
#' power(data, wald_test_nb())
#'
#' # You can specify any parameters in the call. The data argument
#' # will automatically be inserted or overwritten.
#' data |>
#'   power("NB Wald test" = wald_test_nb(equal_dispersion=TRUE, link="log"))
#'
#' # Multiple functions may be used.
#' data |>
#'   power(
#'     wald_test_nb(link='log'),
#'     wald_test_nb(link='sqrt'),
#'     wald_test_nb(link='squared'),
#'     wald_test_nb(link='identity')
#'   )
#'
#' # Just like functions in a pipe, the parentheses are required.
#' # This will error because wald_test_nb is missing parentheses.
#' try(power(data, wald_test_nb))
#' ```
#'
#' In most cases*, any user created test function may be utilized in `...` if the
#' following conditions are satisfied:
#'
#' 1. The function contains argument `data` which is defined as a list with the
#'    first and second elements for simulated data.
#' 2. The return object is a list with element `p` for the p-value of the
#'    hypothesis test.
#'
#' Validate with test cases beforehand.
#'
#' *Simulated data of class `log_lognormal_mixed_two_sample` has both independent
#' and dependent data. To ensure the appropriate test function is used,
#' `power.log_lognormal_mixed_two_sample()` allows only
#' [depower::t_test_welch()] and [depower::t_test_paired()] in `...`. Each will
#' be evaluated on the simulated data according to column `data$cor`. If one or
#' both of these functions are not included in `...`, the corresponding default
#' function will be used automatically. If any other test function is included,
#' an error will be returned.
#'
#' ## Argument `alpha`
#' \eqn{\alpha} is known as the type I assertion probability and is defined as
#' the expected probability of rejecting a null hypothesis when it was actually
#' true. \eqn{\alpha} is compared with the p-value and used as the decision
#' boundary for rejecting or not rejecting the null hypothesis.
#'
#' The family-wise error rate is the expected probability of making one or more
#' type I assertions among a family of hypotheses. Using Bonferroni's method,
#' \eqn{\alpha} is chosen for the family of hypotheses then divided by the
#' number of tests performed (\eqn{m}). Each individual hypothesis is tested at
#' \eqn{\frac{\alpha}{m}}. For example, if you plan to conduct 30 hypothesis
#' tests and want to control the family-wise error rate to no greater than
#' \eqn{\alpha = 0.05}, you would set `alpha = 0.05/30`.
#'
#' The per-family error rate is the expected number of type I assertions among a
#' family of hypotheses. If you calculate power for the scenario where you
#' perform 1,000 hypotheses and want to control the per-family error rate to no
#' greater than 10 type I assertions, you would choose `alpha = 10/1000`. This
#' implicitly assumes all 1,000 hypotheses are truly null. Alternatively, if you
#' assume 800 of these hypotheses are truly null and 200 are not,
#' `alpha = 10/1000` would control the per-family error rate to no greater than
#' 8 type I assertions. If it is acceptable to keep the per-family error rate as
#' 10, setting `alpha = 10/800` would provide greater marginal power than the
#' previous scenario.
#'
#' These two methods assume that the distribution of p-values for the truly null
#' hypotheses are uniform(0,1), but remain valid under various other testing
#' scenarios (such as dependent tests). Other multiple comparison methods, such
#' as FDR control, are common in practice but don't directly fit into this
#' power simulation framework.
#'
#' ## Column `nsims`
#' The final number of valid simulations per unique set of simulation parameters
#' may be less than the original number requested. This may occur when the test
#' results in a missing p-value. For `wald_test_bnb()`, pathological MLE
#' estimates, generally from small sample sizes and very small dispersions, may
#' result in a negative estimated standard deviation of the ratio. Thus the test
#' statistic and p-value would not be calculated. Note that simulated data from
#' `sim_nb()` and `sim_bnb()` may also reduce `nsims` during the data simulation
#' phase.
#'
#' The `nsims` column in the return data frame is the effective number of
#' simulations for power results.
#'
#' @references
#' \insertRef{yu_2017}{depower}
#'
#' \insertRef{yu_2020}{depower}
#'
#' \insertRef{rettiganti_2012}{depower}
#'
#' \insertRef{aban_2009}{depower}
#'
#' \insertRef{julious_2004}{depower}
#'
#' \insertRef{vickerstaff_2019}{depower}
#'
#' @param data (depower)\cr
#'        The simulated data returned by [depower::sim_log_lognormal()],
#'        [depower::sim_nb()], or [depower::sim_bnb()].
#' @param ... (functions)\cr
#'        The function(s) used to perform the test. If empty, a default test
#'        function will be selected based on `class(data)`. Names are used for
#'        labeling and should be unique. If names are empty, the call coerced to
#'        character will be used for name-value pairs. See 'Details'.
#' @param alpha (numeric: `0.05`; `(0,1)`)\cr
#'        The expected probability of rejecting a single null hypothesis when it
#'        is actually true. See 'Details'.
#' @param list_column (Scalar logical: `FALSE`)\cr
#'        If `TRUE`, the `data` and `result` list-columns are included in the
#'        returned data frame. If `FALSE` (default), the `data` and `result`
#'        list-columns are not included in the returned data frame.
#' @param ncores (Scalar integer: `1L`; `[1,Inf)`)\cr
#'        The number of cores (number of worker processes) to use. Do not set
#'        greater than the value returned by [parallel::detectCores()]. May be
#'        helpful when the number of parameter combinations is large and `nsims`
#'        is large.
#'
#' @return A data frame with the following columns appended to the `data` object:
#' \tabular{lll}{
#'   Column \tab Name \tab Description \cr
#'   \tab `alpha`  \tab Type I assertion probability. \cr
#'   \tab `test`   \tab Name-value pair of the function and statistical test: `c(as.character(...) = names(...)`. \cr
#'   \tab `data`   \tab List-column of simulated data. \cr
#'   \tab `result` \tab List-column of test results. \cr
#'   \tab `power`  \tab Power of the test for corresponding parameters.
#' }
#'
#' For `power(list_column = FALSE)`, columns `data`, and `result` are excluded from
#' the data frame.
#'
#' @seealso [depower::plot.depower()]
#'
#' @examples
#' #----------------------------------------------------------------------------
#' # power() examples
#' #----------------------------------------------------------------------------
#' library(depower)
#'
#' # Power for independent two-sample t-test
#' set.seed(1234)
#' data <- sim_log_lognormal(
#'   n1 = 20,
#'   n2 = 20,
#'   ratio = c(1.2, 1.4),
#'   cv1 = 0.4,
#'   cv2 = 0.4,
#'   cor = 0,
#'   nsims = 1000
#' )
#'
#' # Welch's t-test is used by default
#' power(data)
#'
#' # But you can specify anything else that is needed
#' power(
#'   data = data,
#'   "Welch's t-Test" = t_test_welch(alternative = "greater"),
#'   alpha = 0.01
#' )
#'
#' # Power for dependent two-sample t-test
#' set.seed(1234)
#' sim_log_lognormal(
#'   n1 = 20,
#'   n2 = 20,
#'   ratio = c(1.2, 1.4),
#'   cv1 = 0.4,
#'   cv2 = 0.4,
#'   cor = 0.5,
#'   nsims = 1000
#' ) |>
#'   power()
#'
#' # Power for mixed-type two-sample t-test
#' set.seed(1234)
#' sim_log_lognormal(
#'   n1 = 20,
#'   n2 = 20,
#'   ratio = c(1.2, 1.4),
#'   cv1 = 0.4,
#'   cv2 = 0.4,
#'   cor = c(0, 0.5),
#'   nsims = 1000
#' ) |>
#'   power()
#'
#' # Power for one-sample t-test
#' set.seed(1234)
#' sim_log_lognormal(
#'   n1 = 20,
#'   ratio = c(1.2, 1.4),
#'   cv1 = 0.4,
#'   nsims = 1000
#' ) |>
#'   power()
#'
#' # Power for independent two-sample NB test
#' set.seed(1234)
#' sim_nb(
#'   n1 = 10,
#'   mean1 = 10,
#'   ratio = c(1.6, 2),
#'   dispersion1 = 2,
#'   dispersion2 = 2,
#'   nsims = 200
#' ) |>
#'   power()
#'
#' # Power for BNB test
#' set.seed(1234)
#' sim_bnb(
#'   n = 10,
#'   mean1 = 10,
#'   ratio = c(1.4, 1.6),
#'   dispersion = 10,
#'   nsims = 200
#' ) |>
#'   power()
#'
#' @importFrom dplyr rowwise mutate collect cross_join ungroup
#' @importFrom multidplyr partition new_cluster cluster_library
#' @importFrom stats setNames
#'
#' @export
power <- function(
  data,
  ...,
  alpha = 0.05,
  list_column = FALSE,
  ncores = 1L
) {
  if (".funs" %in% ...names()) {
    stop(
      "Argument '.funs' is reserved for use in power.default(). Use power(...) instead."
    )
  }
  UseMethod("power")
}

#' @export
power.nb <- function(
  data,
  ...,
  alpha = 0.05,
  list_column = FALSE,
  ncores = 1L
) {
  #-----------------------------------------------------------------------------
  # Process dots
  #-----------------------------------------------------------------------------
  if (...length() == 0L) {
    .funs <- alist("NB Wald test" = wald_test_nb())
  } else {
    .funs <- dots(...)
  }

  #-----------------------------------------------------------------------------
  # Calculate power
  #-----------------------------------------------------------------------------
  power.default(
    data = data,
    .funs = .funs,
    alpha = alpha,
    list_column = list_column,
    ncores = ncores
  )
}

#' @export
power.bnb <- function(
  data,
  ...,
  alpha = 0.05,
  list_column = FALSE,
  ncores = 1L
) {
  #-----------------------------------------------------------------------------
  # Process dots
  #-----------------------------------------------------------------------------
  if (...length() == 0L) {
    .funs <- alist("BNB Wald test" = wald_test_bnb())
  } else {
    .funs <- dots(...)
  }

  #-----------------------------------------------------------------------------
  # Calculate power
  #-----------------------------------------------------------------------------
  power.default(
    data = data,
    .funs = .funs,
    alpha = alpha,
    list_column = list_column,
    ncores = ncores
  )
}

#' @export
power.log_lognormal_one_sample <- function(
  data,
  ...,
  alpha = 0.05,
  list_column = FALSE,
  ncores = 1L
) {
  #-----------------------------------------------------------------------------
  # Process dots
  #-----------------------------------------------------------------------------
  if (...length() == 0L) {
    .funs <- alist("One-sample t-Test" = t_test_paired())
  } else {
    .funs <- dots(...)
  }

  #-----------------------------------------------------------------------------
  # Calculate power
  #-----------------------------------------------------------------------------
  power.default(
    data = data,
    .funs = .funs,
    alpha = alpha,
    list_column = list_column,
    ncores = ncores
  )
}

#' @export
power.log_lognormal_independent_two_sample <- function(
  data,
  ...,
  alpha = 0.05,
  list_column = FALSE,
  ncores = 1L
) {
  #-----------------------------------------------------------------------------
  # Process dots
  #-----------------------------------------------------------------------------
  if (...length() == 0L) {
    .funs <- alist("Welch's t-Test" = t_test_welch())
  } else {
    .funs <- dots(...)
  }

  #-----------------------------------------------------------------------------
  # Calculate power
  #-----------------------------------------------------------------------------
  power.default(
    data = data,
    .funs = .funs,
    alpha = alpha,
    list_column = list_column,
    ncores = ncores
  )
}

#' @export
power.log_lognormal_dependent_two_sample <- function(
  data,
  ...,
  alpha = 0.05,
  list_column = FALSE,
  ncores = 1L
) {
  #-----------------------------------------------------------------------------
  # Process dots
  #-----------------------------------------------------------------------------
  if (...length() == 0L) {
    .funs <- alist("Paired t-Test" = t_test_paired())
  } else {
    .funs <- dots(...)
  }

  #-----------------------------------------------------------------------------
  # Calculate power
  #-----------------------------------------------------------------------------
  power.default(
    data = data,
    .funs = .funs,
    alpha = alpha,
    list_column = list_column,
    ncores = ncores
  )
}

#' @export
power.log_lognormal_mixed_two_sample <- function(
  data,
  ...,
  alpha = 0.05,
  list_column = FALSE,
  ncores = 1L
) {
  #-----------------------------------------------------------------------------
  # Process dots
  # If a function is not supplied, choose the default functions.
  # If an invalid function is supplied, error
  # If a valid function is supplied, pass that to the corresponding data subset.
  #-----------------------------------------------------------------------------
  if (...length() == 0L) {
    .funs <- alist(
      "Welch's t-Test" = t_test_welch(),
      "Paired t-Test" = t_test_paired()
    )
  } else {
    .funs <- dots(...)

    fun_deparse <- vapply(
      X = .funs,
      FUN = function(x) {
        deparse1(x[[1L]])
      },
      FUN.VALUE = character(1),
      USE.NAMES = FALSE
    )
    if (any(!fun_deparse %in% c("t_test_welch", "t_test_paired"))) {
      stop(
        "Argument '...' must only use t_test_welch() and t_test_paired() for mixed-type data. See ?depower::power()."
      )
    }
    if (!"t_test_welch" %in% fun_deparse) {
      .funs <- c(.funs, "Welch's t-Test" = quote(t_test_welch()))
    }
    if (!"t_test_paired" %in% fun_deparse) {
      .funs <- c(.funs, "Paired t-Test" = quote(t_test_paired()))
    }
  }

  fun_deparse <- vapply(
    X = .funs,
    FUN = function(x) {
      deparse1(x[[1L]])
    },
    FUN.VALUE = character(1),
    USE.NAMES = FALSE
  )
  welch_fun <- .funs[fun_deparse == "t_test_welch"]
  paired_fun <- .funs[fun_deparse == "t_test_paired"]

  #-----------------------------------------------------------------------------
  # Calculate power
  #-----------------------------------------------------------------------------
  res_welch <- power.default(
    data = data[data$cor == 0, ],
    .funs = welch_fun,
    alpha = alpha,
    list_column = list_column,
    ncores = ncores
  )

  res_paired <- power.default(
    data = data[data$cor != 0, ],
    .funs = paired_fun,
    alpha = alpha,
    list_column = list_column,
    ncores = ncores
  )

  rbind(res_welch, res_paired)
}

#' @inherit power
#' @param .funs (Named list of calls: `NULL`)\cr
#'        Internal use argument for `power.default()`. You should probably use
#'        `power(...)` instead. The function(s) used to perform the test. You
#'        may create the named list of calls using `alist()`. Names are used for
#'        labeling each test and should be unique.
#' @keywords Internal
#' @noRd
#' @export
power.default <- function(
  data,
  ...,
  alpha = 0.05,
  list_column = FALSE,
  ncores = 1L,
  .funs = NULL
) {
  #-----------------------------------------------------------------------------
  # Check arguments
  #-----------------------------------------------------------------------------
  if (!is.data.frame(data)) {
    stop("Argument 'data' must be a data frame.")
  }
  if (!inherits(x = data[["data"]], what = "list")) {
    stop(
      "Argument 'data' must contain a column named 'data' that is a list-column of simulated data."
    )
  }

  if (...length() > 0L) {
    stop(
      "Argument '...' must not be used in power.default(). Instead, use argument '.funs'."
    )
  }
  if (is.null(.funs)) {
    stop("Argument '.funs' must be a named list of calls.")
  }
  are_calls <- vapply(
    X = .funs,
    FUN = inherits,
    FUN.VALUE = logical(1L),
    what = "call",
    USE.NAMES = FALSE
  )
  # Should this refer to ... or .funs? ... is what the user sees.
  if (!all(are_calls)) {
    stop(
      "Argument '...' must be name-value pairs of calls. See ?depower::power for an explanation."
    )
  }
  fun_names <- names(.funs)

  if (!is.numeric(alpha) || any(alpha <= 0) || any(alpha >= 1)) {
    stop("Argument 'alpha' must be a numeric vector from (0,1).")
  }

  if (!is.numeric(ncores) || length(ncores) != 1L || ncores < 1L) {
    stop("Argument 'ncores' must be a positive scalar integer.")
  }
  if (ncores > 1L) {
    if (isTRUE(ncores > parallel::detectCores())) {
      max <- parallel::detectCores()
      warning("Argument 'ncores' should not be greater than ", max)
    }
  }

  #-----------------------------------------------------------------------------
  # Add .funs to data
  # Use list2DF as a hack to ensure name=value pairs stick in the data.frame.
  # Don't include alpha because we don't want to re-run the test for each alpha.
  #-----------------------------------------------------------------------------
  res <- dplyr::cross_join(
    x = data,
    y = list2DF(list(test = setNames(fun_names, as.character(.funs))))
  )

  #-----------------------------------------------------------------------------
  # Prepare for evaluation
  #-----------------------------------------------------------------------------
  # Ensure calls have correct data argument
  .funs <- lapply(.funs, function(x) {
    x[["data"]] <- quote(d)
    x
  })

  # A little more efficient than
  # dplyr::mutate(result = list(lapply(data, function(d) {eval(.funs[[test]])})))
  eval2 <- function(data, funs, fun_name) {
    f <- funs[[fun_name]]
    enclos <- parent.frame()
    lapply(
      X = data,
      FUN = function(d) {
        eval(expr = f, envir = list(d = d), enclos = enclos)
      }
    )
  }

  #-----------------------------------------------------------------------------
  # If randomization test, perform parametric simulation of the test statistic
  # distribution under the null hypothesis.
  #-----------------------------------------------------------------------------
  are_simulated <- vapply(
    X = .funs,
    FUN = function(x) {
      isTRUE(eval(x$distribution)$distribution == "simulated")
    },
    FUN.VALUE = logical(1L),
    USE.NAMES = TRUE
  )

  if (any(are_simulated)) {
    # Simulate data under the null hypothesis
    res_null <- sim_null(data, .funs)

    # Use asymptotic method to calculate observed test statistic.
    .funs_null <- lapply(
      .funs,
      function(x) {
        if (isTRUE(eval(x$distribution)$distribution == "simulated")) {
          x[["distribution"]] <- asymptotic()
          x
        } else {
          x
        }
      }
    )

    # Simulate distribution of the null test statistic
    if (ncores > 1L) {
      cluster <- multidplyr::new_cluster(ncores)
      multidplyr::cluster_library(cluster, 'depower')
    }
    res_null <- res_null |>
      dplyr::filter(.data$distribution_test_stat == "simulated") |>
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
        result = list(eval2(
          data = .data$data_null,
          funs = .funs_null,
          fun_name = .data$test
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

    get_test_stat <- function(result, test_stat) {
      vapply(
        X = result,
        FUN = function(x) {
          x[[test_stat]]
        },
        FUN.VALUE = numeric(1L)
      ) |>
        list()
    }

    res_null <- res_null |>
      dplyr::rowwise() |>
      dplyr::mutate(test_stat_null = get_test_stat(.data$result, "chisq")) |>
      dplyr::ungroup() |>
      dplyr::select(dplyr::any_of(names(res)), "test_stat_null")

    joinby <- intersect(names(res), names(res_null))
    res <- dplyr::left_join(res, res_null, by = joinby)

    .funs <- lapply(
      .funs,
      function(x) {
        if (isTRUE(eval(x$distribution)$distribution == "simulated")) {
          x[["distribution"]][["test_stat_null"]] <- quote(test_stat_null)
          x
        } else {
          x
        }
      }
    )
  }

  #-----------------------------------------------------------------------------
  # Run tests on the simulation data
  #-----------------------------------------------------------------------------
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
      result = list(eval2(data = data, funs = .funs, fun_name = .data$test))
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

  # Some tests may have a missing value for the p-value. We need to update the
  # total number of valid simulations.
  if (any_missing_p(res[["result"]])) {
    res <- res |>
      dplyr::rowwise() |>
      dplyr::mutate(nsims = n_valid_pvalue(.data$result)) |>
      dplyr::ungroup()
  }

  #-----------------------------------------------------------------------------
  # Add alpha
  #-----------------------------------------------------------------------------
  res <- dplyr::cross_join(
    x = res,
    y = data.frame(alpha = alpha)
  )

  #-----------------------------------------------------------------------------
  # Calculate power
  #-----------------------------------------------------------------------------
  res <- res |>
    dplyr::rowwise() |>
    dplyr::mutate(power = get_power(.data$result, alpha, .data$nsims)) |>
    dplyr::ungroup()

  #-----------------------------------------------------------------------------
  # Prepare return
  #-----------------------------------------------------------------------------
  # Subset
  if (!list_column) {
    res <- res |>
      dplyr::select(-dplyr::any_of(c("data", "result", "test_stat_null")))
  }

  # Column labels
  vars <- c(
    "n1" = "n1",
    "n2" = "n2",
    "n" = "n",
    "Mean1" = "mean1",
    "Mean2" = "mean2",
    "Mean" = "mean",
    "Ratio" = "ratio",
    "Difference" = "difference",
    "CV1" = "cv1",
    "CV2" = "cv2",
    "CV" = "cv",
    "Correlation" = "cor",
    "Dispersion1" = "dispersion1",
    "Dispersion2" = "dispersion2",
    "Dispersion" = "dispersion",
    "Distribution" = "distribution",
    "N Simulations" = "nsims",
    "Test" = "test",
    "Test statistic under the null" = "test_stat_null",
    "Alpha" = "alpha",
    "Data" = "data",
    "Result" = "result",
    "Power" = "power"
  )
  idx <- match(names(res), vars)
  if (anyNA(idx)) {
    stop("Unknown variable found while labeling data frame.")
  }
  for (i in seq_len(ncol(res))) {
    attr(res[[i]], "label") <- names(vars)[idx][i]
  }

  # Class
  res <- res[order(idx)] # subset will strip attributes, so put before
  class(res) <- class(data)

  #-----------------------------------------------------------------------------
  # Return
  #-----------------------------------------------------------------------------
  res
}

#===============================================================================
# Power helper functions
#===============================================================================
# Get power for a unique set of simulations
get_power <- function(result, alpha, nsims) {
  p <- vapply(
    X = result,
    FUN = function(x) {
      x[["p"]]
    },
    FUN.VALUE = numeric(1L)
  )

  sum(p <= alpha, na.rm = TRUE) / nsims
}

# Creates an index to subset list elements which are not majority zeros in the
# case when we are working in a rowwise dplyr pipeline
not_zeros <- function(x, max_zeros) {
  # Walk the list of data and return whether or not each dataset contains
  # majority nonzeros
  !vapply(
    X = x,
    FUN = function(x) {
      if (is.data.frame(x)) {
        x1 <- x[x[["condition"]] == 1, "value"]
        x2 <- x[x[["condition"]] == 2, "value"]
        xx1 <- sum(x1 == 0, na.rm = TRUE) / length(x1) > max_zeros
        xx2 <- sum(x2 == 0, na.rm = TRUE) / length(x2) > max_zeros
        xx1 || xx2
      } else {
        vapply(
          X = x,
          FUN = function(x) {
            sum(x == 0, na.rm = TRUE) / length(x) > max_zeros
          },
          FUN.VALUE = logical(1),
          USE.NAMES = FALSE
        ) |>
          any(na.rm = TRUE)
      }
    },
    FUN.VALUE = logical(1),
    USE.NAMES = FALSE
  )
}

# Whether or not the res[["data"]] object contains list elements which are
# majority zero.
any_zeros <- function(x, max_zeros) {
  # Walk the list of data and return whether or not each dataset contains
  # majority nonzeros
  lapply(
    X = x,
    FUN = function(x) {
      lapply(
        X = x,
        FUN = function(x) {
          if (is.data.frame(x)) {
            x1 <- x[x[["condition"]] == 1, "value"]
            x2 <- x[x[["condition"]] == 2, "value"]
            xx1 <- sum(x1 == 0, na.rm = TRUE) / length(x1) > max_zeros
            xx2 <- sum(x2 == 0, na.rm = TRUE) / length(x2) > max_zeros
            xx1 || xx2
          } else {
            vapply(
              X = x,
              FUN = function(x) {
                sum(x == 0, na.rm = TRUE) / length(x) > max_zeros
              },
              FUN.VALUE = logical(1),
              USE.NAMES = FALSE
            )
          }
        }
      )
    }
  ) |>
    unlist() |>
    any(na.rm = TRUE)
}

# Whether or not the res[["result"]] object contains missing p-values.
any_missing_p <- function(x) {
  # Walk the list of data and return whether or not each dataset contains
  # majority nonzeros
  lapply(
    X = x,
    FUN = function(x) {
      lapply(
        X = x,
        FUN = function(x) {
          is.na(x[["p"]])
        }
      )
    }
  ) |>
    unlist() |>
    any(na.rm = TRUE)
}

# Calculates the number of non-missing p-values in the case when we are working
# in a rowwise dplyr pipeline
n_valid_pvalue <- function(x) {
  # Walk the list of results and return the number of valid results
  vapply(
    X = x,
    FUN = function(x) {
      is.na(x[["p"]])
    },
    FUN.VALUE = logical(1),
    USE.NAMES = FALSE
  ) |>
    negate() |>
    sum()
}
