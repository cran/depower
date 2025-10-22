#' Test statistic distribution under the null
#'
#' Constructs a list which defines the test statistic reference distribution
#' under the null hypothesis.
#'
#' The default asymptotic test is performed for `distribution = asymptotic()`.
#'
#' When setting argument `distribution = simulated(method = "exact")`, the
#' exact randomization test is defined by:
#'
#' - Independent two-sample tests
#'     1. Calculate the observed test statistic.
#'     2. Check if `length(combn(x=n1+n2, m=n1))<1e6`
#'         1. If `TRUE` continue with the exact randomization test.
#'         2. If `FALSE` revert to the approximate randomization test.
#'     3. For all `combn(x=n1+n2, m=n1)` permutations:
#'         1. Assign corresponding group labels.
#'         2. Calculate the test statistic.
#'     4. Calculate the exact randomization test p-value as the mean of the logical
#'        vector `resampled_test_stats >= observed_test_stat`.
#' - Dependent two-sample tests
#'     1. Calculate the observed test statistic.
#'     2. Check if `npairs < 21` (maximum 2^20 resamples)
#'         1. If `TRUE` continue with the exact randomization test.
#'         2. If `FALSE` revert to the approximate randomization test.
#'     3. For all `2^npairs` permutations:
#'         1. Assign corresponding pair labels.
#'         2. Calculate the test statistic.
#'     4. Calculate the exact randomization test p-value as the mean of the logical
#'        vector `resampled_test_stats >= observed_test_stat`.
#'
#' For argument `distribution = simulated(method = "approximate")`, the
#' approximate randomization test is defined by:
#'
#' - Independent two-sample tests
#'     1. Calculate the observed test statistic.
#'     2. For `nsims` iterations:
#'         1. Randomly assign group labels.
#'         2. Calculate the test statistic.
#'     3. Insert the observed test statistic to the vector of resampled test
#'        statistics.
#'     4. Calculate the approximate randomization test p-value as the mean of
#'        the logical vector `resampled_test_stats >= observed_test_stat`.
#' - Dependent two-sample tests
#'     1. Calculate the observed test statistic.
#'     2. For `nsims` iterations:
#'         1. Randomly assign pair labels.
#'         2. Calculate the test statistic.
#'     3. Insert the observed test statistic to the vector of resampled test
#'        statistics.
#'     4. Calculate the approximate randomization test p-value as the mean of
#'        the logical vector `resampled_test_stats >= observed_test_stat`.
#'
#' In the power analysis setting, [depower::power()], we can simulate data for
#' groups 1 and 2 using their known distributions under the assumptions of the
#' null hypothesis. Unlike above where nonparametric randomization tests
#' are performed, in this setting approximate parametric tests are performed.
#'
#' For example, `power(wald_test_nb(distribution = simulated()))` would result
#' in an approximate parametric Wald test defined by:
#'
#' 1. For each relevant design row in `data`:
#'     1. For `simulated(nsims=integer())` iterations:
#'         1. Simulate new data for group 1 and group 2 under the null hypothesis.
#'         1. Calculate the Wald test statistic, \eqn{\chi^2_{null}}.
#'     1. Collect all \eqn{\chi^2_{null}} into a vector.
#'     1. For each of the `sim_nb(nsims=integer())` simulated datasets:
#'         1. Calculate the Wald test statistic, \eqn{\chi^2_{obs}}.
#'         1. Calculate the p-value based on the empirical null distribution of test statistics, \eqn{\chi^2_{null}}.
#'            (the mean of the logical vector `null_test_stats >= observed_test_stat`)
#'     1. Collect all p-values into a vector.
#'     1. Calculate power as `sum(p <= alpha) / nsims`.
#' 1. Return all results from `power()`.
#'
#' Randomization tests use the positive-biased p-value estimate in the style of
#' \insertCite{davison_1997;textual}{depower}
#' (see also \insertCite{phipson_2010;textual}{depower}):
#'
#' \deqn{
#' \hat{p} = \frac{1 + \sum_{i=1}^B \mathbb{I} \{\chi^2_i \geq \chi^2_{obs}\}}{B + 1}.
#' }
#'
#' The number of resamples defines the minimum observable p-value
#' (e.g. `nsims=1000L` results in min(p-value)=1/1001).
#' It's recommended to set \eqn{\text{nsims} \gg \frac{1}{\alpha}}.
#'
#' @references
#' \insertRef{davison_1997}{depower}
#'
#' \insertRef{phipson_2010}{depower}
#'
#' @param method (Scalar string: `"approximate"`)\cr
#'        The method used to derive the distribution of the test statistic
#'        under the null hypothesis. Must be one of `"approximate"` (default) or
#'        `"exact"`. See 'Details' for additional information.
#' @param nsims (Scalar integer: `1000L`; `[2, Inf)`)\cr
#'        The number of resamples for `method = "approximate"`. Not used for
#'        `method = "exact"`, except for the case when the number of exact
#'        resamples exceeds approximately `1e6` and then `method = "approximate"`
#'        will be used as a fallback. In the [depower::power()] context, `nsims`
#'        defines the number of simulated datasets under the null hypothesis.
#'        For this case you would typically set `nsims` as greater than or equal
#'        to the number of simulated datasets in the design row of the power analysis.
#'        See 'Details' for additional information.
#' @param ncores (Scalar integer: `1L`; `[1, Inf)`)\cr
#'        The number of cores (number of worker processes) to use. Do not set
#'        greater than the value returned by [parallel::detectCores()].
#' @param ... Optional arguments for internal use.
#'
#' @return list
#'
#' @examples
#' #----------------------------------------------------------------------------
#' # simulated() examples
#' #----------------------------------------------------------------------------
#' data |>
#'   wald_test_nb(distribution = simulated(nsims = 200L))
#'
#' @name distribution
NULL

#' @export
#' @rdname distribution
simulated <- function(
  method = "approximate",
  nsims = 1000L,
  ncores = 1L,
  ...
) {
  #-----------------------------------------------------------------------------
  # Check arguments
  #-----------------------------------------------------------------------------
  if (!method %in% c("approximate", "exact")) {
    stop("Argument 'method' must be one of 'approximate' or 'exact'.")
  }
  if (!(is.numeric(nsims) && length(nsims) == 1L && nsims > 1L)) {
    stop("Argument 'nsims' must be a scalar integer greater than 1.")
  }
  if (!(is.numeric(ncores) && length(ncores) == 1L && ncores > 0L)) {
    stop("Argument 'ncores' must be a positive scalar integer.")
  }

  #-----------------------------------------------------------------------------
  # Return
  #-----------------------------------------------------------------------------
  list(
    distribution = "simulated",
    method = method,
    nsims = nsims,
    ncores = ncores,
    ...
  )
}

check_simulated <- function(x) {
  nms <- names(x)
  if (
    !is.list(x) ||
      !all(c("distribution", "method", "nsims", "ncores") %in% nms)
  ) {
    stop(
      "Check argument 'distribution'. You must use 'distribution = depower::simulated()'."
    )
  }
  if (x$distribution != "simulated") {
    stop(
      "Check argument 'distribution'. You must use 'distribution = depower::simulated()'."
    )
  }
  if (!x$method %in% c("approximate", "exact")) {
    stop(
      "Argument 'method' in 'depower::simulated()' must be one of 'approximate' or 'exact'."
    )
  }
  if (!(is.numeric(x$nsims) && length(x$nsims) == 1L && x$nsims > 1L)) {
    stop(
      "Argument 'nsims' in 'depower::simulated()' must be a scalar integer greater than 1."
    )
  }
  if (!(is.numeric(x$ncores) && length(x$ncores) == 1L && x$ncores > 0L)) {
    stop(
      "Argument 'ncores' in 'depower::simulated()' must be a positive scalar integer."
    )
  }

  invisible(x)
}
