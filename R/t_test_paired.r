#' Paired and one-sample t-Tests
#'
#' Performs paired and one-sample t-Tests.
#'
#' @details
#' This function is primarily designed for speed in simulation. Missing values
#' are silently excluded.
#'
#' The one-sample test is used for both the true one-sample scenario and for the
#' paired differences from a dependent two-sample scenario. Below we use paired
#' difference language as that is the most common case. The hypotheses for the
#' paired t-test are
#'
#' \deqn{
#' \begin{aligned}
#' H_{null} &: \mu_{diff} = \mu_{null} \\
#' H_{alt} &: \begin{cases}
#'   \mu_{diff} \neq \mu_{null} & \text{two-sided}\\
#'   \mu_{diff} > \mu_{null} & \text{greater than}\\
#'   \mu_{diff} < \mu_{null} & \text{less than}
#' \end{cases}
#' \end{aligned}
#' }
#'
#' where \eqn{\mu_{diff} = AM(X_2 - X_1)} is the arithmetic mean of the paired
#' differences (sample 2 - sample 1) and \eqn{\mu_{null}} is a constant for the
#' assumed population mean difference (usually \eqn{\mu_{null} = 0}).
#'
#' The test statistic is
#'
#' \deqn{
#' T = \frac{\bar{x}_{diff} - \mu_{null}}{\sqrt{\frac{s^2}{n}}}
#' }
#'
#' where \eqn{\bar{x}_{diff}} is the sample mean of the differences, \eqn{\mu_{null}}
#' is the population mean difference assumed under the null hypothesis, \eqn{n}
#' is the sample size of the differences, and \eqn{s^2} is the sample variance.
#'
#' The critical value of the test statistic has degrees of freedom
#'
#' \deqn{
#' df = n-1
#' }
#'
#' and the p-value is calculated as
#'
#' \deqn{
#' \begin{aligned}
#' p &= \begin{cases}
#'   2 \text{min} \{P(T \geq t_{n-1} \mid H_{null}), P(T \leq t_{n-1} \mid H_{null})\} & \text{two-sided}\\
#'   P(T \geq t_{n-1} \mid H_{null}) & \text{greater than}\\
#'   P(T \leq t_{n-1} \mid H_{null}) & \text{less than}
#' \end{cases}
#' \end{aligned}
#' }
#'
#' Let \eqn{GM(\cdot)} be the geometric mean and \eqn{AM(\cdot)} be the
#' arithmetic mean. For dependent lognormal samples \eqn{X_1} and \eqn{X_2} it
#' follows that \eqn{\ln X_1} and \eqn{\ln X_2} are dependent normally
#' distributed variables. Setting \eqn{\mu_{diff} = AM(\ln X_2 - \ln X_1)}
#' we have
#'
#' \deqn{
#' e^{\mu_{diff}} = GM\left( \frac{X_2}{X_1} \right)
#' }
#'
#' This forms the basis for making inference about the geometric mean ratio of
#' the original lognormal data using the mean difference of the log transformed
#' normal data.
#'
#' @references
#' \insertRef{julious_2004}{depower}
#'
#' \insertRef{hauschke_1992}{depower}
#'
#' \insertRef{johnson_1994}{depower}
#'
#' @param data (list)\cr
#'        A list whose first element is the vector of normal values from sample
#'        1 and the second element is the vector of normal values from sample 2.
#'        Both vectors must be the same sample size and sorted by the
#'        subject/item index. If `length(data) == 1L`, the one-sample test is
#'        used. \link[base]{NA}s are silently excluded. The default output
#'        from [depower::sim_log_lognormal()].
#' @param alternative (string: `"two.sided"`)\cr
#'        The alternative hypothesis. Must be one of `"two.sided"`, `"greater"`,
#'        or `"less"`. See 'Details' for additional information.
#' @param ci_level (Scalar numeric: `NULL`; `(0, 1)`)\cr
#'        If `NULL`, confidence intervals are set as `NA`. If in `(0, 1)`,
#'        confidence intervals are calculated at the specified level.
#' @param mean_null (Scalar numeric: `0`; `(-Inf, Inf)`)\cr
#'        The mean or mean difference assumed under the null hypothesis. See
#'        'Details' for additional information.
#'
#' @return
#' A list with the following elements:
#' \tabular{llll}{
#'   Slot \tab Subslot \tab Name \tab Description \cr
#'   1 \tab   \tab `t`  \tab Value of the t-statistic. \cr
#'   2 \tab   \tab `df` \tab Degrees of freedom for the t-statistic. \cr
#'   3 \tab   \tab `p`  \tab p-value. \cr
#'
#'   4 \tab   \tab `mean_diff` \tab Estimated mean or mean of the differences
#'                             (sample 2 – sample 1). \cr
#'   4 \tab 1 \tab `estimate`  \tab Point estimate. \cr
#'   4 \tab 2 \tab `lower`     \tab Confidence interval lower bound. \cr
#'   4 \tab 3 \tab `upper`     \tab Confidence interval upper bound. \cr
#'
#'   5 \tab   \tab `n`           \tab Number of paired observations. \cr
#'   6 \tab   \tab `method`      \tab Method used for the results. \cr
#'   7 \tab   \tab `alternative` \tab The alternative hypothesis. \cr
#'   8 \tab   \tab `ci_level`    \tab The confidence level. \cr
#'   9 \tab   \tab `mean_null`   \tab Assumed population mean of the differences
#'                               under the null hypothesis.
#' }
#'
#' @seealso [depower::t_test_welch()]
#'
#' @examples
#' #----------------------------------------------------------------------------
#' # t_test_paired() examples
#' #----------------------------------------------------------------------------
#' library(depower)
#'
#' # One-sample t-test
#' set.seed(1234)
#' t_test1 <- sim_log_lognormal(
#'   n1 = 40,
#'   ratio = 1.5,
#'   cv1 = 0.4
#' ) |>
#'   t_test_paired(ci_level = 0.95)
#'
#' t_test1
#'
#' # Paired t-test using two dependent samples
#' set.seed(1234)
#' t_test2 <- sim_log_lognormal(
#'   n1 = 40,
#'   n2 = 40,
#'   ratio = 1.5,
#'   cv1 = 0.4,
#'   cv2 = 0.2,
#'   cor = 0.3
#' ) |>
#'   t_test_paired(ci_level = 0.95)
#'
#' t_test2
#'
#' @importFrom stats pt qt complete.cases
#'
#' @export
t_test_paired <- function(
  data,
  alternative = "two.sided",
  ci_level = NULL,
  mean_null = 0
) {
  #-----------------------------------------------------------------------------
  # Check arguments
  #-----------------------------------------------------------------------------
  len <- length(data)
  if (!is.list(data) || len == 0L || len > 2L) {
    stop("Argument 'data' must be a list of length 1 or 2.")
  }

  has_ci <- !is.null(ci_level)
  if (has_ci) {
    if (!(is.numeric(ci_level) && length(ci_level) == 1L)) {
      stop("Argument 'ci_level' must be a scalar numeric.")
    }
  } else {
    mean_diff_lower <- NA_real_
    mean_diff_upper <- NA_real_
    ci_level <- NA_real_
  }

  if (!(is.numeric(mean_null) && length(mean_null) == 1L)) {
    stop("Argument 'mean_null' must be a scalar numeric.")
  }

  #-----------------------------------------------------------------------------
  # t-test
  #-----------------------------------------------------------------------------
  if (len == 2L) {
    d1 <- data[[1L]]
    d2 <- data[[2L]]
    if (length(d1) != length(d2)) {
      stop("Argument 'data' must have the same sample size for both samples.")
    }
    if (anyNA(d1) || anyNA(d2)) {
      not_na <- complete.cases(d1, d2)
      diff <- d2[not_na] - d1[not_na]
    } else {
      diff <- d2 - d1
    }
  } else {
    diff <- data[[1L]]
    if (anyNA(diff)) {
      diff <- diff[!is.na(diff)]
    }
  }

  n <- length(diff)
  mean_diff <- fmean(x = diff, n = n)
  var_diff <- fvar(x = diff, n = n, mean = mean_diff)
  df <- n - 1L
  stderr <- sqrt(var_diff / n)
  tstat <- (mean_diff - mean_null) / stderr

  if (alternative == "two.sided") {
    # A bit faster than pt(abs(tstat), df, lower.tail = FALSE)?
    p <- 2 * pt(-abs(tstat), df)
    method <- "One-sample t-test for 'two-sided' alternative"
    if (has_ci) {
      alpha <- 1 - (1 - ci_level) / 2
      crit <- qt(p = alpha, df = df)
      mean_diff_lower <- mean_null + (tstat - crit) * stderr
      mean_diff_upper <- mean_null + (tstat + crit) * stderr
    }
  } else if (alternative == "greater") {
    # A bit faster than pt(tstat, df, lower.tail = FALSE)?
    p <- pt(-tstat, df)
    method <- "One-sample t-test for 'greater than' alternative"
    if (has_ci) {
      crit <- qt(p = ci_level, df = df)
      mean_diff_lower <- mean_null + (tstat - crit) * stderr
      mean_diff_upper <- Inf
    }
  } else if (alternative == "less") {
    p <- pt(tstat, df)
    method <- "One-sample t-test for 'less than' alternative"
    if (has_ci) {
      crit <- qt(p = ci_level, df = df)
      mean_diff_lower <- -Inf
      mean_diff_upper <- mean_null + (tstat + crit) * stderr
    }
  } else {
    stop(
      "Argument 'alternative' must be one of 'two.sided', 'greater' or 'less'."
    )
  }

  #-----------------------------------------------------------------------------
  # Return
  #-----------------------------------------------------------------------------
  list(
    t = tstat,
    df = df,
    p = p,
    mean_diff = list(
      estimate = mean_diff,
      lower = mean_diff_lower,
      upper = mean_diff_upper
    ),
    n = n,
    method = method,
    alternative = alternative,
    ci_level = ci_level,
    mean_null = mean_null
  )
}
