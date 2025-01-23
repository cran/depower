#' Welch's t-Test
#'
#' Performs Welch's independent two-sample t-test.
#'
#' @details
#' This function is primarily designed for speed in simulation. Missing values
#' are silently excluded.
#'
#' The hypotheses for Welch's independent two-sample t-test are
#'
#' \deqn{
#' \begin{aligned}
#' H_{null} &: \mu_2 - \mu_1 = \mu_{null} \\
#' H_{alt} &: \begin{cases}
#'   \mu_2 - \mu_1 \neq \mu_{null} & \text{two-sided}\\
#'   \mu_2 - \mu_1 > \mu_{null} & \text{greater than}\\
#'   \mu_2 - \mu_1 < \mu_{null} & \text{less than}
#' \end{cases}
#' \end{aligned}
#' }
#'
#' where \eqn{\mu_1} is the population mean of group 1, \eqn{\mu_2} is the
#' population mean of group 2, and \eqn{\mu_{null}} is a constant for the assumed
#' difference of population means (usually \eqn{\mu_{null} = 0}).
#'
#' The test statistic is
#'
#' \deqn{
#' T = \frac{(\bar{x}_2 - \bar{x}_1) - \mu_{null}}{\sqrt{\frac{s_1^2}{n_1} + \frac{s_2^2}{n_2}}}
#' }
#'
#' where \eqn{\bar{x}_1} and \eqn{\bar{x}_2} are the sample means, \eqn{\mu_{null}}
#' is the difference of population means assumed under the null hypothesis,
#' \eqn{n_1} and \eqn{n_2} are the sample sizes, and \eqn{s_1^2} and \eqn{s_2^2}
#' are the sample variances.
#'
#' The critical value of the test statistic uses the Welch–Satterthwaite degrees
#' of freedom
#'
#' \deqn{
#' v = \frac{\left( \frac{s_1^2}{n_1} + \frac{s_2^2}{n_2} \right)^2}
#'          {(N_1 - 1)^{-1}\left( \frac{s_1^2}{n_1} \right)^2 +
#'          (N_2 - 1)^{-1}\left( \frac{s_2^2}{n_2} \right)^2}
#' }
#'
#' and the p-value is calculated as
#'
#' \deqn{
#' \begin{aligned}
#' p &= \begin{cases}
#'   2 \text{min} \{P(T \geq t_{v} \mid H_{null}), P(T \leq t_{v} \mid H_{null})\} & \text{two-sided}\\
#'   P(T \geq t_{v} \mid H_{null}) & \text{greater than}\\
#'   P(T \leq t_{v} \mid H_{null}) & \text{less than}
#' \end{cases}
#' \end{aligned}
#' }
#'
#' Let \eqn{GM(\cdot)} be the geometric mean and \eqn{AM(\cdot)} be the
#' arithmetic mean. For independent lognormal variables \eqn{X_1} and \eqn{X_2}
#' it follows that \eqn{\ln X_1} and \eqn{\ln X_2} are independent normally
#' distributed variables. Defining
#' \eqn{\mu_{X_2} - \mu_{X_1} = AM(\ln X_2) - AM(\ln X_1)}
#' we have
#'
#' \deqn{
#' e^{\mu_{X_2} - \mu_{X_1}} = \frac{GM(X_2)}{GM(X_1)}
#' }
#'
#' This forms the basis for making inference about the ratio of geometric means
#' of the original lognormal data using the difference of means of the log
#' transformed normal data.
#'
#' @references
#' \insertRef{julious_2004}{depower}
#'
#' \insertRef{hauschke_1992}{depower}
#'
#' \insertRef{johnson_1994}{depower}
#'
#' @param data (list)\cr
#'        A list whose first element is the vector of normal values from group
#'        1 and the second element is the vector of normal values from group 2.
#'        \link[base]{NA}s are silently excluded. The default output from
#'        [depower::sim_log_lognormal()].
#' @param alternative (string: `"two.sided"`)\cr
#'        The alternative hypothesis. Must be one of `"two.sided"`, `"greater"`,
#'        or `"less"`. See 'Details' for additional information.
#' @param ci_level (Scalar numeric: `NULL`; `(0, 1)`)\cr
#'        If `NULL`, confidence intervals are set as `NA`. If in `(0, 1)`,
#'        confidence intervals are calculated at the specified level.
#' @param mean_null (Scalar numeric: `0`; `(-Inf, Inf)`)\cr
#'        The difference of means assumed under the null hypothesis. See
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
#'   4 \tab   \tab `diff_mean` \tab Estimated difference of means
#'                             (group 2 – group 1). \cr
#'   4 \tab 1 \tab `estimate`  \tab Point estimate. \cr
#'   4 \tab 2 \tab `lower`     \tab Confidence interval lower bound. \cr
#'   4 \tab 3 \tab `upper`     \tab Confidence interval upper bound. \cr
#'
#'   5 \tab   \tab `mean1` \tab Estimated mean of group 1. \cr
#'   6 \tab   \tab `mean2` \tab Estimated mean of group 2. \cr
#'   7 \tab   \tab `n1` \tab Sample size of group 1. \cr
#'   8 \tab   \tab `n2` \tab Sample size of group 2. \cr
#'   9  \tab   \tab `method`      \tab Method used for the results. \cr
#'   10 \tab   \tab `alternative` \tab The alternative hypothesis. \cr
#'   11 \tab   \tab `ci_level`    \tab The confidence level. \cr
#'   12 \tab   \tab `mean_null`   \tab Assumed population difference of the means
#'                                under the null hypothesis.
#' }
#'
#' @seealso [stats::t.test()]
#'
#' @examples
#' #----------------------------------------------------------------------------
#' # t_test_welch() examples
#' #----------------------------------------------------------------------------
#' library(depower)
#'
#' # Welch's t-test
#' set.seed(1234)
#' sim_log_lognormal(
#'   n1 = 40,
#'   n2 = 40,
#'   ratio = 1.5,
#'   cv1 = 0.4,
#'   cv2 = 0.4
#' ) |>
#'   t_test_welch(ci_level = 0.95)
#'
#' @importFrom stats pt qt
#'
#' @export
t_test_welch <- function(
    data,
    alternative = "two.sided",
    ci_level = NULL,
    mean_null = 0
) {
  #-----------------------------------------------------------------------------
  # Check arguments
  #-----------------------------------------------------------------------------
  if(!(is.list(data) && length(data) == 2L)) {
    stop("Argument 'data' must be a list of length 2.")
  }

  has_ci <- !is.null(ci_level)
  if(has_ci) {
    if(!(is.numeric(ci_level) && length(ci_level == 1L))) {
      stop("Argument 'ci_level' must be a scalar numeric.")
    }
  } else {
    diff_mean_lower <- NA_real_
    diff_mean_upper <- NA_real_
    ci_level <- NA_real_
  }

  if(!(is.numeric(mean_null) && length(mean_null) == 1L)) {
    stop("Argument 'mean_null' must be a scalar numeric.")
  }

  #-----------------------------------------------------------------------------
  # t-test
  #-----------------------------------------------------------------------------
  d1 <- data[[1L]]
  if(anyNA(d1)) {
    d1 <- d1[!is.na(d1)]
  }
  d2 <- data[[2L]]
  if(anyNA(d2)) {
    d2 <- d2[!is.na(d2)]
  }

  n1 <- length(d1)
  n2 <- length(d2)
  mean1 <- fmean(x = d1, n = n1)
  mean2 <- fmean(x = d2, n = n2)
  diff_mean <- mean2 - mean1
  var_n1 <- fvar(x = d1, n = n1, mean = mean1) / n1
  var_n2 <- fvar(x = d2, n = n2, mean = mean2) / n2
  sum_var_n <- var_n1 + var_n2
  stderr <- sqrt(sum_var_n)

  tstat <- (diff_mean - mean_null) / stderr
  df <- sum_var_n^2 / (var_n1^2 / (n1 - 1L) + var_n2^2 / (n2 - 1L))

  if(alternative == "two.sided") {
    # A bit faster than pt(abs(tstat), df, lower.tail = FALSE)?
    p <- 2 * pt(-abs(tstat), df)
    method <- "Welch's two-sample t-test for 'two-sided' alternative"
    if(has_ci) {
      alpha <- 1 - (1 - ci_level) / 2
      crit <- qt(p = alpha, df = df)
      diff_mean_lower <- mean_null + (tstat - crit) * stderr
      diff_mean_upper <- mean_null + (tstat + crit) * stderr
    }
  } else if(alternative == "greater") {
    # A bit faster than pt(tstat, df, lower.tail = FALSE)?
    p <- pt(-tstat, df)
    method <- "Welch's two-sample t-test for 'greater than' alternative"
    if(has_ci) {
      crit <- qt(p = ci_level, df = df)
      diff_mean_lower <- mean_null + (tstat - crit) * stderr
      diff_mean_upper <- Inf
    }
  } else if(alternative == "less") {
    p <- pt(tstat, df)
    method <- "Welch's two-sample t-test for 'less than' alternative"
    if(has_ci) {
      crit <- qt(p = ci_level, df = df)
      diff_mean_lower <- -Inf
      diff_mean_upper <- mean_null + (tstat + crit) * stderr
    }
  } else {
    stop("Argument 'alternative' must be one of 'two.sided', 'greater' or 'less'.")
  }

  #-----------------------------------------------------------------------------
  # Return
  #-----------------------------------------------------------------------------
  list(
    t = tstat,
    df = df,
    p = p,
    diff_mean = list(
      estimate = diff_mean,
      lower = diff_mean_lower,
      upper = diff_mean_upper
    ),
    mean1 = mean1,
    mean2 = mean2,
    n1 = n1,
    n2 = n2,
    method = method,
    alternative = alternative,
    ci_level = ci_level,
    mean_null = mean_null
  )
}
