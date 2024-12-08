#' Negative log-likelihood for BNB
#'
#' The negative log-likelihood for bivariate negative binomial outcomes.
#'
#' @details
#' These functions are primarily designed for speed in simulation. Arguments are
#' not checked.
#'
#' Suppose \eqn{X_1 \mid G = g \sim \text{Poisson}(\mu g)} and
#' \eqn{X_2 \mid G = g \sim \text{Poisson}(r \mu g)} where
#' \eqn{G \sim \text{Gamma}(\theta, \theta^{-1})} is the random item (subject) effect.
#' Then \eqn{X_1, X_2 \sim \text{BNB}(\mu, r, \theta)} is the joint distribution where
#' \eqn{X_1} and \eqn{X_2} are dependent (though conditionally independent),
#' \eqn{X_1} is the count outcome for sample 1 of the items (subjects),
#' \eqn{X_2} is the count outcome for sample 2 of the items (subjects),
#' \eqn{\mu} is the conditional mean of sample 1, \eqn{r} is the ratio of the
#' conditional means of sample 2 with respect to sample 1, and \eqn{\theta} is
#' the gamma distribution shape parameter which controls the dispersion and the
#' correlation between sample 1 and 2.
#'
#' The likelihood is
#'
#' \deqn{
#' \begin{aligned}
#' L(r, \mu, \theta \mid X_1, X_2) = & \left( \frac{\theta^{\theta}}{\Gamma(\theta)} \right)^{n} \times \\
#'   & \frac{\mu^{\sum{x_{1i}} + \sum{x_{2i}}}}{\prod_{i=1}^{n} x_{1i}!} \frac{r^{\sum{x_{2i}}}}{\prod_{i=1}^{n} x_{2i}!} \times \\
#'   & \frac{\prod_{i = 1}^{n} \Gamma(x_{1i} + x_{2i} + \theta)}{(\mu + r \mu + \theta)^{\sum x_{1i} + x_{2i} + \theta}}
#' \end{aligned}
#' }
#'
#' and the parameter space is
#' \eqn{\Theta = \left\{ (r, \mu, \theta) : r, \mu, \theta > 0 \right\}}.
#' The log-likelihood is
#'
#' \deqn{
#' \begin{aligned}
#' l(r, \mu, \theta) = \ & n \left[ \theta \ln \theta - \ln \Gamma(\theta) \right] + \\
#'   & n (\bar{x}_1 + \bar{x}_2) \ln(\mu) + n \bar{x}_2 \ln r + \\
#'   & \sum_{i=1}^{n}{\ln \Gamma(x_{1i} + x_{2i} + \theta)} - \\
#'   & n (\bar{x}_1 + \bar{x}_2 + \theta) \ln(\mu + r\mu + \theta) - \\
#'   & \sum_{i = 1}^{n}{\ln x_{1i}!} - \sum_{i = 1}^{n}{\ln x_{2i}!}
#' \end{aligned}
#' }
#'
#' @references
#' \insertRef{rettiganti_2012}{depower}
#'
#' \insertRef{aban_2009}{depower}
#'
#' @param param (numeric)\cr
#'        The vector of initial values for BNB parameters. Must be in the
#'        following order for each scenario:
#' - Null: `c(mean1, dispersion)`
#' - Alternative: `c(mean1, mean2, dispersion)`
#'
#' for samples 1 and 2.
#' @param value1 (numeric)\cr
#'        The vector of BNB values from sample 1. Must be sorted by the
#'        subject/item index. Must not contain \link[base]{NA}s.
#' @param value2 (numeric)\cr
#'        The vector of BNB values from sample 2. Must be sorted by the
#'        subject/item index. Must not contain \link[base]{NA}s.
#' @param ratio_null (Scalar numeric: `(0, Inf)`)\cr
#'        The ratio of means assumed under the null hypothesis (sample 2 / sample 1).
#'        Typically `ratio_null = 1` (no difference).
#'
#' @return Scalar numeric negative log-likelihood.
#'
#' @seealso [depower::sim_nb()], [stats::nlminb()], [stats::nlm()],
#'          [stats::optim()]
#'
#' @examples
#' #----------------------------------------------------------------------------
#' # nll_bnb*() examples
#' #----------------------------------------------------------------------------
#' library(depower)
#'
#' set.seed(1234)
#' d <- sim_bnb(
#'   n = 40,
#'   mean1 = 10,
#'   ratio = 1.2,
#'   dispersion = 2
#' )
#'
#' nll_bnb_alt(
#'   param = c(mean1 = 10, mean2 = 12, dispersion = 2),
#'   value1 = d[[1L]],
#'   value2 = d[[2L]]
#' )
#'
#' nll_bnb_null(
#'   param = c(mean = 10, dispersion = 2),
#'   value1 = d[[1L]],
#'   value2 = d[[2L]],
#'   ratio_null = 1
#' )
#'
#' @name nll_bnb
NULL

#' @export
#' @rdname nll_bnb
nll_bnb_null <- function(param, value1, value2, ratio_null) {
  mu1 <- param[1L]
  dispersion <- param[2L]
  n <- length(value1)
  mean1 <- fmean(value1, n = n)
  mean2 <- fmean(value2, n = n)
  sum1 <- sum(lgamma(value1 + value2 + dispersion))
  sum2 <- sum(lgamma(value1 + 1L))
  sum3 <- sum(lgamma(value2 + 1L))

  ll <- n * (dispersion * log(dispersion) - lgamma(dispersion)) +
        n * (mean1 + mean2) * log(mu1) +
        n * mean2 * log(ratio_null) +
        sum1 -
        n * (mean1 + mean2 + dispersion) * log(mu1 + ratio_null * mu1 + dispersion) -
        sum2 -
        sum3

  # Remove names and return negative log-likelihood
  as.numeric(-ll)
}

#' @export
#' @rdname nll_bnb
nll_bnb_alt <- function(param, value1, value2) {
  mu1 <- param[1L]
  mu2 <- param[2L]
  dispersion <- param[3L]
  n <- length(value1)
  mean1 <- fmean(x = value1, n = n)
  mean2 <- fmean(x = value2, n = n)
  sum1 <- sum(lgamma(value1 + value2 + dispersion))
  sum2 <- sum(lgamma(value1 + 1L))
  sum3 <- sum(lgamma(value2 + 1L))

  ll <- n * (dispersion * log(dispersion) - lgamma(dispersion)) +
        n * mean1 * log(mu1) +
        n * mean2 * log(mu2) +
        sum1 -
        n * (mean1 + mean2 + dispersion) * log(mu1 + mu2 + dispersion) -
        sum2 -
        sum3

  # Remove names and return negative log-likelihood
  as.numeric(-ll)
}
