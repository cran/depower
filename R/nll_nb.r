#' Negative log-likelihood for NB
#'
#' The negative log-likelihood for two independent samples of negative binomial
#' distributions.
#'
#' @details
#' These functions are primarily designed for speed in simulation. Limited
#' argument validation is performed.
#'
#' Suppose \eqn{X_1 \sim \text{NB}(\mu, \theta_1)} and
#' \eqn{X_2 \sim \text{NB}(r\mu, \theta_2)} where \eqn{X_1} and \eqn{X_2} are
#' independent, \eqn{X_1} is the count outcome for items in group 1, \eqn{X_2}
#' is the count outcome for items in group 2, \eqn{\mu} is the arithmetic mean
#' count in group 1, \eqn{r} is the ratio of arithmetic means for group 2 with
#' respect to group 1, \eqn{\theta_1} is the dispersion parameter of group 1,
#' and \eqn{\theta_2} is the dispersion parameter of group 2.
#'
#' ## Unequal dispersion parameters
#'
#' When the dispersion parameters are not equal, the likelihood is
#'
#' \deqn{
#' \begin{aligned}
#' L(r, \mu, \theta_1, \theta_2 \mid X_1, X_2) = & \left( \frac{\theta_1^{\theta_1}}{\Gamma(\theta_1)} \right)^{n_1} \frac{\mu^{\sum{x_{1i}}}}{(\mu + \theta_1)^{\sum{x_{1i} + n_1 \theta_1}}} \times \\
#'   & \left( \frac{\theta_2^{\theta_2}}{\Gamma(\theta_2)} \right)^{n_2} \frac{(r \mu)^{\sum{x_{2j}}}}{(r \mu + \theta_2)^{\sum{x_{2j} + n_2 \theta_2}}} \times \\
#'   & \prod_{i = 1}^{n_1}{\frac{\Gamma(x_{1i} + \theta_1)}{x_{1i}!}} \prod_{j = 1}^{n_2}{\frac{\Gamma(x_{2j} + \theta_2)}{x_{2j}!}}
#' \end{aligned}
#' }
#'
#' and the parameter space is
#' \eqn{\Theta = \left\{ (r, \mu, \theta_1, \theta_2) : r, \mu, \theta_1, \theta_2 > 0 \right\}}.
#' The log-likelihood is
#'
#' \deqn{
#' \begin{aligned}
#' l(r, \mu, \theta_1, \theta_2) = \ &n_1 \left[ \theta_1 \ln \theta_1 - \ln \Gamma(\theta_1) \right] + \\
#'   &n_2 \left[ \theta_2 \ln \theta_2 - \ln \Gamma(\theta_2) \right] + \\
#'   &(n_1 \bar{x}_1 + n_2 \bar{x}_2) \ln(\mu) - n_1 (\bar{x}_1 + \theta_1) \ln(\mu + \theta_1) + \\
#'   &n_2 \bar{x}_2 \ln(r) - n_2 (\bar{x}_2 + \theta_2) \ln(r \mu + \theta_2) + \\
#'   &\sum_{i = 1}^{n_1}{\left( \ln \Gamma(x_{1i} + \theta_1) - \ln(x_{1i}!) \right)} + \\
#'   &\sum_{j = 1}^{n_2}{\left( \ln \Gamma(x_{2j} + \theta_2) - \ln(x_{2j}!) \right)}
#' \end{aligned}
#' }
#'
#' ## Equal dispersion parameters
#'
#' When the dispersion parameters are equal, the likelihood is
#'
#' \deqn{
#' \begin{aligned}
#' L(r, \mu, \theta \mid X_1, X_2) = & \left( \frac{\theta^{\theta}}{\Gamma(\theta)} \right)^{n_1 + n_2} \times \\
#'   & \frac{\mu^{\sum{x_{1i}}}}{(\mu + \theta)^{\sum{x_{1i} + n_1 \theta}}} \frac{(r \mu)^{\sum{x_{2j}}}}{(r \mu + \theta)^{\sum{x_{2j} + n_2 \theta}}} \times \\
#'   & \prod_{i = 1}^{n_1}{\frac{\Gamma(x_{1i} + \theta)}{x_{1i}!}} \prod_{j = 1}^{n_2}{\frac{\Gamma(x_{2j} + \theta)}{x_{2j}!}}
#' \end{aligned}
#' }
#'
#' and the parameter space is
#' \eqn{\Theta = \left\{ (r, \mu, \theta) : r, \mu, \theta > 0 \right\}}.
#' The log-likelihood is
#'
#' \deqn{
#' \begin{aligned}
#' l(r, \mu, \theta) = \ &(n_1 + n_2) \left[ \theta \ln \theta - \ln \Gamma(\theta) \right] + \\
#'   &(n_1 \bar{x}_1 + n_2 \bar{x}_2) \ln(\mu) - n_1 (\bar{x}_1 + \theta) \ln(\mu + \theta) + \\
#'   &n_2 \bar{x}_2 \ln(r) - n_2 (\bar{x}_2 + \theta) \ln(r \mu + \theta) + \\
#'   &\sum_{i = 1}^{n_1}{\left( \ln \Gamma(x_{1i} + \theta) - \ln(x_{1i}!) \right)} + \\
#'   &\sum_{j = 1}^{n_2}{\left( \ln \Gamma(x_{2j} + \theta) - \ln(x_{2j}!) \right)}
#' \end{aligned}
#' }
#'
#' @references
#' \insertRef{rettiganti_2012}{depower}
#'
#' \insertRef{aban_2009}{depower}
#'
#' @param param (numeric: `(0, Inf)`)\cr
#'        A vector of NB parameters. Must be in the following order for each scenario:
#' - Null and unequal dispersion: `c(mean, dispersion1, dispersion2)`
#' - Alternative and unequal dispersion: `c(mean1, mean2, dispersion1, dispersion2)`
#' - Null and equal dispersion: `c(mean, dispersion)`
#' - Alternative and equal dispersion: `c(mean1, mean2, dispersion)`
#'
#' for groups 1 and 2.
#' @param value1 (integer: `(0, Inf)`)\cr
#'        The vector of NB values from group 1. Must not contain \link[base]{NA}s.
#' @param value2 (integer: `(0, Inf)`)\cr
#'        The vector of NB values from group 2. Must not contain \link[base]{NA}s.
#' @param equal_dispersion (Scalar logical)\cr
#'        If `TRUE`, the log-likelihood is calculated assuming both groups have
#'        the same population dispersion parameter. If `FALSE` (default), the
#'        log-likelihood is calculated assuming different dispersions.
#' @param ratio_null (Scalar numeric: `(0, Inf)`)\cr
#'        The ratio of means assumed under the null hypothesis (group 2 / group 1).
#'        Typically `ratio_null = 1` (no difference).
#'
#' @return Scalar numeric negative log-likelihood.
#'
#' @seealso [depower::mle_nb]
#'
#' @examples
#' #----------------------------------------------------------------------------
#' # nll_nb_*() examples
#' #----------------------------------------------------------------------------
#' library(depower)
#'
#' set.seed(1234)
#' d <- sim_nb(
#'   n1 = 60,
#'   n2 = 40,
#'   mean1 = 10,
#'   ratio = 1.5,
#'   dispersion1 = 2,
#'   dispersion2 = 8
#' )
#'
#' nll_nb_alt(
#'   param = c(mean1 = 10, mean2 = 15, dispersion1 = 2, dispersion2 = 8),
#'   value1 = d[[1L]],
#'   value2 = d[[2L]],
#'   equal_dispersion = FALSE
#' )
#'
#' nll_nb_null(
#'   param = c(mean = 10, dispersion1 = 2, dispersion2 = 8),
#'   value1 = d[[1L]],
#'   value2 = d[[2L]],
#'   equal_dispersion = FALSE,
#'   ratio_null = 1
#' )
#'
#' @importFrom stats dnbinom
#'
#' @name nll_nb
NULL

#' @export
#' @rdname nll_nb
nll_nb_null <- function(param, value1, value2, equal_dispersion, ratio_null) {
  mean1 <- param[1L]
  mean2 <- mean1 * ratio_null

  if (equal_dispersion) {
    if (length(param) != 2L) {
      stop("Argument 'param' must be length 2 when 'equal_dispersion = TRUE'.")
    }
    dispersion1 <- param[2L]
    dispersion2 <- dispersion1
  } else {
    if (length(param) != 3L) {
      stop("Argument 'param' must be length 3 when 'equal_dispersion = FALSE'.")
    }
    dispersion1 <- param[2L]
    dispersion2 <- param[3L]
  }

  ll <- sum(dnbinom(x = value1, mu = mean1, size = dispersion1, log = TRUE)) +
    sum(dnbinom(x = value2, mu = mean2, size = dispersion2, log = TRUE))

  # Return negative log-likelihood
  -ll
}

#' @export
#' @rdname nll_nb
nll_nb_alt <- function(param, value1, value2, equal_dispersion) {
  mean1 <- param[1L]
  mean2 <- param[2L]

  if (equal_dispersion) {
    if (length(param) != 3L) {
      stop("Argument 'param' must be length 3 when 'equal_dispersion = TRUE'.")
    }
    dispersion1 <- param[3L]
    dispersion2 <- dispersion1
  } else {
    if (length(param) != 4L) {
      stop("Argument 'param' must be length 4 when 'equal_dispersion = FALSE'.")
    }
    dispersion1 <- param[3L]
    dispersion2 <- param[4L]
  }

  ll <- sum(dnbinom(x = value1, mu = mean1, size = dispersion1, log = TRUE)) +
    sum(dnbinom(x = value2, mu = mean2, size = dispersion2, log = TRUE))

  # Return negative log-likelihood
  -ll
}
