#' @title
#' Add Bayesian posterior predictive intervals for power estimates
#'
#' @description
#' Calculates and adds Bayesian posterior predictive intervals for power estimates in objects returned by [depower::power()].
#' The posterior predictive interval quantifies the expected range of power estimates from a future simulation study.
#'
#' @details
#' Power estimation via simulation is a binomial proportion problem.
#' The posterior predictive interval answers: "If I run a new simulation study with \eqn{m} simulations, what range of power estimates might I observe?"
#'
#' Let \eqn{\pi} denote the true power value, \eqn{\hat{\pi} = x/n} denote the observed power value, \eqn{n} denote the number of simulations, and \eqn{x = \text{round}(\hat{\pi} \cdot n)} denote the number of rejections.
#' With a \eqn{\text{Beta}(\alpha, \beta)} prior on the true power \eqn{\pi}, the posterior after observing \eqn{x} successes in \eqn{n} trials is:
#'
#' \deqn{
#'   \pi \mid X = x \sim \text{Beta}(\alpha + x, \beta + n - x).
#' }
#'
#' The posterior predictive distribution for \eqn{Y}, the number of successes in a future study with \eqn{m} trials, is Beta-Binomial:
#'
#' \deqn{
#'   Y \mid X = x \sim \text{BetaBinomial}(m, \alpha + x, \beta + n - x).
#' }
#'
#' The posterior predictive interval is constructed from quantiles of this distribution, expressed as proportions \eqn{Y/m}.
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
#' ## Argument `future_nsims`
#'
#' The argument `future_nsims` allows you to estimate prediction interval bounds for a hypothetical future study with different number of simulations.
#' Note that a small initial number for `nsims` results in substantial uncertainty about the true power.
#' A correspondingly large number of future simulations `future_nsims` will more precisely estimate the true power, but the past large uncertainty is still carried forward.
#' Therefore you still need an adequate number of simulations `nsims` in the original study, not just more in the replication `future_nsims`, to ensure narrow prediction intervals.
#'
#' ## Approximate parametric tests
#'
#' When power is computed using approximate parametric tests (see [depower::simulated()]), the power estimate and confidence/prediction intervals apply to the Monte Carlo test power \eqn{\mu_K = P(\hat{p} \leq \alpha)} rather than the exact test power \eqn{\pi = P(p \leq \alpha)}.
#' These quantities converge as the number of datasets simulated under the null hypothesis \eqn{K} increases.
#' The minimum observable p-value is \eqn{1/(K+1)}, so \eqn{K > 1/\alpha - 1} is required to observe any rejections.
#' For practical accuracy, we recommend choosing \eqn{\text{max}(5000, K \gg 1/\alpha - 1)} for most scenarios.
#' For example, if \eqn{\alpha = 0.05}, use `simulated(nsims = 5000)`.
#'
#' @references
#' \insertRef{gelman_2013}{depower}
#'
#' @param x
#' (data.frame)\cr
#' A data frame returned by [depower::power()], containing columns `power` and `nsims`.
#'
#' @param future_nsims
#' (Scalar integer or `NULL`: `NULL`; `[2, Inf)`)\cr
#' Number of simulated data sets in the future power estimate study.
#' If `NULL` (default), uses the same number as the original study (`nsims`).
#'
#' @param pi_level
#' (Scalar numeric: `0.95`; `(0,1)`)\cr
#' The posterior predictive interval level.
#'
#' @param prior
#' (Numeric vector of length 2: `c(1, 1)`; each `(0, Inf)`)\cr
#' Parameters \eqn{(\alpha, \beta)} for the Beta prior on true power.
#' Default `c(1, 1)` is the uniform prior.
#' Use `c(0.5, 0.5)` for the Jeffreys prior.
#'
#' @returns
#' The input data frame with additional columns:
#' \tabular{ll}{
#'   Name \tab Description \cr
#'   `power_pi_mean`  \tab Predictive mean of future power estimate. \cr
#'   `power_pi_lower` \tab Lower bound of posterior predictive interval. \cr
#'   `power_pi_upper` \tab Upper bound of posterior predictive interval.
#' }
#' and added attribute `"pi_info"` containing the method description, method name, level, prior values, and future simulation count.
#'
#' @seealso
#' [depower::power()],
#' [depower::eval_power_pi()],
#' [depower::add_power_ci()]
#'
#' @examples
#' #----------------------------------------------------------------------------
#' # add_power_pi() examples
#' #----------------------------------------------------------------------------
#' library(depower)
#'
#' set.seed(1234)
#' x <- sim_nb(
#'   n1 = 10,
#'   mean1 = 10,
#'   ratio = c(1.4, 1.6),
#'   dispersion1 = 2,
#'   nsims = 200
#' ) |>
#'   power(wald_test_nb())
#'
#' # Add posterior predictive intervals
#' # default: predict for same number of simulations
#' add_power_pi(x)
#'
#' # Compare posterior predictive interval width across different future
#' # study sizes
#' add_power_pi(x, future_nsims = 100)  # wider
#' add_power_pi(x, future_nsims = 1000) # narrower
#'
#' # Use Jeffreys prior instead of uniform
#' add_power_pi(x, prior = c(0.5, 0.5))
#'
#' # Plot with shaded region for prediction interval of the power estimate.
#' add_power_pi(x) |>
#'   plot()
#'
#' @export
add_power_pi <- function(
  x,
  future_nsims = NULL,
  pi_level = 0.95,
  prior = c(1, 1)
) {
  #-----------------------------------------------------------------------------
  # Check arguments
  #-----------------------------------------------------------------------------
  if (!is.data.frame(x)) {
    stop("Argument 'x' must be a data frame.")
  }
  if (!all(c("power", "nsims") %in% names(x))) {
    stop("Argument 'x' must contain columns 'power' and 'nsims'.")
  }

  if (!is.null(future_nsims)) {
    if (
      !is.numeric(future_nsims) ||
        length(future_nsims) != 1L ||
        future_nsims < 2L ||
        future_nsims != round(future_nsims)
    ) {
      msg <- "Argument 'future_nsims' must be a scalar integer greater than 1 or NULL."
      stop(msg)
    }
  }

  if (
    !is.numeric(pi_level) ||
      length(pi_level) != 1L ||
      pi_level <= 0 ||
      pi_level >= 1
  ) {
    stop("Argument 'pi_level' must be a scalar numeric in (0, 1).")
  }

  if (!is.numeric(prior) || length(prior) != 2L || any(prior <= 0)) {
    stop("Argument 'prior' must be a positive numeric vector of length 2.")
  }

  #-----------------------------------------------------------------------------
  # Compute posterior predictive intervals
  #-----------------------------------------------------------------------------
  n_rejects <- round(x$power * x$nsims)

  pi_result <- binom_pi_bayes(
    x = n_rejects,
    n = x$nsims,
    future_n = future_nsims,
    pred_level = pi_level,
    prior = prior
  )

  #-----------------------------------------------------------------------------
  # Add columns to result
  #-----------------------------------------------------------------------------
  x$power_pi_mean <- pi_result$mean
  x$power_pi_lower <- pi_result$lower
  x$power_pi_upper <- pi_result$upper

  #-----------------------------------------------------------------------------
  # Return
  #-----------------------------------------------------------------------------
  attr(x, "pi_info") <- list(
    description = "Bayesian posterior predictive interval",
    method = "Beta-Binomial predictive interval",
    level = pi_level,
    prior = prior,
    future_nsims = pi_result$future_n
  )
  x
}
