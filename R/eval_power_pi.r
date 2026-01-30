#' @title
#' Evaluate Bayesian posterior predictive intervals for power estimates
#'
#' @description
#' Calculates the Bayesian posterior predictive interval for a power estimate from a simulation study.
#' The posterior predictive interval quantifies the expected range of power estimates from a future simulation study.
#'
#' When the number of simulations used to calculate a test's power is too small, the power estimate will have high uncertainty (wide confidence/prediction intervals).
#' When the number of simulations used to calculate a test's power is too large, computational time may be prohibitive.
#' This function allows you to determine the appropriate number of simulated datasets to reach your desired precision for power before spending computational time on simulations.
#'
#' @details
#' Power estimation via simulation is a binomial proportion problem.
#' The posterior predictive interval answers: "If I run a new simulation study with \eqn{m} simulations, what range of power estimates might I observe?"
#'
#' Let \eqn{\pi} denote the hypothetical true power value, \eqn{\hat{\pi} = x/n} denote the hypothetical observed power value, \eqn{n} denote the number of simulations, and \eqn{x = \text{round}(\hat{\pi} \cdot n)} denote the number of rejections.
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
#' @param power
#' (numeric: `(0, 1)`)\cr
#' Hypothetical power value(s).
#'
#' @param nsims
#' (integer: `[2, Inf)`)\cr
#' Number of simulations.
#'
#' @param future_nsims
#' (integer or `NULL`: `NULL`; `[2, Inf)`)\cr
#' Number of simulations in the future study.
#' If `NULL` (default), uses the same number as `nsims`.
#'
#' @param pi_level
#' (Scalar numeric: `0.95`; `(0,1)`)\cr
#' The posterior predictive interval level.
#'
#' @param prior
#' (numeric vector of length 2: `c(1, 1)`; each `(0, Inf)`)\cr
#' Parameters \eqn{(\alpha, \beta)} for the Beta prior on true power.
#' Default `c(1, 1)` is the uniform prior.
#' Use `c(0.5, 0.5)` for the Jeffreys prior.
#'
#' @returns
#' A list with elements:
#' \tabular{ll}{
#'   Name \tab Description \cr
#'   `mean` \tab Predictive mean of future power estimate. \cr
#'   `lower` \tab Lower bound of posterior predictive interval. \cr
#'   `upper` \tab Upper bound of posterior predictive interval.
#' }
#'
#' @seealso
#' [depower::add_power_pi()],
#' [depower::eval_power_ci()]
#'
#' @examples
#' #----------------------------------------------------------------------------
#' # eval_power_pi() examples
#' #----------------------------------------------------------------------------
#' library(depower)
#'
#' # Expected PI for 80% power with 1000 simulations
#' eval_power_pi(power = 0.80, nsims = 1000)
#'
#' # Compare precision across different simulation counts
#' eval_power_pi(power = 0.80, nsims = c(100, 500, 1000, 5000))
#'
#' # Predict for a larger future study (narrower interval)
#' eval_power_pi(power = 0.80, nsims = 1000, future_nsims = 5000)
#'
#' # Predict for a smaller future study (wider interval)
#' eval_power_pi(power = 0.80, nsims = 1000, future_nsims = 100)
#'
#' # Vectorized over power values
#' eval_power_pi(power = c(0.70, 0.80, 0.90), nsims = 1000)
#'
#' # Use Jeffreys prior instead of uniform
#' eval_power_pi(power = 0.80, nsims = 1000, prior = c(0.5, 0.5))
#'
#' # 99% predictive interval
#' eval_power_pi(power = 0.80, nsims = 1000, pi_level = 0.99)
#'
#' @export
eval_power_pi <- function(
  power,
  nsims,
  future_nsims = NULL,
  pi_level = 0.95,
  prior = c(1, 1)
) {
  #-----------------------------------------------------------------------------
  # Check arguments
  #-----------------------------------------------------------------------------
  if (!is.numeric(power) || any(power <= 0) || any(power >= 1)) {
    stop("Argument 'power' must be numeric in (0, 1).")
  }

  if (!is.numeric(nsims) || any(nsims <= 1) || any(nsims != round(nsims))) {
    stop("Argument 'nsims' must be integers greater than 1.")
  }

  if (
    length(power) != length(nsims) && length(power) != 1L && length(nsims) != 1L
  ) {
    msg <- "Arguments 'power' and 'nsims' must be the same length, or one must be length 1."
    stop(msg)
  }

  if (!is.null(future_nsims)) {
    if (
      !is.numeric(future_nsims) ||
        any(future_nsims <= 1) ||
        any(future_nsims != round(future_nsims))
    ) {
      stop("Argument 'future_nsims' must be integers greater than 1 or NULL.")
    }
    if (length(nsims) != length(future_nsims)) {
      stop("Arguments 'nsims' and 'future_nsims' must be the same length.")
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
  # Recycle to common length
  max_len <- max(length(power), length(nsims))
  power <- rep_len(power, max_len)
  nsims <- rep_len(nsims, max_len)

  n_rejects <- round(power * nsims)

  res <- binom_pi_bayes(
    x = n_rejects,
    n = nsims,
    future_n = future_nsims,
    pred_level = pi_level,
    prior = prior
  )
  res$future_n <- NULL

  #-----------------------------------------------------------------------------
  # Return
  #-----------------------------------------------------------------------------
  res
}
