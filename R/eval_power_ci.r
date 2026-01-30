#' @title
#' Evaluate confidence intervals for power estimates
#'
#' @description
#' Calculates the confidence interval for a power estimate from a simulation study.
#' The confidence interval quantifies uncertainty about the true power parameter.
#'
#' When the number of simulations used to calculate a test's power is too small, the power estimate will have high uncertainty (wide confidence/prediction intervals).
#' When the number of simulations used to calculate a test's power is too large, computational time may be prohibitive.
#' This function allows you to determine the appropriate number of simulated datasets to reach your desired precision for power before spending computational time on simulations.
#'
#' @details
#' Power estimation via simulation is a binomial proportion problem.
#' The confidence interval answers: "What is the plausible range of true power values given my simulation results?"
#'
#' Let \eqn{\pi} denote the hypothetical true power value, \eqn{\hat{\pi} = x/n} denote the hypothetical observed power value, \eqn{n} denote the number of simulations, and \eqn{x = \text{round}(\hat{\pi} \cdot n)} denote the number of rejections.
#' Two methods are available.
#'
#' ## Wilson Score Interval
#'
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
#' ## Clopper-Pearson Interval
#'
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
#' ## Approximate parametric tests
#'
#' When power is computed using approximate parametric tests (see [depower::simulated()]), the power estimate and confidence/prediction intervals apply to the Monte Carlo test power \eqn{\mu_K = P(\hat{p} \leq \alpha)} rather than the exact test power \eqn{\pi = P(p \leq \alpha)}.
#' These quantities converge as the number of datasets simulated under the null hypothesis \eqn{K} increases.
#' The minimum observable p-value is \eqn{1/(K+1)}, so \eqn{K > 1/\alpha - 1} is required to observe any rejections.
#' For practical accuracy, we recommend choosing \eqn{\text{max}(5000, K \gg 1/\alpha - 1)} for most scenarios.
#' For example, if \eqn{\alpha = 0.05}, use `simulated(nsims = 5000)`.
#'
#' @references
#' \insertRef{newcombe_1998}{depower},
#'
#' \insertRef{wilson_1927}{depower},
#'
#' \insertRef{clopper_1934}{depower}
#'
#' @param power
#' (numeric: `(0, 1)`)\cr
#' Hypothetical observed power value(s).
#'
#' @param nsims
#' (integer: `[2, Inf)`)\cr
#' Number of simulations.
#'
#' @param ci_level
#' (Scalar numeric: `0.95`; `(0,1)`)\cr
#' The confidence level.
#'
#' @param method
#' (Scalar character: `"wilson"`; `c("wilson", "exact")`)\cr
#' Method for computing confidence intervals.
#' One of `"wilson"` (default) or `"exact"`.
#' See 'Details' for more information.
#'
#' @returns
#' A list with elements:
#' \tabular{ll}{
#'   Name \tab Description \cr
#'   `lower` \tab Lower bound of confidence interval. \cr
#'   `upper` \tab Upper bound of confidence interval.
#' }
#'
#' @seealso
#' [depower::add_power_ci()],
#' [depower::eval_power_pi()]
#'
#' @examples
#' #----------------------------------------------------------------------------
#' # eval_power_ci() examples
#' #----------------------------------------------------------------------------
#' library(depower)
#'
#' # Expected CI for 80% power with 1000 simulations
#' eval_power_ci(power = 0.80, nsims = 1000)
#'
#' # Compare precision across different simulation counts
#' eval_power_ci(power = 0.80, nsims = c(100, 500, 1000, 5000))
#'
#' # Compare Wilson vs exact method
#' eval_power_ci(power = 0.80, nsims = 1000, method = "wilson")
#' eval_power_ci(power = 0.80, nsims = 1000, method = "exact")
#'
#' # Vectorized over power values
#' eval_power_ci(power = c(0.70, 0.80, 0.90), nsims = 1000)
#'
#' # 99% confidence interval
#' eval_power_ci(power = 0.80, nsims = 1000, ci_level = 0.99)
#'
#' @export
eval_power_ci <- function(
  power,
  nsims,
  ci_level = 0.95,
  method = c("wilson", "exact")
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

  if (
    !is.numeric(ci_level) ||
      length(ci_level) != 1L ||
      ci_level <= 0 ||
      ci_level >= 1
  ) {
    stop("Argument 'ci_level' must be a scalar numeric in (0, 1).")
  }

  method <- match.arg(method)

  #-----------------------------------------------------------------------------
  # Compute confidence intervals
  #-----------------------------------------------------------------------------
  # Recycle to common length
  max_len <- max(length(power), length(nsims))
  power <- rep_len(power, max_len)
  nsims <- rep_len(nsims, max_len)

  n_rejects <- round(power * nsims)

  res <- if (method == "wilson") {
    binom_ci_wilson(
      x = n_rejects,
      n = nsims,
      conf_level = ci_level
    )
  } else {
    binom_ci_clopper_pearson(
      x = n_rejects,
      n = nsims,
      conf_level = ci_level
    )
  }

  #-----------------------------------------------------------------------------
  # Return
  #-----------------------------------------------------------------------------
  res
}
