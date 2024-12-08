#' Wald test for BNB ratio of means
#'
#' Wald test for the ratio of means from bivariate negative binomial outcomes.
#'
#' This function is primarily designed for speed in simulation. Missing values
#' are silently excluded.
#'
#' Suppose \eqn{X_1 \mid G = g \sim \text{Poisson}(\mu g)} and
#' \eqn{X_2 \mid G = g \sim \text{Poisson}(r \mu g)} where
#' \eqn{G \sim \text{Gamma}(\theta, \theta^{-1})} is the random item (subject)
#' effect. Then \eqn{X_1, X_2 \sim \text{BNB}(\mu, r, \theta)} is the joint
#' distribution where \eqn{X_1} and \eqn{X_2} are dependent (though conditionally
#' independent), \eqn{X_1} is the count outcome for sample 1 of the items
#' (subjects), \eqn{X_2} is the count outcome for sample 2 of the items (subjects),
#' \eqn{\mu} is the conditional mean of sample 1, \eqn{r} is the ratio of the
#' conditional means of sample 2 with respect to sample 1, and \eqn{\theta} is
#' the gamma distribution shape parameter which controls the dispersion and the
#' correlation between sample 1 and 2.
#'
#' The hypotheses for the Wald test of \eqn{r} are
#'
#' \deqn{
#' \begin{aligned}
#' H_{null} &: f(r) = f(r_{null}) \\
#' H_{alt} &: f(r) \neq f(r_{null})
#' \end{aligned}
#' }
#'
#' where \eqn{f(\cdot)} is a one-to-one link function with nonzero derivative,
#' \eqn{r = \frac{\bar{X}_2}{\bar{X}_1}} is the population ratio of arithmetic
#' means for sample 2 with respect to sample 1, and \eqn{r_{null}} is a constant
#' for the assumed null population ratio of means (typically \eqn{r_{null} = 1}).
#'
#' \insertCite{rettiganti_2012;textual}{depower} found that \eqn{f(r) = r^2},
#' \eqn{f(r) = r}, and \eqn{f(r) = r^{0.5}} had greatest power when \eqn{r < 1}.
#' However, when \eqn{r > 1}, \eqn{f(r) = \ln r}, the likelihood ratio test, and
#' \eqn{f(r) = r^{0.5}} had greatest power. \eqn{f(r) = r^2} was biased when
#' \eqn{r > 1}. Both \eqn{f(r) = \ln r} and \eqn{f(r) = r^{0.5}} produced
#' acceptable results for any \eqn{r} value. These results depend on the use of
#' asymptotic vs. exact critical values.
#'
#' The Wald test statistic is
#'
#' \deqn{
#' W(f(\hat{r})) = \left( \frac{f \left( \frac{\bar{x}_2}{\bar{x}_1} \right) - f(r_{null})}{f^{\prime}(\hat{r}) \hat{\sigma}_{\hat{r}}} \right)^2
#' }
#'
#' where
#'
#' \deqn{
#' \hat{\sigma}^{2}_{\hat{r}} = \frac{\hat{r} (1 + \hat{r}) (\hat{\mu} + \hat{r}\hat{\mu} + \hat{\theta})}{n \left[ \hat{\mu} (1 + \hat{r}) (\hat{\mu} + \hat{\theta}) - \hat{\theta}\hat{r} \right]}
#' }
#'
#' Under \eqn{H_{null}}, the Wald test statistic is asymptotically distributed
#' as \eqn{\chi^2_1}. The approximate level \eqn{\alpha} test rejects
#' \eqn{H_{null}} if \eqn{W(f(\hat{r})) \geq \chi^2_1(1 - \alpha)}. Note that
#' the asymptotic critical value is known to underestimate the exact critical
#' value. Hence, the nominal significance level may not be achieved for small
#' sample sizes (possibly \eqn{n \leq 10} or \eqn{n \leq 50}). The level of
#' significance inflation also depends on \eqn{f(\cdot)} and is most severe for
#' \eqn{f(r) = r^2}, where only the exact critical value is recommended.
#'
#' @references
#' \insertRef{rettiganti_2012}{depower}
#'
#' \insertRef{aban_2009}{depower}
#'
#' @param data (list)\cr
#'        A list whose first element is the vector of negative binomial values
#'        from sample 1 and the second element is the vector of negative
#'        binomial values from sample 2.
#'        Each vector must be sorted by the subject/item index and must be the
#'        same sample size. \link[base]{NA}s are silently excluded. The default
#'        output from [depower::sim_bnb()].
#' @param ci_level (Scalar numeric: `NULL`; `(0, 1)`)\cr
#'        If `NULL`, confidence intervals are set as `NA`. If in `(0, 1)`,
#'        confidence intervals are calculated at the specified level.
#' @param link (Scalar string: `"log"`)\cr
#'        The one-to-one link function for transformation of the ratio in the
#'        test hypotheses. Must be one of `"log"` (default), `"sqrt"`,
#'        `"squared"`, or "`identity"`. See 'Details' for additional information.
#' @param ratio_null (Scalar numeric: `1`; `(0, Inf)`)\cr
#'        The (pre-transformation) ratio of means assumed under the null
#'        hypothesis (sample 2 / sample 1). Typically `ratio_null = 1`
#'        (no difference). See 'Details' for additional information.
#' @param ... Optional arguments passed to the MLE function [depower::mle_bnb()].
#'
#' @return
#' A list with the following elements:
#' \tabular{llll}{
#'   Slot \tab Subslot \tab Name \tab Description \cr
#'   1 \tab \tab `chisq` \tab Asymptotic \eqn{\chi^2} test statistic for the ratio of means. \cr
#'   2 \tab \tab `df` \tab Degrees of freedom. \cr
#'   3 \tab \tab `p` \tab p-value. \cr
#'
#'   4 \tab \tab `ratio`      \tab Estimated ratio of means (group 2 / group 1). \cr
#'   4 \tab 1 \tab `estimate` \tab Point estimate. \cr
#'   4 \tab 2 \tab `lower`    \tab Confidence interval lower bound. \cr
#'   4 \tab 3 \tab `upper`    \tab Confidence interval upper bound. \cr
#'
#'   5 \tab \tab `mean1` \tab Estimated mean of sample 1. \cr
#'   6 \tab \tab `mean2` \tab Estimated mean of sample 2. \cr
#'   7 \tab \tab `dispersion` \tab Estimated dispersion. \cr
#'   8 \tab \tab `n1` \tab The sample size of sample 1. \cr
#'   9 \tab \tab `n2` \tab The sample size of sample 2. \cr
#'   10 \tab \tab `method` \tab Method used for the results. \cr
#'   11 \tab \tab `ci_level` \tab The confidence level. \cr
#'   12 \tab \tab `link` \tab Link function used to transform the ratio of means in the test hypotheses. \cr
#'   13 \tab \tab `ratio_null` \tab Assumed ratio of means under the null hypothesis. \cr
#'   14 \tab \tab `mle_code` \tab Integer indicating why the optimization process terminated. \cr
#'   15 \tab \tab `mle_message` \tab Information from the optimizer.
#' }
#'
#' @examples
#' #----------------------------------------------------------------------------
#' # wald_test_bnb() examples
#' #----------------------------------------------------------------------------
#' library(depower)
#'
#' set.seed(1234)
#' sim_bnb(
#'   n = 40,
#'   mean1 = 10,
#'   ratio = 1.2,
#'   dispersion = 2
#' ) |>
#'   wald_test_bnb()
#'
#' @importFrom stats pchisq qnorm
#'
#' @export
wald_test_bnb <- function(
    data,
    ci_level = NULL,
    link = "log",
    ratio_null = 1,
    ...
) {
  #-----------------------------------------------------------------------------
  # Check args
  #-----------------------------------------------------------------------------
  if(!(is.list(data) && length(data) == 2L)) {
    stop("Argument 'data' must be a list with 2 elements.")
  }
  if(!is.null(ci_level)) {
    if(!is.numeric(ci_level) || length(ci_level) != 1L || ci_level <= 0 || ci_level >= 1) {
      stop("Argument 'ci_level' must be a scalar numeric from (0, 1).")
    }
  }
  if(!(length(ratio_null) == 1L && ratio_null > 0)) {
    stop("Argument 'ratio_null' must be a positive scalar numeric.")
  }

  #-----------------------------------------------------------------------------
  # Get MLE
  #-----------------------------------------------------------------------------
  mle_alt <- mle_bnb_alt(data, ...)

  #-----------------------------------------------------------------------------
  # Wald test
  #-----------------------------------------------------------------------------
  # Get components
  mle_mean1 <- mle_alt[["mean1"]]
  mle_mean2 <- mle_alt[["mean2"]]
  mle_ratio <- mle_alt[["ratio"]]
  mle_dispersion <- mle_alt[["dispersion"]]
  n1 <- mle_alt[["n1"]]
  n2 <- mle_alt[["n2"]]

  sigma <- sqrt(
    mle_ratio * (1 + mle_ratio) * (mle_mean1 + mle_mean2 + mle_dispersion) /
      (n1 * ((mle_mean1 + mle_mean2) * (mle_mean1 + mle_dispersion) - mle_dispersion * mle_ratio))
  )

  # Link functions
  linkf <- switch(link,
    "log" = log,
    "sqrt" = sqrt,
    "squared" = function(x) x^2,
    "identity" = function(x) x,
    stop("Argument 'link' must be one of 'log', 'sqrt', 'squared', or 'identity'.")
  )
  # See ?D and ?deriv
  derivative <- switch(link,
    "log" = function(x) 1/x,
    "sqrt" = function(x) 1 / (2 * sqrt(x)),
    "squared" = function(x) 2 * x,
    "identity" = function(x) 1,
    stop("Argument 'link' must be one of 'log', 'sqrt', 'squared', or 'identity'.")
  )
  inverse <- switch(link,
    "log" = exp,
    "sqrt" = function(x) x^2,
    "squared" = sqrt,
    "identity" = function(x) x,
    stop("Argument 'link' must be one of 'log', 'sqrt', 'squared', or 'identity'.")
  )

  # Test statistic
  wald_stat <- ((linkf(mle_ratio) - linkf(ratio_null)) / (derivative(mle_ratio) * sigma))^2
  wald_df <- 1L
  wald_p <- pchisq(wald_stat, df = wald_df, lower.tail = FALSE)

  # CI
  if(is.null(ci_level)) {
    r_lower <- NA_real_
    r_upper <- NA_real_
  } else {
    alpha <- (1 - ci_level) / 2
    crit <- qnorm(p = alpha, lower.tail = FALSE)
    r_lower <- inverse(linkf(mle_ratio) - crit * derivative(mle_ratio) * sigma)
    r_upper <- inverse(linkf(mle_ratio) + crit * derivative(mle_ratio) * sigma)
  }

  #-----------------------------------------------------------------------------
  # Prepare result
  #-----------------------------------------------------------------------------
  method <- "Wald test for bivariate negative binomial ratio of means"
  mle_code <- mle_alt[["mle_code"]]
  mle_message <- mle_alt[["mle_message"]]

  #-----------------------------------------------------------------------------
  # Return
  #-----------------------------------------------------------------------------
  list(
    chisq = wald_stat,
    df = wald_df,
    p = wald_p,
    ratio = list(
      estimate = mle_ratio,
      lower = r_lower,
      upper = r_upper
    ),
    mean1 = mle_mean1,
    mean2 = mle_mean2,
    dispersion = mle_dispersion,
    n1 = n1,
    n2 = n2,
    method = method,
    ci_level = ci_level,
    link = link,
    ratio_null = ratio_null,
    mle_code = mle_code,
    mle_message = mle_message
  )
}
