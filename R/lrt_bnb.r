#' Likelihood ratio test for BNB ratio of means
#'
#' Likelihood ratio test for the ratio of means from bivariate negative binomial
#' outcomes.
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
#' The hypotheses for the LRT of \eqn{r} are
#'
#' \deqn{
#' \begin{aligned}
#' H_{null} &: r = r_{null} \\
#' H_{alt} &: r \neq r_{null}
#' \end{aligned}
#' }
#'
#' where \eqn{r = \frac{\bar{X}_2}{\bar{X}_1}} is the population ratio of
#' arithmetic means for sample 2 with respect to sample 1 and \eqn{r_{null}} is
#' a constant for the assumed null population ratio of means (typically
#' \eqn{r_{null} = 1}).
#'
#' The LRT statistic is
#'
#' \deqn{
#' \begin{aligned}
#' \lambda &= -2 \ln \frac{\text{sup}_{\Theta_{null}} L(r, \mu, \theta)}{\text{sup}_{\Theta} L(r, \mu, \theta)} \\
#'   &= -2 \left[ \ln \text{sup}_{\Theta_{null}} L(r, \mu, \theta) - \ln \text{sup}_{\Theta} L(r, \mu, \theta) \right] \\
#'   &= -2(l(r_{null}, \tilde{\mu}, \tilde{\theta}) - l(\hat{r}, \hat{\mu}, \hat{\theta}))
#' \end{aligned}
#' }
#'
#' Under \eqn{H_{null}}, the LRT test statistic is asymptotically distributed
#' as \eqn{\chi^2_1}. The approximate level \eqn{\alpha} test rejects
#' \eqn{H_{null}} if \eqn{\lambda \geq \chi^2_1(1 - \alpha)}. However,
#' the asymptotic critical value is known to underestimate the exact critical
#' value and the nominal significance level may not be achieved for small sample
#' sizes. Argument `distribution` allows control of the distribution of
#' the \eqn{\chi^2_1} test statistic under the null hypothesis by use of
#' functions [depower::asymptotic()] and [depower::simulated()].
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
#' @param ratio_null (Scalar numeric: `1`; `(0, Inf)`)\cr
#'        The ratio of means assumed under the null hypothesis (sample 2 / sample 1).
#'        Typically, `ratio_null = 1` (no difference). See 'Details' for
#'        additional information.
#' @param distribution (function: [depower::asymptotic()] or
#'        [depower::simulated()])\cr
#'        The method used to define the distribution of the \eqn{\chi^2}
#'        likelihood ratio test statistic under the null hypothesis. See
#'        'Details' and [depower::asymptotic()] or [depower::simulated()] for
#'        additional information.
#' @param ... Optional arguments passed to the MLE function [depower::mle_bnb()].
#'
#' @return
#' A list with the following elements:
#' \tabular{llll}{
#'   Slot \tab Subslot \tab Name \tab Description \cr
#'   1 \tab \tab `chisq` \tab \eqn{\chi^2} test statistic for the ratio of means. \cr
#'   2 \tab \tab `df`    \tab Degrees of freedom. \cr
#'   3 \tab \tab `p`     \tab p-value. \cr
#'   4 \tab \tab `ratio` \tab Estimated ratio of means (sample 2 / sample 1). \cr
#'
#'   5 \tab   \tab `alternative` \tab Point estimates under the alternative hypothesis. \cr
#'   5 \tab 1 \tab `mean1` \tab Estimated mean of sample 1. \cr
#'   5 \tab 2 \tab `mean2` \tab Estimated mean of sample 2. \cr
#'   5 \tab 3 \tab `dispersion` \tab Estimated dispersion. \cr
#'
#'   6 \tab   \tab `null` \tab Point estimates under the null hypothesis. \cr
#'   6 \tab 1 \tab `mean1` \tab Estimated mean of sample 1. \cr
#'   6 \tab 2 \tab `mean2` \tab Estimated mean of sample 2. \cr
#'   6 \tab 3 \tab `dispersion` \tab Estimated dispersion. \cr
#'
#'   7 \tab \tab `n1` \tab The sample size of sample 1. \cr
#'   8 \tab \tab `n2` \tab The sample size of sample 2. \cr
#'   9 \tab \tab `method` \tab Method used for the results. \cr
#'   10 \tab \tab `ratio_null` \tab Assumed population ratio of means. \cr
#'   11 \tab \tab `mle_code` \tab Integer indicating why the optimization process terminated. \cr
#'   12 \tab \tab `mle_message` \tab Information from the optimizer.
#' }
#'
#' @seealso [depower::wald_test_bnb()]
#'
#' @examples
#' #----------------------------------------------------------------------------
#' # lrt_bnb() examples
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
#'   lrt_bnb()
#'
#' @importFrom stats pchisq
#'
#' @export
lrt_bnb <- function(data, ratio_null = 1, distribution = asymptotic(), ...) {
  #-----------------------------------------------------------------------------
  # Check args
  #-----------------------------------------------------------------------------
  if (!(is.list(data) && length(data) == 2L)) {
    stop("Argument 'data' must be a list with 2 elements.")
  }
  if (any(data[[1L]] < 0, na.rm = TRUE) || any(data[[2L]] < 0, na.rm = TRUE)) {
    stop("Argument 'data' must not contain negative numbers.")
  }
  if (!(length(ratio_null) == 1L && ratio_null > 0)) {
    stop("Argument 'ratio_null' must be a positive scalar numeric.")
  }
  if (
    !is.list(distribution) ||
      length(distribution$distribution) != 1L ||
      !distribution$distribution %in% c("asymptotic", "simulated")
  ) {
    stop(
      "Argument 'distribution' must be a function from 'depower::asymptotic()' or 'depower::simulated()'."
    )
  }

  #-----------------------------------------------------------------------------
  # Get MLEs
  #-----------------------------------------------------------------------------
  mle_null <- mle_bnb_null(data, ratio_null = ratio_null, ...)
  mle_alt <- mle_bnb_alt(data, ...)

  #-----------------------------------------------------------------------------
  # LRT
  #-----------------------------------------------------------------------------
  nll_null <- mle_null[["nll"]]
  nll_alt <- mle_alt[["nll"]]

  lrt_stat <- 2 * (nll_null - nll_alt)
  lrt_stat <- max(0, lrt_stat) # This might result in a negative number.
  lrt_df <- 1L # mle_alt[["nparams"]] - mle_null[["nparams"]]

  # Randomization test if needed
  lrt_stat_null <- if (distribution$distribution == "simulated") {
    if (is.null(distribution$test_stat_null)) {
      test_stat_method <- "randomization"
      call <- match.call()
      randomize_tests(
        data = data,
        distribution = distribution,
        paired = TRUE,
        call = call
      ) |>
        extract(column = "result", element = "chisq", value = numeric(1L))
    } else {
      test_stat_method <- "parametric"
      distribution$test_stat_null
    }
  } else {
    test_stat_method <- NULL
    NULL
  }

  # p-value
  lrt_p <- pchisq2(
    q = lrt_stat,
    df = lrt_df,
    lower.tail = FALSE,
    q_null = lrt_stat_null
  )

  #-----------------------------------------------------------------------------
  # Prepare result
  #-----------------------------------------------------------------------------
  mean1_alt = mle_alt[["mean1"]]
  mean2_alt = mle_alt[["mean2"]]
  ratio_alt <- mle_alt[["ratio"]]
  dispersion_alt <- mle_alt[["dispersion"]]

  mean1_null = mle_null[["mean1"]]
  mean2_null = mle_null[["mean2"]]
  dispersion_null <- mle_null[["dispersion"]]

  n1 <- mle_null[["n1"]]
  n2 <- mle_null[["n2"]]
  mle_code <- mle_alt[["mle_code"]]
  mle_message <- mle_alt[["mle_message"]]

  method <- paste0(
    str_to_title(distribution$method),
    if (is.null(test_stat_method)) NULL else paste0(" ", test_stat_method),
    " LRT for bivariate negative binomial ratio of means"
  )

  #-----------------------------------------------------------------------------
  # Return
  #-----------------------------------------------------------------------------
  list(
    chisq = lrt_stat,
    df = lrt_df,
    p = lrt_p,
    ratio = ratio_alt,
    alternative = list(
      mean1 = mean1_alt,
      mean2 = mean2_alt,
      dispersion = dispersion_alt
    ),
    null = list(
      mean1 = mean1_null,
      mean2 = mean2_null,
      dispersion = dispersion_null
    ),
    n1 = n1,
    n2 = n2,
    method = method,
    ratio_null = ratio_null,
    mle_code = mle_code,
    mle_message = mle_message
  )
}
