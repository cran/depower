#' Likelihood ratio test for NB ratio of means
#'
#' Likelihood ratio test for the ratio of means from two independent negative
#' binomial outcomes.
#'
#' This function is primarily designed for speed in simulation. Missing values
#' are silently excluded.
#'
#' Suppose \eqn{X_1 \sim NB(\mu, \theta_1)} and
#' \eqn{X_2 \sim NB(r\mu, \theta_2)} where \eqn{X_1} and \eqn{X_2} are
#' independent, \eqn{X_1} is the count outcome for items in group 1, \eqn{X_2}
#' is the count outcome for items in group 2, \eqn{\mu} is the arithmetic mean
#' count in group 1, \eqn{r} is the ratio of arithmetic means for group 2 with
#' respect to group 1, \eqn{\theta_1} is the dispersion parameter of group 1,
#' and \eqn{\theta_2} is the dispersion parameter of group 2.
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
#' arithmetic means for group 2 with respect to group 1 and \eqn{r_{null}} is a
#' constant for the assumed null population ratio of means (typically
#' \eqn{r_{null} = 1}).
#'
#' The LRT statistic is
#'
#' \deqn{
#' \begin{aligned}
#' \lambda &= -2 \ln \frac{\text{sup}_{\Theta_{null}} L(r, \mu, \theta_1, \theta_2)}{\text{sup}_{\Theta} L(r, \mu, \theta_1, \theta_2)} \\
#'   &= -2 \left[ \ln \text{sup}_{\Theta_{null}} L(r, \mu, \theta_1, \theta_2) - \ln \text{sup}_{\Theta} L(r, \mu, \theta_1, \theta_2) \right] \\
#'   &= -2(l(r_{null}, \tilde{\mu}, \tilde{\theta}_1, \tilde{\theta}_2) - l(\hat{r}, \hat{\mu}, \hat{\theta}_1, \hat{\theta}_2))
#' \end{aligned}
#' }
#'
#' Under \eqn{H_{null}}, the LRT test statistic is asymptotically distributed
#' as \eqn{\chi^2_1}. The approximate level \eqn{\alpha} test rejects
#' \eqn{H_{null}} if \eqn{\lambda \geq \chi^2_1(1 - \alpha)}. Note that
#' the asymptotic critical value is known to underestimate the exact critical
#' value. Hence, the nominal significance level may not be achieved for small
#' sample sizes (possibly \eqn{n \leq 10} or \eqn{n \leq 50}).
#'
#' @references
#' \insertRef{rettiganti_2012}{depower}
#'
#' \insertRef{aban_2009}{depower}
#'
#' @param data (list)\cr
#'        A list whose first element is the vector of negative binomial values
#'        from group 1 and the second element is the vector of negative binomial
#'        values from group 2.
#'        \link[base]{NA}s are silently excluded. The default output from
#'        [depower::sim_nb()].
#' @param equal_dispersion (Scalar logical: `FALSE`)\cr
#'        If `TRUE`, the LRT is calculated assuming both groups have the same
#'        population dispersion parameter. If `FALSE` (default), the LRT is
#'        calculated assuming different dispersions.
#' @param ratio_null (Scalar numeric: `1`; `(0, Inf)`)\cr
#'        The ratio of means assumed under the null hypothesis (group 2 / group 1).
#'        Typically `ratio_null = 1` (no difference). See 'Details' for
#'        additional information.
#' @param ... Optional arguments passed to the MLE function [depower::mle_nb()].
#'
#' @return
#' A list with the following elements:
#' \tabular{llll}{
#'   Slot \tab Subslot \tab Name \tab Description \cr
#'   1 \tab  \tab `chisq` \tab \eqn{\chi^2} test statistic for the ratio of means. \cr
#'   2 \tab  \tab `df`    \tab Degrees of freedom. \cr
#'   3 \tab  \tab `p`     \tab p-value. \cr
#'   4 \tab  \tab `ratio` \tab Estimated ratio of means (group 2 / group 1). \cr
#'
#'   5 \tab   \tab `alternative` \tab Point estimates under the alternative hypothesis. \cr
#'   5 \tab 1 \tab `mean1`       \tab Estimated mean of group 1. \cr
#'   5 \tab 2 \tab `mean2`       \tab Estimated mean of group 2. \cr
#'   5 \tab 3 \tab `dispersion1` \tab Estimated dispersion of group 1. \cr
#'   5 \tab 4 \tab `dispersion2` \tab Estimated dispersion of group 2. \cr
#'
#'   6 \tab   \tab `null` \tab Point estimates under the null hypothesis. \cr
#'   6 \tab 1 \tab `mean1`       \tab Estimated mean of group 1. \cr
#'   6 \tab 2 \tab `mean2`       \tab Estimated mean of group 2. \cr
#'   6 \tab 3 \tab `dispersion1` \tab Estimated dispersion of group 1. \cr
#'   6 \tab 4 \tab `dispersion2` \tab Estimated dispersion of group 2. \cr
#'
#'   7 \tab \tab `n1` \tab Sample size of group 1. \cr
#'   8 \tab \tab `n2` \tab Sample size of group 2. \cr
#'   9 \tab \tab `method` \tab Method used for the results. \cr
#'   10 \tab \tab `equal_dispersion` \tab Whether or not equal dispersions were assumed. \cr
#'   11 \tab \tab `ratio_null` \tab Assumed population ratio of means. \cr
#'   12 \tab \tab `mle_code` \tab Integer indicating why the optimization process terminated. \cr
#'   13 \tab \tab `mle_message` \tab Information from the optimizer.
#' }
#'
#' @examples
#' #----------------------------------------------------------------------------
#' # lrt_nb() examples
#' #----------------------------------------------------------------------------
#' library(depower)
#'
#' set.seed(1234)
#' sim_nb(
#'   n1 = 60,
#'   n2 = 40,
#'   mean1 = 10,
#'   ratio = 1.5,
#'   dispersion1 = 2,
#'   dispersion2 = 8
#' ) |>
#'   lrt_nb()
#'
#' @importFrom stats pchisq
#'
#' @export
lrt_nb <- function(data, equal_dispersion = FALSE, ratio_null = 1, ...) {
  #-----------------------------------------------------------------------------
  # Check args
  #-----------------------------------------------------------------------------
  if(!(is.list(data) && length(data) == 2L)) {
    stop("Argument 'data' must be a list with 2 elements.")
  }
  if(!(is.logical(equal_dispersion) && length(equal_dispersion) == 1L)) {
    stop("Argument 'equal_dispersion' must be a scalar logical.")
  }
  if(!(length(ratio_null) == 1L && ratio_null > 0)) {
    stop("Argument 'ratio_null' must be a positive scalar numeric.")
  }

  #-----------------------------------------------------------------------------
  # Get MLEs
  #-----------------------------------------------------------------------------
  mle_null <- mle_nb_null(data, equal_dispersion = FALSE, ratio_null = ratio_null, ...)
  mle_alt <- mle_nb_alt(data, equal_dispersion = FALSE, ...)

  #-----------------------------------------------------------------------------
  # LRT
  #-----------------------------------------------------------------------------
  nll_null <- mle_null[["nll"]]
  nll_alt <- mle_alt[["nll"]]

  lrt_stat <- 2 * (nll_null - nll_alt)
  lrt_df <- 1L # mle_alt[["nparams"]] - mle_null[["nparams"]]
  lrt_p <- pchisq(lrt_stat, df = lrt_df, lower.tail = FALSE)

  #-----------------------------------------------------------------------------
  # Prepare result
  #-----------------------------------------------------------------------------
  mean1_alt = mle_alt[["mean1"]]
  mean2_alt = mle_alt[["mean2"]]
  ratio_alt = mle_alt[["ratio"]]

  mean1_null = mle_null[["mean1"]]
  mean2_null = mle_null[["mean2"]]

  dispersion1_alt = mle_alt[["dispersion1"]]
  dispersion2_alt = mle_alt[["dispersion2"]]
  dispersion1_null = mle_null[["dispersion1"]]
  dispersion2_null = mle_null[["dispersion2"]]

  n1 <- mle_alt[["n1"]]
  n2 <- mle_alt[["n2"]]

  mle_code <- mle_alt[["mle_code"]]
  mle_message <- mle_alt[["mle_message"]]

  method <- "LRT for independent negative binomial ratio of means"

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
      dispersion1 = dispersion1_alt,
      dispersion2 = dispersion2_alt
    ),
    null = list(
      mean1 = mean1_null,
      mean2 = mean2_null,
      dispersion1 = dispersion1_null,
      dispersion2 = dispersion2_null
    ),
    n1 = n1,
    n2 = n2,
    method = method,
    equal_dispersion = equal_dispersion,
    ratio_null = ratio_null,
    mle_code = mle_code,
    mle_message = mle_message
  )
}
