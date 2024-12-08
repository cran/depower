#' MLE for BNB
#'
#' Maximum likelihood estimates (MLE) for bivariate negative binomial outcomes.
#'
#' These functions are primarily designed for speed in simulation. Missing values
#' are silently excluded.
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
#' The MLEs of \eqn{r} and \eqn{\mu} are \eqn{\hat{r} = \frac{\bar{x}_2}{\bar{x}_1}}
#' and \eqn{\hat{\mu} = \bar{x}_1}. The MLE of \eqn{\theta} is found by
#' maximizing the profile log-likelihood
#' \eqn{l(\hat{r}, \hat{\mu}, \theta)} with respect to \eqn{\theta}. When
#' \eqn{r = r_{null}} is known, the MLE of \eqn{\mu} is
#' \eqn{\tilde{\mu} = \frac{\bar{x}_1 + \bar{x}_2}{1 + r_{null}}} and
#' \eqn{\tilde{\theta}} is obtained by maximizing the profile log-likelihood
#' \eqn{l(r_{null}, \tilde{\mu}, \theta)} with respect to \eqn{\theta}.
#'
#' The backend method for numerical optimization is controlled by argument
#' `method` which refers to [stats::nlm()], [stats::nlminb()], or
#' [stats::optim()]. If you would like to see warnings from the optimizer,
#' include argument `warnings = TRUE`.
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
#'        Typically `ratio_null = 1` (no difference).
#' @param method (string: `"nlm_constrained"`)\cr
#'        The optimization method. Must choose one of `"nlm"`,
#'        `"nlm_constrained"`, `"optim"`, or `"optim_constrained"`. The default
#'        bounds for constrained optimization are `[1e-03, 1e06]`.
#' @param ... Optional arguments passed to the optimization method.
#'
#' @return
#' - For `mle_bnb_alt`, a list with the following elements:
#' \tabular{lll}{
#'   Slot \tab Name \tab Description \cr
#'   1 \tab `mean1` \tab MLE for mean of sample 1. \cr
#'   2 \tab `mean2` \tab MLE for mean of sample 2. \cr
#'   3 \tab `ratio` \tab MLE for ratio of means. \cr
#'   4 \tab `dispersion` \tab MLE for BNB dispersion. \cr
#'   5 \tab `nll` \tab Minimum of negative log-likelihood. \cr
#'   6 \tab `nparams` \tab Number of estimated parameters. \cr
#'   7 \tab `n1` \tab Sample size of sample 1. \cr
#'   8 \tab `n2` \tab Sample size of sample 2. \cr
#'   9 \tab `method` \tab Method used for the results. \cr
#'   10 \tab `mle_method` \tab Method used for optimization. \cr
#'   11 \tab `mle_code` \tab Integer indicating why the optimization process terminated. \cr
#'   12 \tab `mle_message` \tab Additional information from the optimizer.
#' }
#' - For `mle_bnb_null`, a list with the following elements:
#' \tabular{lll}{
#'   Slot \tab Name \tab Description \cr
#'   1 \tab `mean1` \tab MLE for mean of sample 1. \cr
#'   2 \tab `mean2` \tab MLE for mean of sample 2. \cr
#'   3 \tab `ratio_null` \tab Population ratio of means assumed for null hypothesis.
#'                            `mean2 = mean1 * ratio_null`. \cr
#'   4 \tab `dispersion` \tab MLE for BNB dispersion. \cr
#'   5 \tab `nll` \tab Minimum of negative log-likelihood. \cr
#'   6 \tab `nparams` \tab Number of estimated parameters. \cr
#'   7 \tab `n1` \tab Sample size of sample 1. \cr
#'   8 \tab `n2` \tab Sample size of sample 2. \cr
#'   9 \tab `method` \tab Method used for the results. \cr
#'   10 \tab `mle_method` \tab Method used for optimization. \cr
#'   11 \tab `mle_code` \tab Integer indicating why the optimization process terminated. \cr
#'   12 \tab `mle_message` \tab Additional information from the optimizer.
#' }
#'
#' @seealso [depower::sim_bnb()], [depower::nll_bnb]
#'
#' @examples
#' #----------------------------------------------------------------------------
#' # mle_bnb() examples
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
#' mle_alt <- d |>
#'   mle_bnb_alt()
#'
#' mle_null <- d |>
#'   mle_bnb_null()
#'
#' mle_alt
#' mle_null
#'
#' @importFrom stats nlm nlminb optim
#'
#' @name mle_bnb
NULL

#' @export
#' @rdname mle_bnb
mle_bnb_null <- function(
    data,
    ratio_null = 1,
    method = "nlm_constrained",
    ...
) {
  #-----------------------------------------------------------------------------
  # Check args
  #-----------------------------------------------------------------------------
  if(!(is.list(data) && length(data) == 2L)) {
    stop("Argument 'data' must be a list with 2 elements.")
  }
  if(!(length(ratio_null) == 1L && ratio_null > 0)) {
    stop("Argument 'ratio_null' must be a positive scalar numeric.")
  }

  #-----------------------------------------------------------------------------
  # Check arguments
  #-----------------------------------------------------------------------------
  value1 <- data[[1L]]
  value2 <- data[[2L]]

  n1 <- length(value1)
  n2 <- length(value2)

  if(n1 != n2) {
    stop("Argument 'data' must have the same sample size for both samples.")
  }

  if(anyNA(value1) || anyNA(value2)) {
    not_na <- complete.cases(value1, value2)
    value1 <- value1[not_na]
    value2 <- value2[not_na]

    n1 <- length(value1)
    n2 <- length(value2)
  }

  #-----------------------------------------------------------------------------
  # Prepare initial values
  #-----------------------------------------------------------------------------
  init_mean1 <- fmean(value1)

  # dispersion = (mean * p) / (1 - p)
  # Var(X) = mean + mean^2/dispersion
  init_dispersion <- 2

  #-----------------------------------------------------------------------------
  # MLEs
  #-----------------------------------------------------------------------------
  mle <- mle(
    method = method,
    nll = nll_bnb_null,
    parameters = c(mu1 = init_mean1, dispersion = init_dispersion),
    value1 = value1,
    value2 = value2,
    ratio_null = ratio_null,
    ...
  )

  details = "Null hypothesis MLEs for bivariate negative binomial data"

  #-----------------------------------------------------------------------------
  # Return
  #-----------------------------------------------------------------------------
  list(
    mean1 = as.numeric(mle[["estimate"]][1L]),
    mean2 = as.numeric(mle[["estimate"]][1L]) * ratio_null,
    ratio_null = ratio_null,
    dispersion = as.numeric(mle[["estimate"]][2L]),
    n1 = n1,
    n2 = n2,
    nll = mle[["minimum"]],
    nparams = 2L,
    method = details,
    mle_method = method,
    mle_code = mle[["code"]],
    mle_message = mle[["message"]]
  )
}

#' @export
#' @rdname mle_bnb
mle_bnb_alt <- function(data, method = "nlm_constrained", ...) {
  #-----------------------------------------------------------------------------
  # Check args
  #-----------------------------------------------------------------------------
  if(!(is.list(data) && length(data) == 2L)) {
    stop("Argument 'data' must be a list with 2 elements.")
  }

  #-----------------------------------------------------------------------------
  # Check arguments
  #-----------------------------------------------------------------------------
  value1 <- data[[1L]]
  value2 <- data[[2L]]

  n1 <- length(value1)
  n2 <- length(value2)

  if(n1 != n2) {
    stop("Argument 'data' must have the same sample size for both samples.")
  }

  if(anyNA(value1) || anyNA(value2)) {
    not_na <- complete.cases(value1, value2)
    value1 <- value1[not_na]
    value2 <- value2[not_na]

    n1 <- length(value1)
    n2 <- length(value2)
  }

  #-----------------------------------------------------------------------------
  # Prepare initial values
  #-----------------------------------------------------------------------------
  init_mean1 <- fmean(value1)
  init_mean2 <- fmean(value2)

  # dispersion = (mean * p) / (1 - p)
  # Var(X) = mean + mean^2/dispersion
  init_dispersion <- 2

  #-----------------------------------------------------------------------------
  # MLEs
  #-----------------------------------------------------------------------------
  mle <- mle(
    method = method,
    nll = nll_bnb_alt,
    parameters = c(mean1 = init_mean1, mean2 = init_mean2, dispersion = init_dispersion),
    value1 = value1,
    value2 = value2,
    ...
  )

  details = "Alternative hypothesis MLEs for bivariate negative binomial data"

  mean1 = as.numeric(mle[["estimate"]][1L])
  mean2 = as.numeric(mle[["estimate"]][2L])

  #-----------------------------------------------------------------------------
  # Return
  #-----------------------------------------------------------------------------
  list(
    mean1 = mean1,
    mean2 = mean2,
    ratio = mean2 / mean1,
    dispersion = as.numeric(mle[["estimate"]][3L]),
    nll = mle[["minimum"]],
    nparams = 3L,
    n1 = n1,
    n2 = n2,
    method = details,
    mle_method = method,
    mle_code = mle[["code"]],
    mle_message = mle[["message"]]
  )
}
