#' MLE for NB
#'
#' Maximum likelihood estimates (MLE) for two independent negative binomial
#' outcomes.
#'
#' These functions are primarily designed for speed in simulation. Missing values
#' are silently excluded.
#'
#' Suppose \eqn{X_1 \sim \text{NB}(\mu, \theta_1)} and
#' \eqn{X_2 \sim \text{NB}(r\mu, \theta_2)}, where \eqn{X_1} and \eqn{X_2} are
#' independent, \eqn{X_1} is the count outcome for items in group 1, \eqn{X_2}
#' is the count outcome for items in group 2, \eqn{\mu} is the arithmetic mean
#' count in group 1, \eqn{r} is the ratio of arithmetic means for group 2 with
#' respect to group 1, \eqn{\theta_1} is the dispersion parameter of group 1,
#' and \eqn{\theta_2} is the dispersion parameter of group 2.
#'
#' The MLEs of \eqn{r} and \eqn{\mu} are \eqn{\hat{r} = \frac{\bar{x}_2}{\bar{x}_1}}
#' and \eqn{\hat{\mu} = \bar{x}_1}. The MLEs of \eqn{\theta_1} and \eqn{\theta_2}
#' are found by maximizing the profile log-likelihood
#' \eqn{l(\hat{r}, \hat{\mu}, \theta_1, \theta_2)} with respect to
#' \eqn{\theta_1} and \eqn{\theta_2}. When \eqn{r = r_{null}} is known, the MLE
#' of \eqn{\mu} is
#' \eqn{\tilde{\mu} = \frac{n_1 \bar{x}_1 + n_2 \bar{x}_2}{n_1 + n_2}} and
#' \eqn{\tilde{\theta}_1} and \eqn{\tilde{\theta}_2} are obtained by maximizing
#' the profile log-likelihood \eqn{l(r_{null}, \tilde{\mu}, \theta_1, \theta_2)}.
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
#'        from group 1 and the second element is the vector of negative binomial
#'        values from group 2.
#'        \link[base]{NA}s are silently excluded. The default output from
#'        [depower::sim_nb()].
#' @param equal_dispersion (Scalar logical: `FALSE`)\cr
#'        If `TRUE`, the MLEs are calculated assuming both groups have the same
#'        population dispersion parameter. If `FALSE` (default), the MLEs are
#'        calculated assuming different dispersions.
#' @param ratio_null (Scalar numeric: `1`; `(0, Inf)`)\cr
#'        The ratio of means assumed under the null hypothesis (group 2 / group 1).
#'        Typically `ratio_null = 1` (no difference).
#' @param method (string: `"nlm_constrained"`)\cr
#'        The optimization method. Must choose one of `"nlm"`,
#'        `"nlm_constrained"`, `"optim"`, or `"optim_constrained"`. The default
#'        bounds for constrained optimization are `[1e-03, 1e06]`.
#' @param ... Optional arguments passed to the optimization method.
#'
#' @return
#' - For `mle_nb_alt()`, a list with the following elements:
#' \tabular{lll}{
#'   Slot \tab Name \tab Description \cr
#'   1 \tab `mean1` \tab MLE for mean of group 1. \cr
#'   2 \tab `mean2` \tab MLE for mean of group 2. \cr
#'   3 \tab `ratio` \tab MLE for ratio of means. \cr
#'   4 \tab `dispersion1` \tab MLE for dispersion of group 1. \cr
#'   5 \tab `dispersion2` \tab MLE for dispersion of group 2. \cr
#'   6 \tab `equal_dispersion` \tab Were equal dispersions assumed. \cr
#'   7 \tab `n1` \tab Sample size of group 1. \cr
#'   8 \tab `n2` \tab Sample size of group 2. \cr
#'   9 \tab `nll` \tab Minimum of negative log-likelihood. \cr
#'   10 \tab `nparams` \tab Number of estimated parameters. \cr
#'   11 \tab `method` \tab Method used for the results. \cr
#'   12 \tab `mle_method` \tab Method used for optimization. \cr
#'   13 \tab `mle_code` \tab Integer indicating why the optimization process terminated. \cr
#'   14 \tab `mle_message` \tab Additional information from the optimizer.
#' }
#' - For `mle_nb_null()`, a list with the following elements:
#' \tabular{lll}{
#'   Slot \tab Name \tab Description \cr
#'   1 \tab `mean1` \tab MLE for mean of group 1. \cr
#'   2 \tab `mean2` \tab MLE for mean of group 2. \cr
#'   3 \tab `ratio_null` \tab Population ratio of means assumed for null hypothesis.
#'                            `mean2 = mean1 * ratio_null`. \cr
#'   4 \tab `dispersion1` \tab MLE for dispersion of group 1. \cr
#'   5 \tab `dispersion2` \tab MLE for dispersion of group 2. \cr
#'   6 \tab `equal_dispersion` \tab Were equal dispersions assumed. \cr
#'   7 \tab `n1` \tab Sample size of group 1. \cr
#'   8 \tab `n2` \tab Sample size of group 2. \cr
#'   9 \tab `nll` \tab Minimum of negative log-likelihood. \cr
#'   10 \tab `nparams` \tab Number of estimated parameters. \cr
#'   11 \tab `method` \tab Method used for the results. \cr
#'   12 \tab `mle_method` \tab Method used for optimization. \cr
#'   13 \tab `mle_code` \tab Integer indicating why the optimization process terminated. \cr
#'   14 \tab `mle_message` \tab Additional information from the optimizer.
#' }
#'
#' @seealso [depower::sim_nb()], [depower::nll_nb]
#'
#' @examples
#' #----------------------------------------------------------------------------
#' # mle_nb() examples
#' #----------------------------------------------------------------------------
#' library(depower)
#'
#' d <- sim_nb(
#'   n1 = 60,
#'   n2 = 40,
#'   mean1 = 10,
#'   ratio = 1.5,
#'   dispersion1 = 2,
#'   dispersion2 = 8
#' )
#'
#' mle_alt <- d |>
#'   mle_nb_alt()
#'
#' mle_null <- d |>
#'   mle_nb_null()
#'
#' mle_alt
#' mle_null
#'
#' @importFrom stats nlm nlminb optim
#'
#' @name mle_nb
NULL

#' @export
#' @rdname mle_nb
mle_nb_null <- function(
    data,
    equal_dispersion = FALSE,
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
  if(!(is.logical(equal_dispersion) && length(equal_dispersion) == 1L)) {
    stop("Argument 'equal_dispersion' must be a scalar logical.")
  }
  if(!(length(ratio_null) == 1L && ratio_null > 0)) {
    stop("Argument 'ratio_null' must be a positive scalar numeric.")
  }

  #-----------------------------------------------------------------------------
  # Prepare data and initial values
  #-----------------------------------------------------------------------------
  value1 <- data[[1L]]
  value2 <- data[[2L]]

  if(anyNA(value1)) {
    value1 <- value1[!is.na(value1)]
  }
  if(anyNA(value2)) {
    value2 <- value2[!is.na(value2)]
  }

  init_mean1 <- fmean(value1)

  if(equal_dispersion) {
    # dispersion = (mean * p) / (1 - p)
    # Var(X) = mean + mean^2/dispersion
    init_dispersion <- 2
    param = c(mean1 = init_mean1, dispersion = init_dispersion)
    nparams <- 2L
  } else {
    # dispersion = (mean * p) / (1 - p)
    # Var(X) = mean + mean^2/dispersion
    init_dispersion1 <- 2
    init_dispersion2 <- 2
    param = c(mean1 = init_mean1, dispersion1 = init_dispersion1, dispersion2 = init_dispersion2)
    nparams <- 3L
  }

  #-----------------------------------------------------------------------------
  # MLEs
  #-----------------------------------------------------------------------------
  mle <- mle(
    method = method,
    nll = nll_nb_null,
    parameters = param,
    value1 = value1,
    value2 = value2,
    equal_dispersion = equal_dispersion,
    ratio_null = ratio_null,
    ...
  )

  #-----------------------------------------------------------------------------
  # Prepare result
  #-----------------------------------------------------------------------------
  details <- "Null hypothesis MLEs for independent negative binomial data"

  mean1 <- as.numeric(mle[["estimate"]][1L])
  mean2 <- mean1 * ratio_null

  if(equal_dispersion) {
    dispersion1 <- mle[["estimate"]][2L]
    dispersion2 <- dispersion1
  } else {
    dispersion1 <- mle[["estimate"]][2L]
    dispersion2 <- mle[["estimate"]][3L]
  }

  #-----------------------------------------------------------------------------
  # Return
  #-----------------------------------------------------------------------------
  list(
    mean1 = mean1,
    mean2 = mean2,
    ratio_null = ratio_null,
    dispersion1 = as.numeric(dispersion1),
    dispersion2 = as.numeric(dispersion2),
    equal_dispersion = equal_dispersion,
    n1 = length(value1),
    n2 = length(value2),
    nll = mle[["minimum"]],
    nparams = nparams,
    method = details,
    mle_method = method,
    mle_code = mle[["code"]],
    mle_message = mle[["message"]]
  )
}

#' @export
#' @rdname mle_nb
mle_nb_alt <- function(data, equal_dispersion = FALSE, method = "nlm_constrained", ...) {
  #-----------------------------------------------------------------------------
  # Check args
  #-----------------------------------------------------------------------------
  if(!(is.list(data) && length(data) == 2L)) {
    stop("Argument 'data' must be a list with 2 elements.")
  }
  if(!(is.logical(equal_dispersion) && length(equal_dispersion) == 1L)) {
    stop("Argument 'equal_dispersion' must be a scalar logical.")
  }

  #-----------------------------------------------------------------------------
  # Prepare data and initial values
  #-----------------------------------------------------------------------------
  value1 <- data[[1L]]
  value2 <- data[[2L]]

  if(anyNA(value1)) {
    value1 <- value1[!is.na(value1)]
  }
  if(anyNA(value2)) {
    value2 <- value2[!is.na(value2)]
  }

  init_mean1 <- fmean(value1)
  init_mean2 <- fmean(value2)

  if(equal_dispersion) {
    # dispersion = (mean * p) / (1 - p)
    # Var(X) = mean + mean^2/dispersion
    init_dispersion <- 2
    param = c(mean1 = init_mean1, mean2 = init_mean2, dispersion = init_dispersion)
    nparams <- 3L
  } else {
    # dispersion = (mean * p) / (1 - p)
    # Var(X) = mean + mean^2/dispersion
    init_dispersion1 <- 2
    init_dispersion2 <- 2
    param = c(mean1 = init_mean1, mean2 = init_mean2, dispersion1 = init_dispersion1, dispersion2 = init_dispersion2)
    nparams <- 4L
  }

  #-----------------------------------------------------------------------------
  # MLEs
  #-----------------------------------------------------------------------------
  mle <- mle(
    method = method,
    nll = nll_nb_alt,
    parameters = param,
    value1 = value1,
    value2 = value2,
    equal_dispersion = equal_dispersion,
    ...
  )

  #-----------------------------------------------------------------------------
  # Prepare result
  #-----------------------------------------------------------------------------
  details <- "Alternative hypothesis MLEs for independent negative binomial data"

  mean1 = as.numeric(mle[["estimate"]][1L])
  mean2 = as.numeric(mle[["estimate"]][2L])
  ratio = mean2 / mean1

  if(equal_dispersion) {
    dispersion1 <- mle[["estimate"]][3L]
    dispersion2 <- dispersion1
  } else {
    dispersion1 <- mle[["estimate"]][3L]
    dispersion2 <- mle[["estimate"]][4L]
  }

  #-----------------------------------------------------------------------------
  # Return
  #-----------------------------------------------------------------------------
  list(
    mean1 = mean1,
    mean2 = mean2,
    ratio = ratio,
    dispersion1 = as.numeric(dispersion1),
    dispersion2 = as.numeric(dispersion2),
    equal_dispersion = equal_dispersion,
    n1 = length(value1),
    n2 = length(value2),
    nll = mle[["minimum"]],
    nparams = nparams,
    method = details,
    mle_method = method,
    mle_code = mle[["code"]],
    mle_message = mle[["message"]]
  )
}
