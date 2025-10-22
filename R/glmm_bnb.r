#' GLMM for BNB ratio of means
#'
#' Generalized linear mixed model for bivariate negative binomial outcomes.
#'
#' Uses [glmmTMB::glmmTMB()] in the form
#'
#' ```{r, eval = FALSE}
#' glmmTMB(
#'   formula = value ~ condition + (1 | item),
#'   data = data,
#'   dispformula = ~ 1,
#'   family = nbinom2
#' )
#' ````
#'
#' to model dependent negative binomial outcomes
#' \eqn{X_1, X_2 \sim \text{BNB}(\mu, r, \theta)} where \eqn{\mu} is the mean of
#' sample 1, \eqn{r} is the ratio of the means of sample 2 with respect to
#' sample 1, and \eqn{\theta} is the dispersion parameter.
#'
#' The hypotheses for the LRT and Wald test of \eqn{r} are
#'
#' \deqn{
#' \begin{aligned}
#' H_{null} &: log(r) = 0 \\
#' H_{alt} &: log(r) \neq 0
#' \end{aligned}
#' }
#'
#' where \eqn{r = \frac{\bar{X}_2}{\bar{X}_1}} is the population ratio of
#' arithmetic means for sample 2 with respect to sample 1 and
#' \eqn{log(r_{null}) = 0} assumes the population means are identical.
#'
#' When simulating data from [depower::sim_bnb()], the mean is a function of the
#' item (subject) random effect which in turn is a function of the dispersion
#' parameter. Thus, `glmm_bnb()` has biased mean and dispersion estimates. The
#' bias increases as the dispersion parameter gets smaller and decreases as
#' the dispersion parameter gets larger. However, estimates of the ratio and
#' standard deviation of the random intercept tend to be accurate. The p-value
#' for `glmm_bnb()` is generally overconservative compared to `glmm_poisson()`,
#' `wald_test_bnb()` and `lrt_bnb()`. In summary, the negative binomial
#' mixed-effects model fit by `glmm_bnb()` is not recommended for the BNB data
#' simulated by `sim_bnb()`. Instead, `wald_test_bnb()` or `lrt_bnb()` should
#' typically be used instead.
#'
#' @references
#' \insertRef{hilbe_2011}{depower}
#'
#' \insertRef{hilbe_2014}{depower}
#'
#' @param data (list)\cr
#'        A list whose first element is the vector of negative binomial values
#'        from sample 1 and the second element is the vector of negative
#'        binomial values from sample 2.
#'        Each vector must be sorted by the subject/item index and must be the
#'        same sample size. \link[base]{NA}s are silently excluded. The default
#'        output from [depower::sim_bnb()].
#' @param test (String: `"wald"`; `"c("wald", "lrt")`)\cr
#'        The statistical method used for the test results. `test = "wald"`
#'        performs a Wald test and optionally the Wald confidence intervals.
#'        `test = "lrt"` performs a likelihood ratio test and optionally
#'        the profile likelihood confidence intervals (means and ratio). The
#'        Wald confidence interval is always used for the limits of the mean of
#'        sample 2, dispersion, and standard deviation of the item (subject)
#'        random intercept.
#' @param ci_level (Scalar numeric: `NULL`; `(0, 1)`)\cr
#'        If `NULL`, confidence intervals are set as `NA`. If in `(0, 1)`,
#'        confidence intervals are calculated at the specified level.
#'        Profile likelihood intervals are computationally intensive, so
#'        intervals from `test = "lrt"` may be slow.
#' @param ... Optional arguments passed to [glmmTMB::glmmTMB()].
#'
#' @return A list with the following elements:
#' \tabular{llll}{
#'   Slot \tab Subslot \tab Name \tab Description \cr
#'
#'   1 \tab \tab `chisq` \tab \eqn{\chi^2} test statistic for the ratio of means. \cr
#'   2 \tab \tab `df`    \tab Degrees of freedom. \cr
#'   3 \tab \tab `p`     \tab p-value. \cr
#'
#'   4 \tab   \tab `ratio`    \tab Estimated ratio of means (sample 2 / sample 1). \cr
#'   4 \tab 1 \tab `estimate` \tab Point estimate. \cr
#'   4 \tab 2 \tab `lower`    \tab Confidence interval lower bound. \cr
#'   4 \tab 3 \tab `upper`    \tab Confidence interval upper bound. \cr
#'
#'   5 \tab   \tab `mean1`    \tab Estimated mean of sample 1. \cr
#'   5 \tab 1 \tab `estimate` \tab Point estimate. \cr
#'   5 \tab 2 \tab `lower`    \tab Confidence interval lower bound. \cr
#'   5 \tab 3 \tab `upper`    \tab Confidence interval upper bound. \cr
#'
#'   6 \tab   \tab `mean2`    \tab Estimated mean of sample 2. \cr
#'   6 \tab 1 \tab `estimate` \tab Point estimate. \cr
#'   6 \tab 2 \tab `lower`    \tab Wald confidence interval lower bound. \cr
#'   6 \tab 3 \tab `upper`    \tab Wald confidence interval upper bound. \cr
#'
#'   7 \tab   \tab `dispersion` \tab Estimated dispersion. \cr
#'   7 \tab 1 \tab `estimate`   \tab Point estimate. \cr
#'   7 \tab 2 \tab `lower`      \tab Confidence interval lower bound. \cr
#'   7 \tab 3 \tab `upper`      \tab Confidence interval upper bound. \cr
#'
#'   8 \tab   \tab `item_sd`  \tab Estimated standard deviation of the item
#'                            (subject) random intercept. \cr
#'   8 \tab 1 \tab `estimate` \tab Point estimate. \cr
#'   8 \tab 2 \tab `lower`    \tab Confidence interval lower bound. \cr
#'   8 \tab 3 \tab `upper`    \tab Confidence interval upper bound. \cr
#'
#'   9 \tab \tab `n1`          \tab Sample size of sample 1. \cr
#'   10 \tab \tab `n2`          \tab Sample size of sample 2. \cr
#'   11 \tab \tab `method`      \tab Method used for the results. \cr
#'   12 \tab \tab `test`        \tab Type of hypothesis test. \cr
#'   13 \tab \tab `alternative` \tab The alternative hypothesis. \cr
#'   14 \tab \tab `ci_level`    \tab Confidence level of the interval. \cr
#'   15 \tab \tab `hessian`     \tab Information about the Hessian matrix. \cr
#'   16 \tab \tab `convergence` \tab Information about convergence.
#' }
#'
#' @seealso [depower::wald_test_bnb()],
#' [depower::lrt_bnb()],
#' [depower::glmm_poisson()]
#'
#' @examples
#' #----------------------------------------------------------------------------
#' # glmm_bnb() examples
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
#' lrt <- glmm_bnb(d, test = "lrt")
#' lrt
#'
#' wald <- glmm_bnb(d, test = "wald", ci_level = 0.95)
#' wald
#'
#' #----------------------------------------------------------------------------
#' # Compare results to manual calculation of chi-square statistic
#' #----------------------------------------------------------------------------
#' # Use the same data, but as a data frame instead of list
#' set.seed(1234)
#' d <- sim_bnb(
#'   n = 40,
#'   mean1 = 10,
#'   ratio = 1.2,
#'   dispersion = 2,
#'   return_type = "data.frame"
#' )
#'
#' mod_alt <- glmmTMB::glmmTMB(
#'   formula = value ~ condition + (1 | item),
#'   data = d,
#'   dispformula = ~ 1,
#'   family = glmmTMB::nbinom2
#' )
#' mod_null <- glmmTMB::glmmTMB(
#'   formula = value ~ 1 + (1 | item),
#'   data = d,
#'   dispformula = ~ 1,
#'   family = glmmTMB::nbinom2
#' )
#'
#' lrt_chisq <- as.numeric(-2 * (logLik(mod_null) - logLik(mod_alt)))
#' lrt_chisq
#' wald_chisq <- summary(mod_alt)$coefficients$cond["condition2", "z value"]^2
#' wald_chisq
#'
#' anova(mod_null, mod_alt)
#'
#' @importFrom glmmTMB glmmTMB nbinom2
#' @importFrom stats logLik pchisq confint vcov qnorm
#'
#' @export
glmm_bnb <- function(data, test = "wald", ci_level = NULL, ...) {
  #-----------------------------------------------------------------------------
  # Check arguments
  #-----------------------------------------------------------------------------
  if (!(is.list(data) && length(data) == 2L)) {
    stop("Argument 'data' must be a list with 2 elements.")
  }

  n1 <- length(data[[1L]])
  n2 <- length(data[[2L]])
  if (n1 != n2) {
    stop("Argument 'data' must have the same sample size for both samples.")
  }
  if (anyNA(data[[1L]]) || anyNA(data[[2L]])) {
    not_na <- complete.cases(data[[1L]], data[[2L]])
    data[[1L]] <- data[[1L]][not_na]
    data[[2L]] <- data[[2L]][not_na]
    n1 <- length(data[[1L]])
    n2 <- length(data[[2L]])
  }
  if (n1 < 2L) {
    stop("Argument 'data' must have sample size greater than 1.")
  }

  test <- tolower(test)
  if (length(test) != 1L || !test %in% c("wald", "lrt")) {
    stop("Argument 'test' must be a string from 'wald' or 'lrt'.")
  }
  is_lrt <- test == "lrt"

  has_ci <- !is.null(ci_level)
  if (has_ci) {
    if (
      !is.numeric(ci_level) ||
        length(ci_level) != 1L ||
        ci_level <= 0 ||
        ci_level >= 1
    ) {
      stop("Argument 'ci_level' must be a scalar numeric from (0, 1).")
    }
    ci_method <- if (is_lrt) "profile" else "wald"
  }

  #-----------------------------------------------------------------------------
  # Fit models
  #-----------------------------------------------------------------------------
  data <- list_to_df(data)

  mod_alt <- glmmTMB(
    formula = value ~ condition + (1 | item),
    data = data,
    dispformula = ~1,
    family = nbinom2,
    ...
  )
  if (is_lrt) {
    mod_null <- glmmTMB(
      formula = value ~ 1 + (1 | item),
      data = data,
      dispformula = ~1,
      family = nbinom2,
      ...
    )
  }

  #-----------------------------------------------------------------------------
  # Calculate results
  #-----------------------------------------------------------------------------
  mod_alt_summary <- summary(mod_alt)

  chisq <- if (is_lrt) {
    as.numeric(-2 * (logLik(mod_null) - logLik(mod_alt)))
  } else {
    mod_alt_summary$coefficients$cond["condition2", "z value"]^2
  }
  df <- 1L
  p <- pchisq(chisq, df = df, lower.tail = FALSE)

  b0 <- mod_alt_summary$coefficients$cond["(Intercept)", "Estimate"]
  b1 <- mod_alt_summary$coefficients$cond["condition2", "Estimate"]
  mean1 <- exp(b0)
  mean2 <- exp(b0 + b1)
  ratio <- exp(b1)

  dispersion <- mod_alt_summary$sigma

  item_sd <- sqrt(as.numeric(mod_alt_summary$varcor$cond$item))

  if (has_ci) {
    beta_ci <- confint(
      object = mod_alt,
      parm = "beta_",
      level = ci_level,
      method = ci_method
    )
    # Currently 'sigma' so no transformation is needed.
    dispersion_ci <- confint(
      object = mod_alt,
      parm = "disp_",
      level = ci_level,
      method = "wald"
    )
    item_sd_ci <- confint(
      object = mod_alt,
      parm = "theta_",
      level = ci_level,
      method = "wald" # Unwanted results for profile...
    )

    mean1_lower <- exp(beta_ci["(Intercept)", 1L])
    mean1_upper <- exp(beta_ci["(Intercept)", 2L])

    # Wald CI for mean2 (profile likelihood not possible?)
    Vcond <- vcov(mod_alt)$cond
    se_m2 <- sqrt(Vcond[1, 1] + Vcond[2, 2] + 2 * Vcond[1, 2])
    mean2_lower <- as.numeric(exp(
      (b0 + b1) + qnorm((1 - ci_level) / 2) * se_m2
    ))
    mean2_upper <- as.numeric(exp(
      (b0 + b1) + qnorm((1 + ci_level) / 2) * se_m2
    ))

    ratio_lower <- exp(beta_ci["condition2", 1L])
    ratio_upper <- exp(beta_ci["condition2", 2L])

    dispersion_lower <- dispersion_ci[1L, 1L]
    dispersion_upper <- dispersion_ci[1L, 2L]

    # Unwanted results for profile...
    item_sd_lower <- item_sd_ci["Std.Dev.(Intercept)|item", 1L]
    item_sd_upper <- item_sd_ci["Std.Dev.(Intercept)|item", 2L]
  } else {
    mean1_lower <- NA_real_
    mean1_upper <- NA_real_

    mean2_lower <- NA_real_
    mean2_upper <- NA_real_

    ratio_lower <- NA_real_
    ratio_upper <- NA_real_

    dispersion_lower <- NA_real_
    dispersion_upper <- NA_real_

    item_sd_lower <- NA_real_
    item_sd_upper <- NA_real_
  }

  hessian <- if (!mod_alt$sdr$pdHess) {
    "Warning: Hessian of fixed effects was not positive definite."
  } else {
    "Hessian appears to be positive definite."
  }
  convergence <- mod_alt$fit$message

  method <- "GLMM for bivariate negative binomial ratio of means"
  alternative <- "two.sided"

  #-----------------------------------------------------------------------------
  # Return
  #-----------------------------------------------------------------------------
  list(
    chisq = chisq,
    df = df,
    p = p,
    ratio = list(estimate = ratio, lower = ratio_lower, upper = ratio_upper),
    mean1 = list(estimate = mean1, lower = mean1_lower, upper = mean1_upper),
    mean2 = list(estimate = mean2, lower = mean2_lower, upper = mean2_upper),
    dispersion = list(
      estimate = dispersion,
      lower = dispersion_lower,
      upper = dispersion_upper
    ),
    item_sd = list(
      estimate = item_sd,
      lower = item_sd_lower,
      upper = item_sd_upper
    ),
    n1 = n1,
    n2 = n2,
    method = method,
    test = test,
    alternative = alternative,
    ci_level = ci_level,
    hessian = hessian,
    convergence = convergence
  )
}
