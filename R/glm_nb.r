#' GLM for NB ratio of means
#'
#' Generalized linear model for two independent negative binomial outcomes.
#'
#' Uses [glmmTMB::glmmTMB()] in the form
#'
#' ```{r, eval = FALSE}
#' glmmTMB(
#'   formula = value ~ condition,
#'   data = data,
#'   dispformula = ~ condition,
#'   family = nbinom2
#' )
#' ````
#'
#' to model independent negative binomial outcomes
#' \eqn{X_1 \sim \text{NB}(\mu, \theta_1)} and \eqn{X_2 \sim \text{NB}(r\mu, \theta_2)}
#' where \eqn{\mu} is the mean of group 1, \eqn{r} is the ratio of the means of
#' group 2 with respect to group 1, \eqn{\theta_1} is the dispersion parameter
#' of group 1, and \eqn{\theta_2} is the dispersion parameter of group 2.
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
#' arithmetic means for group 2 with respect to group 1 and
#' \eqn{log(r_{null}) = 0} assumes the population means are identical.
#'
#' @references
#' \insertRef{hilbe_2011}{depower}
#'
#' \insertRef{hilbe_2014}{depower}
#'
#' @param data (list)\cr
#'        A list whose first element is the vector of negative binomial values
#'        from group 1 and the second element is the vector of negative binomial
#'        values from group 2.
#'        \link[base]{NA}s are silently excluded. The default output from
#'        [depower::sim_nb()].
#' @param equal_dispersion (Scalar logical: `FALSE`)\cr
#'        If `TRUE`, the model is fit assuming both groups have the same
#'        population dispersion parameter. If `FALSE` (default), the model is
#'        fit assuming different dispersions.
#' @param test (String: `"wald"`; `"c("wald", "lrt")`)\cr
#'        The statistical method used for the test results. `test = "wald"`
#'        performs a Wald test and optionally the Wald confidence intervals.
#'        `test = "lrt"` performs a likelihood ratio test and optionally
#'        the profile likelihood confidence intervals.
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
#'   1 \tab \tab `chisq` \tab Asymptotic \eqn{\chi^2} test statistic for the ratio of means. \cr
#'   2 \tab \tab `df`    \tab Degrees of freedom. \cr
#'   3 \tab \tab `p`     \tab p-value. \cr
#'
#'   4 \tab   \tab `ratio`    \tab Estimated ratio of means (group 2 / group 1). \cr
#'   4 \tab 1 \tab `estimate` \tab Point estimate. \cr
#'   4 \tab 2 \tab `lower`    \tab Confidence interval lower bound. \cr
#'   4 \tab 3 \tab `upper`    \tab Confidence interval upper bound. \cr
#'
#'   5 \tab   \tab `mean1`    \tab Estimated mean of group 1. \cr
#'   5 \tab 1 \tab `estimate` \tab Point estimate. \cr
#'   5 \tab 2 \tab `lower`    \tab Confidence interval lower bound. \cr
#'   5 \tab 3 \tab `upper`    \tab Confidence interval upper bound. \cr
#'
#'   6 \tab   \tab `mean2`    \tab Estimated mean of group 2. \cr
#'   6 \tab 1 \tab `estimate` \tab Point estimate. \cr
#'   6 \tab 2 \tab `lower`    \tab Confidence interval lower bound. \cr
#'   6 \tab 3 \tab `upper`    \tab Confidence interval upper bound. \cr
#'
#'   7 \tab   \tab `dispersion1` \tab Estimated dispersion of group 1. \cr
#'   7 \tab 1 \tab `estimate`    \tab Point estimate. \cr
#'   7 \tab 2 \tab `lower`       \tab Confidence interval lower bound. \cr
#'   7 \tab 3 \tab `upper`       \tab Confidence interval upper bound. \cr
#'
#'   8 \tab   \tab `dispersion2` \tab Estimated dispersion of group 2. \cr
#'   8 \tab 1 \tab `estimate`    \tab Point estimate. \cr
#'   8 \tab 2 \tab `lower`       \tab Confidence interval lower bound. \cr
#'   8 \tab 3 \tab `upper`       \tab Confidence interval upper bound. \cr
#'
#'   9 \tab \tab `n1` \tab Sample size of group 1. \cr
#'   10 \tab \tab `n2` \tab Sample size of group 2. \cr
#'
#'   11 \tab \tab `method`      \tab Method used for the results. \cr
#'   12 \tab \tab `test`        \tab Type of hypothesis test. \cr
#'   13 \tab \tab `alternative` \tab The alternative hypothesis. \cr
#'   14 \tab \tab `equal_dispersion` \tab Whether or not equal dispersions were assumed. \cr
#'   15 \tab \tab `ci_level`    \tab Confidence level of the intervals. \cr
#'   16 \tab \tab `hessian`     \tab Information about the Hessian matrix. \cr
#'   17 \tab \tab `convergence` \tab Information about convergence.
#' }
#'
#' @seealso [glmmTMB::glmmTMB()]
#'
#' @examples
#' #----------------------------------------------------------------------------
#' # glm_nb() examples
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
#' lrt <- glm_nb(d, equal_dispersion = FALSE, test = "lrt", ci_level = 0.95)
#' lrt
#'
#' wald <- glm_nb(d, equal_dispersion = FALSE, test = "wald", ci_level = 0.95)
#' wald
#'
#' #----------------------------------------------------------------------------
#' # Compare results to manual calculation of chi-square statistic
#' #----------------------------------------------------------------------------
#' # Use the same data, but as a data frame instead of list
#' set.seed(1234)
#' d <- sim_nb(
#'   n1 = 60,
#'   n2 = 40,
#'   mean1 = 10,
#'   ratio = 1.5,
#'   dispersion1 = 2,
#'   dispersion2 = 8,
#'   return_type = "data.frame"
#' )
#'
#' mod_alt <- glmmTMB::glmmTMB(
#'   formula = value ~ condition,
#'   data = d,
#'   dispformula = ~ condition,
#'   family = glmmTMB::nbinom2,
#' )
#' mod_null <- glmmTMB::glmmTMB(
#'   formula = value ~ 1,
#'   data = d,
#'   dispformula = ~ condition,
#'   family = glmmTMB::nbinom2,
#' )
#'
#' lrt_chisq <- as.numeric(-2 * (logLik(mod_null) - logLik(mod_alt)))
#' lrt_chisq
#' wald_chisq <- summary(mod_alt)$coefficients$cond["condition2", "z value"]^2
#' wald_chisq
#'
#' anova(mod_null, mod_alt)
#'
#' @importFrom stats logLik pchisq confint
#'
#' @export
glm_nb <- function(
    data,
    equal_dispersion = FALSE,
    test = "wald",
    ci_level = NULL,
    ...
) {
  #-----------------------------------------------------------------------------
  # Check arguments
  #-----------------------------------------------------------------------------
  if(!(is.list(data) && length(data) == 2L)) {
    stop("Argument 'data' must be a list with 2 elements.")
  }

  if(anyNA(data[[1L]]) || anyNA(data[[2L]])) {
    data[[1L]] <- data[[1L]][!is.na(data[[1L]])]
    data[[2L]] <- data[[2L]][!is.na(data[[2L]])]
  }

  n1 <- length(data[[1L]])
  n2 <- length(data[[2L]])
  if(n1 < 2L || n2 < 2L) {
    stop("Argument 'data' must have sample size greater than 1.")
  }

  if(!(is.logical(equal_dispersion) && length(equal_dispersion) == 1L)) {
    stop("Argument 'equal_dispersion' must be a scalar logical.")
  }

  test <- tolower(test)
  if(length(test) != 1L || !test %in% c("wald", "lrt")) {
    stop("Argument 'test' must be a string from 'wald' or 'lrt'.")
  }
  is_lrt <- test == "lrt"

  has_ci <- !is.null(ci_level)
  if(has_ci) {
    if(!is.numeric(ci_level) || length(ci_level) != 1L || ci_level <= 0 || ci_level >= 1) {
      stop("Argument 'ci_level' must be a scalar numeric from (0, 1).")
    }
    ci_method <- if(is_lrt) "profile" else "wald"
  }

  #-----------------------------------------------------------------------------
  # Fit models
  #-----------------------------------------------------------------------------
  data <- list_to_df(data)

  if(equal_dispersion) {
    dispformula <- ~ 1
  } else {
    dispformula <- ~ condition
  }

  mod_alt <- glmmTMB(
    formula = value ~ condition,
    data = data,
    dispformula = dispformula,
    family = nbinom2,
    ...
  )
  if(is_lrt) {
    mod_null <- glmmTMB(
      formula = value ~ 1,
      data = data,
      dispformula = dispformula,
      family = nbinom2,
      ...
    )
  }

  #-----------------------------------------------------------------------------
  # Calculate results
  #-----------------------------------------------------------------------------
  mod_alt_summary <- summary(mod_alt)

  chisq <- if(is_lrt) {
    as.numeric(-2 * (logLik(mod_null) - logLik(mod_alt)))
  } else {
    mod_alt_summary$coefficients$cond["condition2", "z value"]^2
  }
  df <- 1L
  p <- pchisq(chisq, df = df, lower.tail = FALSE)

  mean1 <- exp(mod_alt_summary$coefficients$cond["(Intercept)", "Estimate"])
  mean2 <- exp(
    mod_alt_summary$coefficients$cond["(Intercept)", "Estimate"] +
      mod_alt_summary$coefficients$cond["condition2", "Estimate"]
  )

  ratio <- exp(mod_alt_summary$coefficients$cond["condition2", "Estimate"])

  if(equal_dispersion) {
    dispersion1 <- mod_alt_summary$sigma
    dispersion2 <- mod_alt_summary$sigma
  } else {
    dispersion1 <- exp(mod_alt_summary$coefficients$disp["(Intercept)", "Estimate"])
    dispersion2 <- exp(mod_alt_summary$coefficients$disp["(Intercept)", "Estimate"] +
                         mod_alt_summary$coefficients$disp["condition2", "Estimate"])
  }

  if(has_ci) {
    if(equal_dispersion) {
      beta_ci <- confint(
        object = mod_alt,
        parm = "beta_",
        level = ci_level,
        method = ci_method
      )
      dispersion_ci <- confint(
        object = mod_alt,
        parm = "disp_",
        level = ci_level,
        method = ci_method
      )
      trans <- switch(dimnames(dispersion_ci)[[1]],
        "sigma" = function(x) x,
        "d~(Intercept)" = exp,
        stop("Unable to select tranformation function for dispersion.")
      )
    } else {
      out_ci <- confint(
        object = mod_alt,
        parm = "beta_",
        level = ci_level,
        method = ci_method
      )
      beta_ci <- out_ci[1:2, ]
      # Always 'd' so exp is needed
      dispersion_ci <- out_ci[3:4, ]
    }

    mean1_lower <- exp(beta_ci[1L, 1L])
    mean1_upper <- exp(beta_ci[1L, 2L])

    mean2_lower <- exp(beta_ci[1L, 1L] + beta_ci[2L, 1L])
    mean2_upper <- exp(beta_ci[1L, 2L] + beta_ci[2L, 2L])

    ratio_lower <- exp(beta_ci[2L, 1L])
    ratio_upper <- exp(beta_ci[2L, 2L])

    if(equal_dispersion) {
      dispersion1_lower <- trans(dispersion_ci[1L, 1L])
      dispersion1_upper <- trans(dispersion_ci[1L, 2L])

      dispersion2_lower <- dispersion1_lower
      dispersion2_upper <- dispersion1_upper
    } else {
      dispersion1_lower <- exp(dispersion_ci[1L, 1L])
      dispersion1_upper <- exp(dispersion_ci[1L, 2L])

      dispersion2_lower <- exp(dispersion_ci[1L, 1L] + dispersion_ci[2L, 1L])
      dispersion2_upper <- exp(dispersion_ci[1L, 2L] + dispersion_ci[2L, 2L])
    }
  } else {
    mean1_lower <- NA_real_
    mean1_upper <- NA_real_

    mean2_lower <- NA_real_
    mean2_upper <- NA_real_

    ratio_lower <- NA_real_
    ratio_upper <- NA_real_

    dispersion1_lower <- NA_real_
    dispersion1_upper <- NA_real_

    dispersion2_lower <- NA_real_
    dispersion2_upper <- NA_real_
  }

  hessian <- if(!mod_alt$sdr$pdHess) {
    "Warning: Hessian of fixed effects was not positive definite."
  } else {
    "Hessian appears to be positive definite."
  }
  convergence <- mod_alt$fit$message

  method <- "GLM for independent negative binomial ratio of means"
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
    dispersion1 = list(estimate = dispersion1, lower = dispersion1_lower, upper = dispersion1_upper),
    dispersion2 = list(estimate = dispersion2, lower = dispersion2_lower, upper = dispersion2_upper),
    n1 = n1,
    n2 = n2,
    method = method,
    test = test,
    alternative = alternative,
    equal_dispersion = equal_dispersion,
    ci_level = ci_level,
    hessian = hessian,
    convergence = convergence
  )
}
