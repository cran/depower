library(tinytest)

#-------------------------------------------------------------------------------
# Data
#-------------------------------------------------------------------------------
set.seed(1234)
df <- sim_bnb(
  n = 30,
  mean1 = 10,
  ratio = 2,
  dispersion = 2
)

#-------------------------------------------------------------------------------
# Structure
#-------------------------------------------------------------------------------
expect_equal(
  wald_test_bnb(df, ci_level = 0.95),
  list(
    chisq = 102.702520421,
    df = 1L,
    p = 3.89448651999e-24,
    ratio = list(
      estimate = 2.09352020228649,
      lower = 1.81476416139,
      upper = 2.41510592502
    ),
    mean1 = 9.26668400485702,
    mean2 = 19.3999901723732,
    dispersion = 1.3238392543631,
    n1 = 30L,
    n2 = 30L,
    method = "Asymptotic Wald test for bivariate negative binomial ratio of means",
    ci_level = 0.95,
    link = "log",
    ratio_null = 1,
    mle_code = 0L,
    mle_message = "relative convergence (4)"
  ),
  tolerance = 0.001,
  scale = 1
)

#-------------------------------------------------------------------------------
# CI
#-------------------------------------------------------------------------------
res <- wald_test_bnb(df, ci_level = NULL)
expect_true(is.na(res$ratio$lower))
expect_true(is.na(res$ratio$upper))

#-------------------------------------------------------------------------------
# ratio_null
#-------------------------------------------------------------------------------
res_null <- wald_test_bnb(df)
res_alt <- wald_test_bnb(df, ratio_null = 1.5)
expect_true(res_null$chisq > res_alt$chisq)

#-------------------------------------------------------------------------------
# distribution
#-------------------------------------------------------------------------------
set.seed(1234)
df2 <- sim_bnb(
  n = 9,
  mean1 = 10,
  ratio = 2,
  dispersion = 2
)
res_default <- wald_test_bnb(df2)
res_asymptotic <- wald_test_bnb(df2, distribution = asymptotic())
res_approx <- wald_test_bnb(
  df2,
  distribution = simulated(method = "approximate", nsims = 200)
)
res_exact <- wald_test_bnb(df2, distribution = simulated(method = "exact"))

expect_identical(res_default, res_asymptotic)
remove <- grep("^p$|method|distribution", names(res_default))
expect_identical(res_default[-remove], res_approx[-remove])
expect_identical(res_default[-remove], res_exact[-remove])

expect_equal(res_asymptotic$p, 0.0000109102409, tolerance = 0.000001, scale = 1)
expect_equal(res_approx$p, 0.0099, tolerance = 0.0099, scale = 1)
expect_equal(res_exact$p, 0.00389, tolerance = 0.00389, scale = 1)

expect_equal(
  res_asymptotic$method,
  "Asymptotic Wald test for bivariate negative binomial ratio of means"
)
expect_equal(
  res_approx$method,
  "Approximate randomization Wald test for bivariate negative binomial ratio of means"
)
expect_equal(
  res_exact$method,
  "Exact randomization Wald test for bivariate negative binomial ratio of means"
)

#-------------------------------------------------------------------------------
# Variance of the ratio
#-------------------------------------------------------------------------------
# The variance is r * (1 + r) / (n * mean1).
# This value is the element for the ratio in the inverse of the Fisher
# information matrix in Appendix A.2 of Rettiganti and Nagaraja (2012).
# Eq. (10) of that paper is not correct.
x <- c(2L, 4L, 25L, 8L)
y <- c(5L, 9L, 22L, 39L)
res_hand <- wald_test_bnb(list(x, y), ci_level = 0.95)

# These checks use the estimates that the function returns. The calculation has
# no optimizer step, and the tolerance allows only for floating-point rounding.
formula_tolerance <- 1e-12
r_hat <- res_hand$ratio$estimate
sigma_hat <- sqrt(r_hat * (1 + r_hat) / (length(x) * res_hand$mean1))

expected <- list(
  log = (log(r_hat) * r_hat / sigma_hat)^2,
  sqrt = ((sqrt(r_hat) - 1) * 2 * sqrt(r_hat) / sigma_hat)^2,
  identity = ((r_hat - 1) / sigma_hat)^2,
  squared = ((r_hat^2 - 1) / (2 * r_hat * sigma_hat))^2
)
for (link in names(expected)) {
  expect_equal(
    wald_test_bnb(list(x, y), link = link)$chisq,
    expected[[link]],
    tolerance = formula_tolerance,
    info = link
  )
}

expect_equal(
  c(res_hand$ratio$lower, res_hand$ratio$upper),
  exp(log(r_hat) + c(-1, 1) * qnorm(0.975) * sigma_hat / r_hat),
  tolerance = formula_tolerance
)

# These checks use exact totals. nlminb() stops on a change in the objective,
# and it does not bound the error in the means. The tolerance allows for this
# error. It is much smaller than the 10% variance error of Eq. (10) for x, y.
mle_tolerance <- 0.01

# Hand-worked statistic for the log link
expect_equal(
  res_hand$chisq,
  log(sum(y) / sum(x))^2 / (1 / sum(x) + 1 / sum(y)),
  tolerance = mle_tolerance
)

# The statistic depends on the totals and not on the dispersion. These two
# datasets have equal totals and different dispersion estimates.
res_flat <- wald_test_bnb(list(c(9L, 10L, 10L, 10L), c(18L, 19L, 19L, 19L)))
expect_true(res_flat$dispersion > 10 * res_hand$dispersion)
expect_equal(res_flat$chisq, res_hand$chisq, tolerance = mle_tolerance)

# With the log link, a swap of the two samples does not change the statistic.
expect_equal(
  wald_test_bnb(list(y, x))$chisq,
  res_hand$chisq,
  tolerance = mle_tolerance
)

# The statistic equals the Wald statistic of a Poisson GLM with a subject effect.
d_glm <- data.frame(
  count = c(x, y),
  subject = factor(rep(seq_along(x), 2L)),
  sample = rep(0:1, each = length(x))
)
fit_glm <- glm(count ~ subject + sample, family = poisson(), data = d_glm)
expect_equal(
  res_hand$chisq,
  unname(coef(summary(fit_glm))["sample", "z value"]^2),
  tolerance = mle_tolerance
)

# Small means gave a negative variance with Eq. (10), and the p-value was NaN.
res_small <- wald_test_bnb(list(c(1L, 0L, 0L, 0L), c(3L, 2L, 2L, 1L)))
expect_equal(
  res_small$chisq,
  log(8)^2 / (1 + 1 / 8),
  tolerance = mle_tolerance
)
expect_false(is.na(res_small$p))
