library(depower)
library(tinytest)

# Access internal function
binom_ci_clopper_pearson <- depower:::binom_ci_clopper_pearson

#-------------------------------------------------------------------------------
# Input validation
#-------------------------------------------------------------------------------
# x must be non-negative integers
expect_error(
  binom_ci_clopper_pearson(x = -1, n = 10),
  pattern = "non-negative integers"
)
expect_error(
  binom_ci_clopper_pearson(x = 5.5, n = 10),
  pattern = "non-negative integers"
)
expect_error(
  binom_ci_clopper_pearson(x = "a", n = 10),
  pattern = "non-negative integers"
)

# n must be positive integers
expect_error(
  binom_ci_clopper_pearson(x = 5, n = 0),
  pattern = "positive integers"
)
expect_error(
  binom_ci_clopper_pearson(x = 5, n = -10),
  pattern = "positive integers"
)
expect_error(
  binom_ci_clopper_pearson(x = 5, n = 10.5),
  pattern = "positive integers"
)

# x must be <= n
expect_error(
  binom_ci_clopper_pearson(x = 15, n = 10),
  pattern = "less than or equal"
)

# conf_level must be scalar numeric in (0, 1)
expect_error(
  binom_ci_clopper_pearson(x = 5, n = 10, conf_level = 0),
  pattern = "scalar numeric"
)
expect_error(
  binom_ci_clopper_pearson(x = 5, n = 10, conf_level = 1),
  pattern = "scalar numeric"
)
expect_error(
  binom_ci_clopper_pearson(x = 5, n = 10, conf_level = c(0.9, 0.95)),
  pattern = "scalar numeric"
)
expect_error(
  binom_ci_clopper_pearson(x = 5, n = 10, conf_level = "0.95"),
  pattern = "scalar numeric"
)

#-------------------------------------------------------------------------------
# Output structure
#-------------------------------------------------------------------------------
res <- binom_ci_clopper_pearson(x = 80, n = 100)

expect_true(is.list(res))
expect_equal(names(res), c("lower", "upper"))
expect_true(is.numeric(res$lower))
expect_true(is.numeric(res$upper))

#-------------------------------------------------------------------------------
# Correct values
#-------------------------------------------------------------------------------
# Manual calculation for x = 80, n = 100, conf_level = 0.95
# alpha = 0.05
# lower = qbeta(0.025, 80, 21) = 0.7076394
# upper = qbeta(0.975, 81, 20) = 0.8713065
res <- binom_ci_clopper_pearson(x = 80, n = 100)
expect_equal(res$lower, qbeta(0.025, 80, 21), tolerance = 1e-6)
expect_equal(res$upper, qbeta(0.975, 81, 20), tolerance = 1e-6)

# Edge case: x = 0
res <- binom_ci_clopper_pearson(x = 0, n = 100)
expect_equal(res$lower, 0, tolerance = 1e-6)
expect_true(res$upper > 0)
expect_equal(res$upper, qbeta(0.975, 1, 100), tolerance = 1e-6)

# Edge case: x = n
res <- binom_ci_clopper_pearson(x = 100, n = 100)
expect_true(res$lower < 1)
expect_equal(res$lower, qbeta(0.025, 100, 1), tolerance = 1e-6)
expect_equal(res$upper, 1, tolerance = 1e-6)

#-------------------------------------------------------------------------------
# Vectorized inputs
#-------------------------------------------------------------------------------
res <- binom_ci_clopper_pearson(x = c(8, 80), n = c(10, 100))

expect_equal(length(res$lower), 2L)
expect_equal(length(res$upper), 2L)

# Verify first element matches scalar call
res_scalar <- binom_ci_clopper_pearson(x = 8, n = 10)
expect_equal(res$lower[1], res_scalar$lower)
expect_equal(res$upper[1], res_scalar$upper)

#-------------------------------------------------------------------------------
# Confidence level affects width
#-------------------------------------------------------------------------------
res_95 <- binom_ci_clopper_pearson(x = 50, n = 100, conf_level = 0.95)
res_99 <- binom_ci_clopper_pearson(x = 50, n = 100, conf_level = 0.99)

# 99% CI should be wider than 95% CI
width_95 <- res_95$upper - res_95$lower
width_99 <- res_99$upper - res_99$lower
expect_true(width_99 > width_95)

#-------------------------------------------------------------------------------
# Interval properties
#-------------------------------------------------------------------------------
# Lower bound should be <= point estimate <= upper bound
res <- binom_ci_clopper_pearson(x = 50, n = 100)
p_hat <- 50 / 100
expect_true(res$lower <= p_hat)
expect_true(res$upper >= p_hat)

# Bounds should be in [0, 1]
res <- binom_ci_clopper_pearson(x = c(0, 1, 50, 99, 100), n = rep(100, 5))
expect_true(all(res$lower >= 0))
expect_true(all(res$upper <= 1))

#-------------------------------------------------------------------------------
# Clopper-Pearson is wider than Wilson (conservative)
#-------------------------------------------------------------------------------
cp <- binom_ci_clopper_pearson(x = 50, n = 100)
wilson <- depower:::binom_ci_wilson(x = 50, n = 100)

width_cp <- cp$upper - cp$lower
width_wilson <- wilson$upper - wilson$lower
expect_true(width_cp > width_wilson)
