library(depower)
library(tinytest)

#-------------------------------------------------------------------------------
# Input validation
#-------------------------------------------------------------------------------
# power must be numeric in (0, 1)
expect_error(
  eval_power_ci(power = 0, nsims = 100),
  pattern = "numeric in \\(0, 1\\)"
)
expect_error(
  eval_power_ci(power = 1, nsims = 100),
  pattern = "numeric in \\(0, 1\\)"
)
expect_error(
  eval_power_ci(power = -0.5, nsims = 100),
  pattern = "numeric in \\(0, 1\\)"
)
expect_error(
  eval_power_ci(power = 1.5, nsims = 100),
  pattern = "numeric in \\(0, 1\\)"
)
expect_error(
  eval_power_ci(power = "0.8", nsims = 100),
  pattern = "numeric in \\(0, 1\\)"
)

# nsims must be integers greater than 1
expect_error(
  eval_power_ci(power = 0.8, nsims = 0),
  pattern = "integers greater than 1"
)
expect_error(
  eval_power_ci(power = 0.8, nsims = -10),
  pattern = "integers greater than 1"
)
expect_error(
  eval_power_ci(power = 0.8, nsims = 10.5),
  pattern = "integers greater than 1"
)

# power and nsims length compatibility
expect_error(
  eval_power_ci(power = c(0.7, 0.8), nsims = c(100, 200, 300)),
  pattern = "same length"
)

# ci_level must be scalar numeric in (0, 1)
expect_error(
  eval_power_ci(power = 0.8, nsims = 100, ci_level = 0),
  pattern = "scalar numeric"
)
expect_error(
  eval_power_ci(power = 0.8, nsims = 100, ci_level = 1),
  pattern = "scalar numeric"
)
expect_error(
  eval_power_ci(power = 0.8, nsims = 100, ci_level = c(0.9, 0.95)),
  pattern = "scalar numeric"
)

#-------------------------------------------------------------------------------
# Output structure
#-------------------------------------------------------------------------------
res <- eval_power_ci(power = 0.8, nsims = 100)

expect_true(is.list(res))
expect_equal(names(res), c("lower", "upper"))
expect_true(is.numeric(res$lower))
expect_true(is.numeric(res$upper))

#-------------------------------------------------------------------------------
# Correct values - matches backend functions
#-------------------------------------------------------------------------------
# Wilson method
res_wilson <- eval_power_ci(power = 0.8, nsims = 100, method = "wilson")
expected_wilson <- depower:::binom_ci_wilson(x = 80, n = 100, conf_level = 0.95)
expect_equal(res_wilson$lower, expected_wilson$lower)
expect_equal(res_wilson$upper, expected_wilson$upper)

# Exact method
res_exact <- eval_power_ci(power = 0.8, nsims = 100, method = "exact")
expected_exact <- depower:::binom_ci_clopper_pearson(
  x = 80,
  n = 100,
  conf_level = 0.95
)
expect_equal(res_exact$lower, expected_exact$lower)
expect_equal(res_exact$upper, expected_exact$upper)

#-------------------------------------------------------------------------------
# Interval properties
#-------------------------------------------------------------------------------
res <- eval_power_ci(power = 0.8, nsims = 100)

# Lower bound should be <= power <= upper bound
expect_true(res$lower <= 0.8)
expect_true(res$upper >= 0.8)

# Bounds should be in [0, 1]
expect_true(res$lower >= 0)
expect_true(res$upper <= 1)

#-------------------------------------------------------------------------------
# Vectorized inputs
#-------------------------------------------------------------------------------
# Vectorized over power
res <- eval_power_ci(power = c(0.7, 0.8, 0.9), nsims = 100)
expect_equal(length(res$lower), 3L)
expect_equal(length(res$upper), 3L)

# Vectorized over nsims
res <- eval_power_ci(power = 0.8, nsims = c(100, 500, 1000))
expect_equal(length(res$lower), 3L)
expect_equal(length(res$upper), 3L)

# Both vectorized (same length)
res <- eval_power_ci(power = c(0.7, 0.8), nsims = c(100, 200))
expect_equal(length(res$lower), 2L)
expect_equal(length(res$upper), 2L)

# Verify first element matches scalar call
res_vec <- eval_power_ci(power = c(0.7, 0.8), nsims = c(100, 100))
res_scalar <- eval_power_ci(power = 0.7, nsims = 100)
expect_equal(res_vec$lower[1], res_scalar$lower)
expect_equal(res_vec$upper[1], res_scalar$upper)

#-------------------------------------------------------------------------------
# More simulations gives narrower interval
#-------------------------------------------------------------------------------
res_small <- eval_power_ci(power = 0.8, nsims = 100)
res_large <- eval_power_ci(power = 0.8, nsims = 1000)

width_small <- res_small$upper - res_small$lower
width_large <- res_large$upper - res_large$lower
expect_true(width_large < width_small)

#-------------------------------------------------------------------------------
# ci_level affects width
#-------------------------------------------------------------------------------
res_95 <- eval_power_ci(power = 0.8, nsims = 100, ci_level = 0.95)
res_99 <- eval_power_ci(power = 0.8, nsims = 100, ci_level = 0.99)

width_95 <- res_95$upper - res_95$lower
width_99 <- res_99$upper - res_99$lower
expect_true(width_99 > width_95)

#-------------------------------------------------------------------------------
# Exact method is wider than Wilson (conservative)
#-------------------------------------------------------------------------------
res_wilson <- eval_power_ci(power = 0.8, nsims = 100, method = "wilson")
res_exact <- eval_power_ci(power = 0.8, nsims = 100, method = "exact")

width_wilson <- res_wilson$upper - res_wilson$lower
width_exact <- res_exact$upper - res_exact$lower
expect_true(width_exact > width_wilson)
