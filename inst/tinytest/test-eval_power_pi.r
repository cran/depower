library(depower)
library(tinytest)

#-------------------------------------------------------------------------------
# Input validation
#-------------------------------------------------------------------------------
# power must be numeric in (0, 1)
expect_error(
  eval_power_pi(power = 0, nsims = 100),
  pattern = "numeric in \\(0, 1\\)"
)
expect_error(
  eval_power_pi(power = 1, nsims = 100),
  pattern = "numeric in \\(0, 1\\)"
)
expect_error(
  eval_power_pi(power = -0.5, nsims = 100),
  pattern = "numeric in \\(0, 1\\)"
)
expect_error(
  eval_power_pi(power = 1.5, nsims = 100),
  pattern = "numeric in \\(0, 1\\)"
)
expect_error(
  eval_power_pi(power = "0.8", nsims = 100),
  pattern = "numeric in \\(0, 1\\)"
)

# nsims must be integers greater than 1
expect_error(
  eval_power_pi(power = 0.8, nsims = 0),
  pattern = "integers greater than 1"
)
expect_error(
  eval_power_pi(power = 0.8, nsims = -10),
  pattern = "integers greater than 1"
)
expect_error(
  eval_power_pi(power = 0.8, nsims = 10.5),
  pattern = "integers greater than 1"
)

# power and nsims length compatibility
expect_error(
  eval_power_pi(power = c(0.7, 0.8), nsims = c(100, 200, 300)),
  pattern = "same length"
)

# future_nsims must be integers greater than 1 or NULL
expect_error(
  eval_power_pi(power = 0.8, nsims = 100, future_nsims = 0),
  pattern = "integers greater than 1 or NULL"
)
expect_error(
  eval_power_pi(power = 0.8, nsims = 100, future_nsims = -10),
  pattern = "integers greater than 1 or NULL"
)
expect_error(
  eval_power_pi(power = 0.8, nsims = 100, future_nsims = 10.5),
  pattern = "integers greater than 1 or NULL"
)

# pi_level must be scalar numeric in (0, 1)
expect_error(
  eval_power_pi(power = 0.8, nsims = 100, pi_level = 0),
  pattern = "scalar numeric"
)
expect_error(
  eval_power_pi(power = 0.8, nsims = 100, pi_level = 1),
  pattern = "scalar numeric"
)
expect_error(
  eval_power_pi(power = 0.8, nsims = 100, pi_level = c(0.9, 0.95)),
  pattern = "scalar numeric"
)

# prior must be positive numeric vector of length 2
expect_error(
  eval_power_pi(power = 0.8, nsims = 100, prior = c(1)),
  pattern = "positive numeric vector of length 2"
)
expect_error(
  eval_power_pi(power = 0.8, nsims = 100, prior = c(1, 0)),
  pattern = "positive numeric vector of length 2"
)
expect_error(
  eval_power_pi(power = 0.8, nsims = 100, prior = c(-1, 1)),
  pattern = "positive numeric vector of length 2"
)

#-------------------------------------------------------------------------------
# Output structure
#-------------------------------------------------------------------------------
res <- eval_power_pi(power = 0.8, nsims = 100)

expect_true(is.list(res))
expect_equal(names(res), c("mean", "lower", "upper"))
expect_true(is.numeric(res$mean))
expect_true(is.numeric(res$lower))
expect_true(is.numeric(res$upper))

#-------------------------------------------------------------------------------
# Correct values - matches backend function
#-------------------------------------------------------------------------------
res <- eval_power_pi(power = 0.8, nsims = 100)
expected <- depower:::binom_pi_bayes(x = 80, n = 100, pred_level = 0.95)
expect_equal(res$mean, expected$mean)
expect_equal(res$lower, expected$lower)
expect_equal(res$upper, expected$upper)

# With different prior
res_jeffreys <- eval_power_pi(power = 0.8, nsims = 100, prior = c(0.5, 0.5))
expected_jeffreys <- depower:::binom_pi_bayes(
  x = 80,
  n = 100,
  prior = c(0.5, 0.5)
)
expect_equal(res_jeffreys$mean, expected_jeffreys$mean)
expect_equal(res_jeffreys$lower, expected_jeffreys$lower)
expect_equal(res_jeffreys$upper, expected_jeffreys$upper)

#-------------------------------------------------------------------------------
# Interval properties
#-------------------------------------------------------------------------------
res <- eval_power_pi(power = 0.8, nsims = 100)

# Lower bound should be <= mean <= upper bound
expect_true(res$lower <= res$mean)
expect_true(res$upper >= res$mean)

# Bounds should be in [0, 1]
expect_true(res$lower >= 0)
expect_true(res$upper <= 1)

#-------------------------------------------------------------------------------
# Vectorized inputs
#-------------------------------------------------------------------------------
# Vectorized over power
res <- eval_power_pi(power = c(0.7, 0.8, 0.9), nsims = 100)
expect_equal(length(res$mean), 3L)
expect_equal(length(res$lower), 3L)
expect_equal(length(res$upper), 3L)

# Vectorized over nsims
res <- eval_power_pi(power = 0.8, nsims = c(100, 500, 1000))
expect_equal(length(res$mean), 3L)
expect_equal(length(res$lower), 3L)
expect_equal(length(res$upper), 3L)

# Both vectorized (same length)
res <- eval_power_pi(power = c(0.7, 0.8), nsims = c(100, 200))
expect_equal(length(res$mean), 2L)
expect_equal(length(res$lower), 2L)
expect_equal(length(res$upper), 2L)

# Verify first element matches scalar call
res_vec <- eval_power_pi(power = c(0.7, 0.8), nsims = c(100, 100))
res_scalar <- eval_power_pi(power = 0.7, nsims = 100)
expect_equal(res_vec$mean[1], res_scalar$mean)
expect_equal(res_vec$lower[1], res_scalar$lower)
expect_equal(res_vec$upper[1], res_scalar$upper)

#-------------------------------------------------------------------------------
# More simulations gives narrower interval
#-------------------------------------------------------------------------------
res_small <- eval_power_pi(power = 0.8, nsims = 100)
res_large <- eval_power_pi(power = 0.8, nsims = 1000)

width_small <- res_small$upper - res_small$lower
width_large <- res_large$upper - res_large$lower
expect_true(width_large < width_small)

#-------------------------------------------------------------------------------
# future_nsims affects interval width
#-------------------------------------------------------------------------------
# Larger future_nsims should give narrower interval
res_small_future <- eval_power_pi(power = 0.8, nsims = 100, future_nsims = 50)
res_large_future <- eval_power_pi(power = 0.8, nsims = 100, future_nsims = 500)

width_small_future <- res_small_future$upper - res_small_future$lower
width_large_future <- res_large_future$upper - res_large_future$lower
expect_true(width_large_future < width_small_future)

# Mean should be unchanged by future_nsims
expect_equal(res_small_future$mean, res_large_future$mean)

#-------------------------------------------------------------------------------
# pi_level affects width
#-------------------------------------------------------------------------------
res_95 <- eval_power_pi(power = 0.8, nsims = 100, pi_level = 0.95)
res_99 <- eval_power_pi(power = 0.8, nsims = 100, pi_level = 0.99)

width_95 <- res_95$upper - res_95$lower
width_99 <- res_99$upper - res_99$lower
expect_true(width_99 > width_95)

#-------------------------------------------------------------------------------
# Default future_nsims equals nsims
#-------------------------------------------------------------------------------
res_default <- eval_power_pi(power = 0.8, nsims = 100)
res_explicit <- eval_power_pi(power = 0.8, nsims = 100, future_nsims = 100)

expect_equal(res_default$mean, res_explicit$mean)
expect_equal(res_default$lower, res_explicit$lower)
expect_equal(res_default$upper, res_explicit$upper)
