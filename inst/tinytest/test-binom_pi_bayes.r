library(depower)
library(tinytest)

# Access internal function
binom_pi_bayes <- depower:::binom_pi_bayes

#-------------------------------------------------------------------------------
# Input validation
#-------------------------------------------------------------------------------
# x must be non-negative integers
expect_error(
  binom_pi_bayes(x = -1, n = 10),
  pattern = "non-negative integers"
)
expect_error(
  binom_pi_bayes(x = 5.5, n = 10),
  pattern = "non-negative integers"
)
expect_error(
  binom_pi_bayes(x = "a", n = 10),
  pattern = "non-negative integers"
)

# n must be positive integers
expect_error(
  binom_pi_bayes(x = 5, n = 0),
  pattern = "positive integers"
)
expect_error(
  binom_pi_bayes(x = 5, n = -10),
  pattern = "positive integers"
)
expect_error(
  binom_pi_bayes(x = 5, n = 10.5),
  pattern = "positive integers"
)

# x and n must be same length
expect_error(
  binom_pi_bayes(x = c(5, 6), n = 10),
  pattern = "same length"
)

# x must be <= n
expect_error(
  binom_pi_bayes(x = 15, n = 10),
  pattern = "less than or equal"
)

# future_n must be positive integers or NULL
expect_error(
  binom_pi_bayes(x = 5, n = 10, future_n = 0),
  pattern = "positive integers or NULL"
)
expect_error(
  binom_pi_bayes(x = 5, n = 10, future_n = -10),
  pattern = "positive integers or NULL"
)
expect_error(
  binom_pi_bayes(x = 5, n = 10, future_n = 10.5),
  pattern = "positive integers or NULL"
)

# future_n length must be 1 or same as x
expect_error(
  binom_pi_bayes(x = c(5, 6), n = c(10, 10), future_n = c(20, 30, 40)),
  pattern = "length 1 or the same length"
)

# pred_level must be scalar numeric in (0, 1)
expect_error(
  binom_pi_bayes(x = 5, n = 10, pred_level = 0),
  pattern = "scalar numeric"
)
expect_error(
  binom_pi_bayes(x = 5, n = 10, pred_level = 1),
  pattern = "scalar numeric"
)
expect_error(
  binom_pi_bayes(x = 5, n = 10, pred_level = c(0.9, 0.95)),
  pattern = "scalar numeric"
)
expect_error(
  binom_pi_bayes(x = 5, n = 10, pred_level = "0.95"),
  pattern = "scalar numeric"
)

# prior must be positive numeric vector of length 2
expect_error(
  binom_pi_bayes(x = 5, n = 10, prior = c(1)),
  pattern = "positive numeric vector of length 2"
)
expect_error(
  binom_pi_bayes(x = 5, n = 10, prior = c(1, 0)),
  pattern = "positive numeric vector of length 2"
)
expect_error(
  binom_pi_bayes(x = 5, n = 10, prior = c(-1, 1)),
  pattern = "positive numeric vector of length 2"
)

#-------------------------------------------------------------------------------
# Output structure
#-------------------------------------------------------------------------------
res <- binom_pi_bayes(x = 80, n = 100)

expect_true(is.list(res))
expect_equal(names(res), c("mean", "lower", "upper", "future_n"))
expect_true(is.numeric(res$mean))
expect_true(is.numeric(res$lower))
expect_true(is.numeric(res$upper))
expect_true(is.numeric(res$future_n))

#-------------------------------------------------------------------------------
# Correct values - predictive mean
#-------------------------------------------------------------------------------
# With uniform prior Beta(1,1), posterior is Beta(1+x, 1+n-x)
# Predictive mean = (alpha + x) / (alpha + beta + n) = (1 + 80) / (1 + 1 + 100) = 81/102
res <- binom_pi_bayes(x = 80, n = 100)
expect_equal(res$mean, 81 / 102, tolerance = 1e-6)

# With Jeffreys prior Beta(0.5, 0.5)
# Predictive mean = (0.5 + 80) / (0.5 + 0.5 + 100) = 80.5/101
res <- binom_pi_bayes(x = 80, n = 100, prior = c(0.5, 0.5))
expect_equal(res$mean, 80.5 / 101, tolerance = 1e-6)

#-------------------------------------------------------------------------------
# Interval properties
#-------------------------------------------------------------------------------
# Lower bound should be <= mean <= upper bound
res <- binom_pi_bayes(x = 50, n = 100)
expect_true(res$lower <= res$mean)
expect_true(res$upper >= res$mean)

# Bounds should be in [0, 1]
res <- binom_pi_bayes(x = c(0, 1, 50, 99, 100), n = rep(100, 5))
expect_true(all(res$lower >= 0))
expect_true(all(res$upper <= 1))

#-------------------------------------------------------------------------------
# Edge cases
#-------------------------------------------------------------------------------
# x = 0
res <- binom_pi_bayes(x = 0, n = 100)
expect_equal(res$lower, 0, tolerance = 1e-6)
expect_true(res$upper > 0)
expect_true(res$mean > 0)

# x = n
res <- binom_pi_bayes(x = 100, n = 100)
expect_true(res$lower < 1)
expect_equal(res$upper, 1, tolerance = 1e-6)
expect_true(res$mean < 1)

#-------------------------------------------------------------------------------
# Vectorized inputs
#-------------------------------------------------------------------------------
res <- binom_pi_bayes(x = c(8, 80), n = c(10, 100))

expect_equal(length(res$mean), 2L)
expect_equal(length(res$lower), 2L)
expect_equal(length(res$upper), 2L)

# Verify first element matches scalar call
res_scalar <- binom_pi_bayes(x = 8, n = 10)
expect_equal(res$mean[1], res_scalar$mean)
expect_equal(res$lower[1], res_scalar$lower)
expect_equal(res$upper[1], res_scalar$upper)

#-------------------------------------------------------------------------------
# future_n affects interval width
#-------------------------------------------------------------------------------
# Larger future_n should give narrower interval (more precision)
res_small <- binom_pi_bayes(x = 50, n = 100, future_n = 50)
res_large <- binom_pi_bayes(x = 50, n = 100, future_n = 500)

width_small <- res_small$upper - res_small$lower
width_large <- res_large$upper - res_large$lower
expect_true(width_large < width_small)

# Mean should be unchanged by future_n
expect_equal(res_small$mean, res_large$mean)

#-------------------------------------------------------------------------------
# pred_level affects width
#-------------------------------------------------------------------------------
res_95 <- binom_pi_bayes(x = 50, n = 100, pred_level = 0.95)
res_99 <- binom_pi_bayes(x = 50, n = 100, pred_level = 0.99)

# 99% PI should be wider than 95% PI
width_95 <- res_95$upper - res_95$lower
width_99 <- res_99$upper - res_99$lower
expect_true(width_99 > width_95)

#-------------------------------------------------------------------------------
# Default future_n equals n
#-------------------------------------------------------------------------------
res_default <- binom_pi_bayes(x = 80, n = 100)
res_explicit <- binom_pi_bayes(x = 80, n = 100, future_n = 100)

expect_equal(res_default$mean, res_explicit$mean)
expect_equal(res_default$lower, res_explicit$lower)
expect_equal(res_default$upper, res_explicit$upper)

#-------------------------------------------------------------------------------
# Vectorized future_n
#-------------------------------------------------------------------------------
res <- binom_pi_bayes(x = c(50, 50), n = c(100, 100), future_n = c(50, 200))

# Means should be the same
expect_equal(res$mean[1], res$mean[2])

# Second interval should be narrower (larger future_n)
width1 <- res$upper[1] - res$lower[1]
width2 <- res$upper[2] - res$lower[2]
expect_true(width2 < width1)

# Scalar future_n should be recycled
res_scalar_future <- binom_pi_bayes(
  x = c(50, 50),
  n = c(100, 100),
  future_n = 200
)
expect_equal(res_scalar_future$lower[1], res_scalar_future$lower[2])
