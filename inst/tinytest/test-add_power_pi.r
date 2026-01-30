library(depower)
library(tinytest)

#-------------------------------------------------------------------------------
# Setup
#-------------------------------------------------------------------------------
set.seed(1234)
sim_data <- sim_nb(
  n1 = 10,
  mean1 = 10,
  ratio = c(1.4, 1.6),
  dispersion1 = 2,
  nsims = 100
)
power_result <- power(sim_data, wald_test_nb())

#-------------------------------------------------------------------------------
# Basic functionality
#-------------------------------------------------------------------------------
res <- add_power_pi(power_result)
expect_true("power_pi_mean" %in% names(res))
#expect_true("power_pi_se" %in% names(res))
expect_true("power_pi_lower" %in% names(res))
expect_true("power_pi_upper" %in% names(res))

# PI bounds are valid
expect_true(all(res$power_pi_lower >= 0))
expect_true(all(res$power_pi_upper <= 1))
expect_true(all(res$power_pi_lower <= res$power_pi_upper))

#-------------------------------------------------------------------------------
# Predictive mean with uniform prior
#-------------------------------------------------------------------------------
# With uniform prior Beta(1,1): E[p_new] = (x + 1) / (n + 2)
n <- power_result$nsims
x <- round(power_result$power * n)
expected_mean <- (x + 1) / (n + 2)
expect_equal(res$power_pi_mean, expected_mean, check.attributes = FALSE)

#-------------------------------------------------------------------------------
# Predictive variance with uniform prior
#-------------------------------------------------------------------------------
# Var[p_new] = (x+1)(n-x+1)(n+2+m) / [m * (n+2)^2 * (n+3)]
# with future_nsims = n (default)
m <- n
alpha_post <- x + 1
beta_post <- n - x + 1
ab_sum <- alpha_post + beta_post
expected_var <- (alpha_post * beta_post * (ab_sum + m)) /
  (m * ab_sum^2 * (ab_sum + 1))
#expect_equal(
#  res$power_pi_se^2,
#  expected_var,
#  tolerance = 1e-10,
#  check.attributes = FALSE
#)

#-------------------------------------------------------------------------------
# Effect of m on prediction interval width
#-------------------------------------------------------------------------------
res_m100 <- add_power_pi(power_result, future_nsims = 100)
res_m500 <- add_power_pi(power_result, future_nsims = 500)
res_m1000 <- add_power_pi(power_result, future_nsims = 1000)

# Larger m should give narrower prediction intervals
width_100 <- res_m100$power_pi_upper - res_m100$power_pi_lower
width_500 <- res_m500$power_pi_upper - res_m500$power_pi_lower
width_1000 <- res_m1000$power_pi_upper - res_m1000$power_pi_lower

expect_true(all(width_500 <= width_100))
expect_true(all(width_1000 <= width_500))

#-------------------------------------------------------------------------------
# Effect of pi_level
#-------------------------------------------------------------------------------
res_90 <- add_power_pi(power_result, pi_level = 0.90)
res_99 <- add_power_pi(power_result, pi_level = 0.99)

# 99% PI should be wider than 95% PI
width_95 <- res$power_pi_upper - res$power_pi_lower
width_99 <- res_99$power_pi_upper - res_99$power_pi_lower
expect_true(all(width_99 >= width_95))

# 90% PI should be narrower than 95% PI
width_90 <- res_90$power_pi_upper - res_90$power_pi_lower
expect_true(all(width_90 <= width_95))

#-------------------------------------------------------------------------------
# Jeffreys prior
#-------------------------------------------------------------------------------
res_jeffreys <- add_power_pi(power_result, prior = c(0.5, 0.5))
expect_true(all(res_jeffreys$power_pi_lower >= 0))
expect_true(all(res_jeffreys$power_pi_upper <= 1))

# Jeffreys prior predictive mean: (x + 0.5) / (n + 1)
expected_mean_jeffreys <- (x + 0.5) / (n + 1)
expect_equal(
  res_jeffreys$power_pi_mean,
  expected_mean_jeffreys,
  check.attributes = FALSE
)

#-------------------------------------------------------------------------------
# Prediction interval vs confidence interval
#-------------------------------------------------------------------------------
# PI should generally be wider than CI (for same level)
res_ci <- add_power_ci(power_result, ci_level = 0.95)
res_pi <- add_power_pi(power_result, pi_level = 0.95)

ci_width <- res_ci$power_ci_upper - res_ci$power_ci_lower
pi_width <- res_pi$power_pi_upper - res_pi$power_pi_lower

# PI is wider (accounts for both parameter uncertainty and sampling variability)
expect_true(all(pi_width >= ci_width - 0.01)) # small tolerance for discretization

#-------------------------------------------------------------------------------
# Edge cases
#-------------------------------------------------------------------------------
# Power = 0
edge_data <- power_result
edge_data$power <- c(0, 0.5)
edge_data$nsims <- c(100, 100)

res_edge <- add_power_pi(edge_data)
expect_equal(res_edge$power_pi_lower[1], 0)
expect_true(res_edge$power_pi_upper[1] > 0)

# Power = 1
edge_data$power <- c(1, 0.5)
res_edge <- add_power_pi(edge_data)
expect_true(res_edge$power_pi_lower[1] < 1)
expect_equal(res_edge$power_pi_upper[1], 1)

#-------------------------------------------------------------------------------
# Argument validation
#-------------------------------------------------------------------------------
expect_error(add_power_pi("not a data frame"))
expect_error(add_power_pi(data.frame(x = 1)))
expect_error(add_power_pi(power_result, future_nsims = 0))
expect_error(add_power_pi(power_result, future_nsims = -1))
expect_error(add_power_pi(power_result, future_nsims = 1.5))
expect_error(add_power_pi(power_result, future_nsims = "100"))
expect_error(add_power_pi(power_result, pi_level = 0))
expect_error(add_power_pi(power_result, pi_level = 1))
expect_error(add_power_pi(power_result, prior = c(1)))
expect_error(add_power_pi(power_result, prior = c(0, 1)))
expect_error(add_power_pi(power_result, prior = c(-1, 1)))
