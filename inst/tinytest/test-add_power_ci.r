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
# Wilson method (default)
res_wilson <- add_power_ci(power_result)
expect_true("power_ci_lower" %in% names(res_wilson))
expect_true("power_ci_upper" %in% names(res_wilson))

# CI bounds are valid
expect_true(all(res_wilson$power_ci_lower >= 0))
expect_true(all(res_wilson$power_ci_upper <= 1))
expect_true(all(res_wilson$power_ci_lower <= res_wilson$power))
expect_true(all(res_wilson$power_ci_upper >= res_wilson$power))

# Exact method
res_exact <- add_power_ci(power_result, method = "exact")
expect_true(all(res_exact$power_ci_lower >= 0))
expect_true(all(res_exact$power_ci_upper <= 1))

# Exact intervals are wider than Wilson (conservative)
expect_true(all(
  res_exact$power_ci_upper - res_exact$power_ci_lower >=
    res_wilson$power_ci_upper - res_wilson$power_ci_lower - 1e-10
))

#-------------------------------------------------------------------------------
# CI level
#-------------------------------------------------------------------------------
res_90 <- add_power_ci(power_result, ci_level = 0.90)
res_99 <- add_power_ci(power_result, ci_level = 0.99)

# 99% CI should be wider than 95% CI
width_95 <- res_wilson$power_ci_upper - res_wilson$power_ci_lower
width_99 <- res_99$power_ci_upper - res_99$power_ci_lower
expect_true(all(width_99 >= width_95))

# 90% CI should be narrower than 95% CI
width_90 <- res_90$power_ci_upper - res_90$power_ci_lower
expect_true(all(width_90 <= width_95))

#-------------------------------------------------------------------------------
# Edge cases
#-------------------------------------------------------------------------------
# Power = 0
edge_data <- power_result
edge_data$power <- c(0, 0.5)
edge_data$nsims <- c(100, 100)

res_edge <- add_power_ci(edge_data, method = "wilson")
expect_equal(res_edge$power_ci_lower[1], 0, tolerance = 1e-10)
expect_true(res_edge$power_ci_upper[1] > 0)

res_edge_exact <- add_power_ci(edge_data, method = "exact")
expect_equal(res_edge_exact$power_ci_lower[1], 0)

# Power = 1
edge_data$power <- c(1, 0.5)
res_edge <- add_power_ci(edge_data, method = "wilson")
expect_true(res_edge$power_ci_lower[1] < 1)
expect_equal(res_edge$power_ci_upper[1], 1, tolerance = 1e-10)

res_edge_exact <- add_power_ci(edge_data, method = "exact")
expect_equal(res_edge_exact$power_ci_upper[1], 1)

#-------------------------------------------------------------------------------
# Argument validation
#-------------------------------------------------------------------------------
expect_error(add_power_ci("not a data frame"))
expect_error(add_power_ci(data.frame(x = 1)))
expect_error(add_power_ci(power_result, ci_level = 0))
expect_error(add_power_ci(power_result, ci_level = 1))
expect_error(add_power_ci(power_result, ci_level = "0.95"))
expect_error(add_power_ci(power_result, method = "invalid"))
