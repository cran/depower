library(tinytest)
library(depower)

#-------------------------------------------------------------------------------
# method
#-------------------------------------------------------------------------------
data <- sim_bnb(
  n = 5,
  mean1 = 10,
  ratio = 2,
  dispersion = 5,
  nsims = 40
)
res <- power(
  data,
  wald_test_bnb(distribution = simulated(nsims = 40)),
  lrt_bnb(distribution = simulated(nsims = 40)),
  list_column = TRUE
)
expect_equal(
  current = res$result[[1]][[1]]$method,
  target = "Approximate parametric Wald test for bivariate negative binomial ratio of means"
)
expect_equal(
  current = res$result[[2]][[1]]$method,
  target = "Approximate parametric LRT for bivariate negative binomial ratio of means"
)

#-------------------------------------------------------------------------------
# Rettiganti 2012
#-------------------------------------------------------------------------------
# Figure 4
# Exact parametric BNB tests for n=10, mean1=5.9, dispersion=0.49
# Figure 5
# Asymptotic BNB tests for n=10, mean1=5.9, dispersion=0.49
# Below are figure-derived estimates for ratio=c(0.5, 0.7, 1, 1.2, 1.5, 2)
lrt_approx_power <- c(0.825, 0.4125, 0.05, 0.175, 0.645, 0.955)
lrt_asympt_power <- lrt_approx_power
rst_approx_power <- c(0.825, 0.4125, 0.05, 0.175, 0.645, 0.955)
wald_approx_log_power <- c(0.825, 0.41, 0.05, 0.175, 0.645, 0.955)
wald_approx_identity_power <- c(0.88, 0.52, 0.05, 0.08, 0.45, 0.885)
wald_approx_squared_power <- c(0.905, 0.54, 0.05, 0.006, 0.08, 0.45)
wald_asympt_squared_power <- c(0.93, 0.625, 0.075, 0.035, 0.265, 0.725)
wald_approx_sqrt_power <- c(0.86, 0.47, 0.05, 0.125, 0.56, 0.94)

grid <- expand.grid(
  ratio = c(0.5, 0.7, 1, 1.2, 1.5, 2),
  test = c(
    "lrt_bnb(distribution = simulated(nsims = nsims))",
    "lrt_bnb(distribution = asymptotic())",
    "wald_test_bnb(link = \"log\", distribution = simulated(nsims = nsims))",
    "wald_test_bnb(link = \"identity\", distribution = simulated(nsims = nsims))",
    "wald_test_bnb(link = \"squared\", distribution = simulated(nsims = nsims))",
    "wald_test_bnb(link = \"squared\", distribution = asymptotic())",
    "wald_test_bnb(link = \"sqrt\", distribution = simulated(nsims = nsims))"
  )
)
grid$expected_power <- c(
  lrt_approx_power,
  lrt_asympt_power,
  wald_approx_log_power,
  wald_approx_identity_power,
  wald_approx_squared_power,
  wald_asympt_squared_power,
  wald_approx_sqrt_power
)

nsims <- 20000
set.seed(1234)
data <- sim_bnb(
  n = 10,
  mean1 = 5.9,
  ratio = c(0.5, 0.7, 1, 1.2, 1.5, 2),
  dispersion = 0.49,
  nsims = nsims
)
power <- power(
  data,
  lrt_bnb(distribution = simulated(nsims = nsims)),
  lrt_bnb(distribution = asymptotic()),
  wald_test_bnb(link = "log", distribution = simulated(nsims = nsims)),
  wald_test_bnb(link = "identity", distribution = simulated(nsims = nsims)),
  wald_test_bnb(link = "squared", distribution = simulated(nsims = nsims)),
  wald_test_bnb(link = "squared", distribution = asymptotic()),
  wald_test_bnb(link = "sqrt", distribution = simulated(nsims = nsims)),
  alpha = 0.05,
  ncores = 4
)

# should look like Figure 4 and 5 of Rettiganti 2012
plot(power)

res <- grid |>
  dplyr::left_join(power[c("ratio", "test", "power")]) |>
  dplyr::mutate(
    abs_diff = abs(power - expected_power),
    rel_diff = abs(power - expected_power) / expected_power
  )

# At least one of the following must be true:
# 1. relative difference is <=0.1
# 2. absolute difference is <=0.02
plot(res$abs_diff, res$rel_diff)
which(res$rel_diff > 0.1 & res$abs_diff > 0.02)
expect_true(all(res$rel_diff <= 0.1 | res$abs_diff <= 0.02))
