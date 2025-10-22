library(tinytest)
library(depower)

#-------------------------------------------------------------------------------
# method
#-------------------------------------------------------------------------------
data <- sim_nb(
  n1 = 5,
  n2 = 5,
  mean1 = 10,
  ratio = 2,
  dispersion1 = 5,
  dispersion2 = 5,
  nsims = 40
)
res <- power(
  data,
  wald_test_nb(distribution = simulated(nsims = 40)),
  lrt_nb(distribution = simulated(nsims = 40)),
  list_column = TRUE
)
expect_equal(
  current = res$result[[1]][[1]]$method,
  target = "Approximate parametric Wald test for independent negative binomial ratio of means"
)
expect_equal(
  current = res$result[[2]][[1]]$method,
  target = "Approximate parametric LRT for independent negative binomial ratio of means"
)

#-------------------------------------------------------------------------------
# Rettiganti 2012: Equal dispersion figure 2
#-------------------------------------------------------------------------------
# Figure 2
# Exact parametric NB tests for n1=50, n2=50, mean1=5.9,
# dispersion1=0.49, dispersion2=0.49
# Below are figure-derived estimates for ratio=c(0.5, 0.7, 1, 1.2, 1.5, 2)
lrt_power <- c(0.605, 0.205, 0.05, 0.093, 0.27, 0.635)
wald_log_power <- c(0.605, 0.205, 0.05, 0.093, 0.27, 0.635)
wald_identity_power <- c(0.725, 0.315, 0.05, 0.02, 0.03, 0.14)
wald_squared_power <- c(0.735, 0.32, 0.05, 0.011, 0.002, 0.001)
wald_sqrt_power <- c(0.695, 0.275, 0.05, 0.05, 0.16, 0.465)

grid <- expand.grid(
  ratio = c(0.5, 0.7, 1, 1.2, 1.5, 2),
  test = c(
    "lrt_nb(equal_dispersion = TRUE, distribution = simulated(nsims = nsims))",
    "wald_test_nb(equal_dispersion = TRUE, link = \"log\", distribution = simulated(nsims = nsims))",
    "wald_test_nb(equal_dispersion = TRUE, link = \"identity\", distribution = simulated(nsims = nsims))",
    "wald_test_nb(equal_dispersion = TRUE, link = \"squared\", distribution = simulated(nsims = nsims))",
    "wald_test_nb(equal_dispersion = TRUE, link = \"sqrt\", distribution = simulated(nsims = nsims))"
  )
)
grid$expected_power <- c(
  lrt_power,
  wald_log_power,
  wald_identity_power,
  wald_squared_power,
  wald_sqrt_power
)

nsims <- 10000
set.seed(1234)
data <- sim_nb(
  n1 = 50,
  n2 = 50,
  mean1 = 5.9,
  ratio = c(0.5, 0.7, 1, 1.2, 1.5, 2),
  dispersion1 = 0.49,
  dispersion2 = 0.49,
  nsims = nsims
)
power <- power(
  data,
  lrt_nb(
    equal_dispersion = TRUE,
    distribution = simulated(nsims = nsims)
  ),
  wald_test_nb(
    equal_dispersion = TRUE,
    link = "log",
    distribution = simulated(nsims = nsims)
  ),
  wald_test_nb(
    equal_dispersion = TRUE,
    link = "identity",
    distribution = simulated(nsims = nsims)
  ),
  wald_test_nb(
    equal_dispersion = TRUE,
    link = "squared",
    distribution = simulated(nsims = nsims)
  ),
  wald_test_nb(
    equal_dispersion = TRUE,
    link = "sqrt",
    distribution = simulated(nsims = nsims)
  ),
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
# 1. relative difference is <=0.2
# 2. absolute difference is <=0.02
plot(res$abs_diff, res$rel_diff)
which(res$rel_diff > 0.2 & res$abs_diff > 0.02)
expect_true(all(res$rel_diff <= 0.2 | res$abs_diff <= 0.02))

#-------------------------------------------------------------------------------
# Rettiganti 2012: Unequal dispersion figure 3a
#-------------------------------------------------------------------------------
# Figure 3a
# Exact parametric NB tests for n1=50, n2=50, mean1=5.9,
# dispersion1=0.49, dispersion2=0.75*0.49
# Below are figure-derived estimates for ratio=c(0.5, 0.7, 1, 1.2, 1.5, 2)
lrt_75_power <- c(0.54, 0.18, 0.05, 0.085, 0.24, 0.575)
wald_75_log_power <- c(0.55, 0.185, 0.05, 0.085, 0.225, 0.555)
wald_75_identity_power <- c(0.67, 0.28, 0.05, 0.015, 0.01, 0.05)
wald_75_squared_power <- c(0.675, 0.28, 0.05, 0.015, 0.002, 0.0005)

grid <- expand.grid(
  ratio = c(0.5, 0.7, 1, 1.2, 1.5, 2),
  test = c(
    "lrt_nb(distribution = simulated(nsims = nsims))",
    "wald_test_nb(link = \"log\", distribution = simulated(nsims = nsims))",
    "wald_test_nb(link = \"identity\", distribution = simulated(nsims = nsims))",
    "wald_test_nb(link = \"squared\", distribution = simulated(nsims = nsims))"
  )
)
grid$expected_power <- c(
  lrt_75_power,
  wald_75_log_power,
  wald_75_identity_power,
  wald_75_squared_power
)

nsims <- 10000
set.seed(1234)
data <- sim_nb(
  n1 = 50,
  n2 = 50,
  mean1 = 5.9,
  ratio = c(0.5, 0.7, 1, 1.2, 1.5, 2),
  dispersion1 = 0.49,
  dispersion2 = 0.49 * 0.75,
  nsims = nsims
)
power <- power(
  data,
  lrt_nb(distribution = simulated(nsims = nsims)),
  wald_test_nb(link = "log", distribution = simulated(nsims = nsims)),
  wald_test_nb(link = "identity", distribution = simulated(nsims = nsims)),
  wald_test_nb(link = "squared", distribution = simulated(nsims = nsims)),
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

#-------------------------------------------------------------------------------
# Rettiganti 2012: Unequal dispersion figure 3b
#-------------------------------------------------------------------------------
# Figure 3b
# Exact parametric NB tests for n1=50, n2=50, mean1=5.9,
# dispersion1=0.49, dispersion2=1.25*0.49
# Below are figure-derived estimates for ratio=c(0.5, 0.7, 1, 1.2, 1.5, 2)
lrt_125_power <- c(0.655, 0.225, 0.05, 0.1, 0.28, 0.67)
wald_125_log_power <- c(0.65, 0.225, 0.05, 0.1, 0.285, 0.675)
wald_125_identity_power <- c(0.765, 0.335, 0.05, 0.02, 0.05, 0.24)
wald_125_squared_power <- c(0.775, 0.34, 0.05, 0.01, 0.001, 0.0005)

grid <- expand.grid(
  ratio = c(0.5, 0.7, 1, 1.2, 1.5, 2),
  test = c(
    "lrt_nb(distribution = simulated(nsims = nsims))",
    "wald_test_nb(link = \"log\", distribution = simulated(nsims = nsims))",
    "wald_test_nb(link = \"identity\", distribution = simulated(nsims = nsims))",
    "wald_test_nb(link = \"squared\", distribution = simulated(nsims = nsims))"
  )
)
grid$expected_power <- c(
  lrt_125_power,
  wald_125_log_power,
  wald_125_identity_power,
  wald_125_squared_power
)

nsims <- 10000
set.seed(1234)
data <- sim_nb(
  n1 = 50,
  n2 = 50,
  mean1 = 5.9,
  ratio = c(0.5, 0.7, 1, 1.2, 1.5, 2),
  dispersion1 = 0.49,
  dispersion2 = 0.49 * 1.25,
  nsims = nsims
)
power <- power(
  data,
  lrt_nb(distribution = simulated(nsims = nsims)),
  wald_test_nb(link = "log", distribution = simulated(nsims = nsims)),
  wald_test_nb(link = "identity", distribution = simulated(nsims = nsims)),
  wald_test_nb(link = "squared", distribution = simulated(nsims = nsims)),
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
