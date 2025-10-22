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
  lrt_bnb(df),
  list(
    chisq = 109.81838142765,
    df = 1L,
    p = 1.07392901150704e-25,
    ratio = 2.09352020228649,
    alternative = list(
      mean1 = 9.26668400485702,
      mean2 = 19.3999901723732,
      dispersion = 1.3238392543631
    ),
    null = list(
      mean1 = 14.3333374149575,
      mean2 = 14.3333374149575,
      dispersion = 1.32384042150663
    ),
    n1 = 30L,
    n2 = 30L,
    method = "Asymptotic LRT for bivariate negative binomial ratio of means",
    ratio_null = 1,
    mle_code = 0L,
    mle_message = "relative convergence (4)"
  ),
  tolerance = 0.001,
  scale = 1
)

#-------------------------------------------------------------------------------
# ratio_null
#-------------------------------------------------------------------------------
res_null <- lrt_bnb(df)
res_alt <- lrt_bnb(df, ratio_null = 1.5)
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
res_default <- lrt_bnb(df2)
res_asymptotic <- lrt_bnb(df2, distribution = asymptotic())
res_approx <- lrt_bnb(
  df2,
  distribution = simulated(method = "approximate", nsims = 200)
)
res_exact <- lrt_bnb(df2, distribution = simulated(method = "exact"))

expect_identical(res_default, res_asymptotic)
remove <- grep("^p$|method|distribution", names(res_default))
expect_identical(res_default[-remove], res_approx[-remove])
expect_identical(res_default[-remove], res_exact[-remove])

expect_equal(res_asymptotic$p, 0.000005979295, tolerance = 0.000001, scale = 1)
expect_equal(res_approx$p, 0.0099, tolerance = 0.0099, scale = 1)
expect_equal(res_exact$p, 0.003898, tolerance = 0.003898, scale = 1)

expect_equal(
  res_asymptotic$method,
  "Asymptotic LRT for bivariate negative binomial ratio of means"
)
expect_equal(
  res_approx$method,
  "Approximate randomization LRT for bivariate negative binomial ratio of means"
)
expect_equal(
  res_exact$method,
  "Exact randomization LRT for bivariate negative binomial ratio of means"
)
