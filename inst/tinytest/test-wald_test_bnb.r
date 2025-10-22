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
    chisq = 111.168674737738,
    df = 1L,
    p = 5.43444763279816e-26,
    ratio = list(
      estimate = 2.09352020228649,
      lower = 1.824858458574,
      upper = 2.40173522323837
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

expect_equal(res_asymptotic$p, 0.000002103641, tolerance = 0.000001, scale = 1)
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
