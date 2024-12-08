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
    method = "Wald test for bivariate negative binomial ratio of means",
    ci_level = 0.95,
    link = "log",
    ratio_null = 1,
    mle_code = 0L,
    mle_message = "relative convergence (4)"
  ),
  tolerance = 0.0001,
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
