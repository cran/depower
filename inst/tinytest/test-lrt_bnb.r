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
    method = "LRT for bivariate negative binomial ratio of means",
    ratio_null = 1,
    mle_code = 0L,
    mle_message = "relative convergence (4)"
  ),
  tolerance = 0.0001,
  scale = 1
)

#-------------------------------------------------------------------------------
# ratio_null
#-------------------------------------------------------------------------------
res_null <- lrt_bnb(df)
res_alt <- lrt_bnb(df, ratio_null = 1.5)
expect_true(res_null$chisq > res_alt$chisq)
