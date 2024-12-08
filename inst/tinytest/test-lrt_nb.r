library(tinytest)

#-------------------------------------------------------------------------------
# Data
#-------------------------------------------------------------------------------
set.seed(1234)
df <- sim_nb(
  n1 = 30,
  n2 = 30,
  mean1 = 10,
  ratio = 2,
  dispersion1 = 2,
  dispersion2 = 2
)

#-------------------------------------------------------------------------------
# Structure
#-------------------------------------------------------------------------------
expect_equal(
  lrt_nb(df),
  list(
    chisq = 17.6485222363296,
    df = 1L,
    p = 2.65721652816314e-05,
    ratio = 2.8672985640574,
    alternative = list(
      mean1 = 7.03333334185242,
      mean2 = 20.1666665916305,
      dispersion1 = 1.44546942500338,
      dispersion2 = 1.75018905026833
    ),
    null = list(
      mean1 = 15.4987305710827,
      mean2 = 15.4987305710827,
      dispersion1 = 0.823515047968445,
      dispersion2 = 1.56092808992695
    ),
    n1 = 30L,
    n2 = 30L,
    method = "LRT for independent negative binomial ratio of means",
    equal_dispersion = FALSE,
    ratio_null = 1,
    mle_code = 0L,
    mle_message = "relative convergence (4)"
  ),
  tolerance = 0.0001,
  scale = 1
)

#-------------------------------------------------------------------------------
# equal dispersion
#-------------------------------------------------------------------------------
res <- lrt_nb(df, equal_dispersion = TRUE)
expect_equal(res$dispersion1, res$dispersion2)

#-------------------------------------------------------------------------------
# ratio_null
#-------------------------------------------------------------------------------
res_null <- lrt_nb(df)
res_alt <- lrt_nb(df, ratio_null = 1.5)
expect_true(res_null$chisq > res_alt$chisq)
