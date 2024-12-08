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
  wald_test_nb(df, ci_level = 0.95),
  list(
    chisq = 22.878902739329,
    df = 1L,
    p = 1.72535240209498e-06,
    ratio = list(
      estimate = 2.8672985640574,
      lower = 1.86216608006458,
      upper = 4.41496660446127
    ),
    mean1 = 7.03333334185242,
    mean2 = 20.1666665916305,
    dispersion1 = 1.44546942500338,
    dispersion2 = 1.75018905026833,
    n1 = 30L,
    n2 = 30L,
    method = "Wald test for independent negative binomial ratio of means",
    ci_level = 0.95,
    equal_dispersion = FALSE,
    link = "log",
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
res <- wald_test_nb(df, equal_dispersion = TRUE)
expect_equal(res$dispersion1, res$dispersion2)

#-------------------------------------------------------------------------------
# CI
#-------------------------------------------------------------------------------
res <- wald_test_nb(df, ci_level = NULL)
expect_true(is.na(res$ratio$lower))
expect_true(is.na(res$ratio$upper))

#-------------------------------------------------------------------------------
# ratio_null
#-------------------------------------------------------------------------------
res_null <- wald_test_nb(df)
res_alt <- wald_test_nb(df, ratio_null = 1.5)
expect_true(res_null$chisq > res_alt$chisq)
