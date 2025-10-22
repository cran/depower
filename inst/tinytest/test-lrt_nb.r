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
    method = "Asymptotic LRT for independent negative binomial ratio of means",
    equal_dispersion = FALSE,
    ratio_null = 1,
    mle_code = 0L,
    mle_message = "relative convergence (4)"
  ),
  tolerance = 0.001,
  scale = 1
)

#-------------------------------------------------------------------------------
# equal dispersion
#-------------------------------------------------------------------------------
res <- lrt_nb(df, equal_dispersion = TRUE)
expect_equal(res$null$dispersion1, res$null$dispersion2)

res <- lrt_nb(df, equal_dispersion = TRUE)
expect_equal(res$alt$dispersion1, res$alt$dispersion2)

res <- lrt_nb(df, equal_dispersion = FALSE)
expect_true(res$null$dispersion1 != res$null$dispersion2)

res <- lrt_nb(df, equal_dispersion = FALSE)
expect_true(res$alt$dispersion1 != res$alt$dispersion2)

#-------------------------------------------------------------------------------
# ratio_null
#-------------------------------------------------------------------------------
res_null <- lrt_nb(df)
res_alt <- lrt_nb(df, ratio_null = 1.5)
expect_true(res_null$chisq > res_alt$chisq)

#-------------------------------------------------------------------------------
# distribution
#-------------------------------------------------------------------------------
set.seed(1234)
df2 <- sim_nb(
  n1 = 6,
  n2 = 6,
  mean1 = 10,
  ratio = 2,
  dispersion1 = 2,
  dispersion2 = 2
)
res_default <- lrt_nb(df2)
res_asymptotic <- lrt_nb(df2, distribution = asymptotic())
res_approx <- lrt_nb(
  df2,
  distribution = simulated(method = "approximate", nsims = 200)
)
res_exact <- lrt_nb(df2, distribution = simulated(method = "exact"))

expect_identical(res_default, res_asymptotic)
remove <- grep("^p$|method|distribution", names(res_default))
expect_identical(res_default[-remove], res_approx[-remove])
expect_identical(res_default[-remove], res_exact[-remove])

expect_equal(res_asymptotic$p, 0.004121459, tolerance = 0.0001, scale = 1)
expect_equal(res_approx$p, 0.01485, tolerance = 0.01485, scale = 1)
expect_equal(res_exact$p, 0.00324, tolerance = 0.00324, scale = 1)

expect_equal(
  res_asymptotic$method,
  "Asymptotic LRT for independent negative binomial ratio of means"
)
expect_equal(
  res_approx$method,
  "Approximate randomization LRT for independent negative binomial ratio of means"
)
expect_equal(
  res_exact$method,
  "Exact randomization LRT for independent negative binomial ratio of means"
)
