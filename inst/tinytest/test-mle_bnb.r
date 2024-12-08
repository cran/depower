library(tinytest)

#-------------------------------------------------------------------------------
# Data
#-------------------------------------------------------------------------------
set.seed(1234)
d <- sim_bnb(
  n = 100,
  mean1 = 10,
  ratio = 1.5,
  dispersion = 2
)

#-------------------------------------------------------------------------------
# Structure
#-------------------------------------------------------------------------------
mle_alt <- d |>
  mle_bnb_alt()

mle_null <- d |>
  mle_bnb_null()

expect_equal(
  mle_alt,
  list(
    mean1 = 10.7699852853049,
    mean2 = 15.5200160538314,
    ratio = 1.44104338517599,
    dispersion = 2.02796083712345,
    nll = 643.781017465927,
    nparams = 3L,
    n1 = 100L,
    n2 = 100L,
    method = "Alternative hypothesis MLEs for bivariate negative binomial data",
    mle_method = "nlm_constrained",
    mle_code = 0L,
    mle_message = "relative convergence (4)"
  ),
  scale = 1,
  tolerance = 0.001
)

expect_equal(
  mle_null,
  list(
    mean1 = 13.1450001776808,
    mean2 = 13.1450001776808,
    ratio_null = 1,
    dispersion = 2.02796487492354,
    n1 = 100L,
    n2 = 100L,
    nll = 686.928387792057,
    nparams = 2L,
    method = "Null hypothesis MLEs for bivariate negative binomial data",
    mle_method = "nlm_constrained",
    mle_code = 0L,
    mle_message = "relative convergence (4)"
  ),
  scale = 1,
  tolerance = 0.001
)
