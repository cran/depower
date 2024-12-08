library(tinytest)

#-------------------------------------------------------------------------------
# Data
#-------------------------------------------------------------------------------
set.seed(1234)
d <- sim_nb(
  n1 = 100,
  n2 = 100,
  mean1 = 10,
  ratio = 1.5,
  dispersion1 = 2,
  dispersion2 = 3
)

#-------------------------------------------------------------------------------
# Structure
#-------------------------------------------------------------------------------
mle_alt <- d |>
  mle_nb_alt()

mle_null <- d |>
  mle_nb_null()

expect_equal(
  mle_alt,
  list(
    mean1 = 9.69999965859556,
    mean2 = 16.040000011519,
    ratio = 1.65360830681116,
    dispersion1 = 1.69578780654153,
    dispersion2 = 3.05237042745438,
    equal_dispersion = FALSE,
    n1 = 100L,
    n2 = 100L,
    nll = 687.403904834001,
    nparams = 4L,
    method = "Alternative hypothesis MLEs for independent negative binomial data",
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
    mean1 = 13.821053815302,
    mean2 = 13.821053815302,
    ratio_null = 1,
    dispersion1 = 1.40703995637252,
    dispersion2 = 2.86301847547128,
    equal_dispersion = FALSE,
    n1 = 100L,
    n2 = 100L,
    nll = 697.909674155722,
    nparams = 3L,
    method = "Null hypothesis MLEs for independent negative binomial data",
    mle_method = "nlm_constrained",
    mle_code = 0L,
    mle_message = "relative convergence (4)"
  ),
  scale = 1,
  tolerance = 0.001
)
