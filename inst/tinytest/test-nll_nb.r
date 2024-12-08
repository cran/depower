library(tinytest)

#-------------------------------------------------------------------------------
# Unequal disp
#-------------------------------------------------------------------------------
set.seed(1234)
d <- sim_nb(
  n1 = 100,
  n2 = 100,
  mean1 = 10,
  ratio = 1.5,
  dispersion1 = 3,
  dispersion2 = 3
)

alt <- nll_nb_alt(
  param = c(mean1 = 10, mean2 = 15, dispersion = 3),
  value1 = d[[1L]],
  value2 = d[[2L]],
  equal_dispersion = TRUE
)

null <- nll_nb_null(
  param = c(mean1 = 10, dispersion = 3),
  value1 = d[[1L]],
  value2 = d[[2L]],
  equal_dispersion = TRUE,
  ratio_null = 1
)

expect_equal(alt, 674.9149, scale = 1, tolerance = 0.001)
expect_equal(null, 707.5977, scale = 1, tolerance = 0.001)

#-------------------------------------------------------------------------------
# Equal disp
#-------------------------------------------------------------------------------
set.seed(1234)
d <- sim_nb(
  n1 = 100,
  n2 = 100,
  mean1 = 10,
  ratio = 1.5,
  dispersion1 = 3,
  dispersion2 = 100
)

alt <- nll_nb_alt(
  param = c(mean1 = 10, mean2 = 15, dispersion1 = 3, dispersion2 = 100),
  value1 = d[[1L]],
  value2 = d[[2L]],
  equal_dispersion = FALSE
)

null <- nll_nb_null(
  param = c(mean1 = 10, dispersion1 = 3, dispersion2 = 100),
  value1 = d[[1L]],
  value2 = d[[2L]],
  equal_dispersion = FALSE,
  ratio_null = 1
)

expect_equal(alt, 600.1207, scale = 1, tolerance = 0.001)
expect_equal(null, 709.0366, scale = 1, tolerance = 0.001)
