library(tinytest)

#-------------------------------------------------------------------------------
# Data
#-------------------------------------------------------------------------------
set.seed(1234)
d <- sim_bnb(
  n = 100,
  mean1 = 10,
  ratio = 1.5,
  dispersion = 3
)

#-------------------------------------------------------------------------------
# Test
#-------------------------------------------------------------------------------
alt <- nll_bnb_alt(
  param = c(mean1 = 10, mean2 = 15, dispersion = 3),
  value1 = d[[1L]],
  value2 = d[[2L]]
)

null <- nll_bnb_null(
  param = c(mean1 = 10, dispersion = 3),
  value1 = d[[1L]],
  value2 = d[[2L]],
  ratio_null = 1
)

expect_equal(alt, 611.1351, scale = 1, tolerance = 0.01)
expect_equal(null, 668.9654, scale = 1, tolerance = 0.01)
