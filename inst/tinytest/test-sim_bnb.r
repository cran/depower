library(tinytest)

#-------------------------------------------------------------------------------
# Structure
#-------------------------------------------------------------------------------
set.seed(1234)
d <- sim_bnb(
  n = 30,
  mean1 = 10,
  ratio = 1.3,
  dispersion = 3,
  nsims = 1
)
expect_equal(length(d), 2L)
expect_true(inherits(d, what = "list"))

set.seed(1234)
d <- sim_bnb(
  n = 30,
  mean1 = 10,
  ratio = 1.3,
  dispersion = 3,
  nsims = 1,
  return_type = "data.frame"
)
expect_equal(length(d), 3L)
expect_true(is.data.frame(d))

set.seed(1234)
d <- sim_bnb(
  n = 30,
  mean1 = 10,
  ratio = 1.3,
  dispersion = 3,
  nsims = 200
)
d$data <- NULL

expect_equal(
  d,
  structure(
    list(
      n1 = structure(30, label = "n1"),
      n2 = structure(30, label = "n2"),
      mean1 = structure(10, label = "Mean1"),
      mean2 = structure(13, label = "Mean2"),
      ratio = structure(1.3, label = "Ratio"),
      dispersion1 = structure(3, label = "Dispersion"),
      nsims = structure(200, label = "N Simulations"),
      distribution = structure("Dependent two-sample BNB", label = "Distribution")),
    row.names = c(NA, -1L),
    class = c("depower", "bnb", "tbl_df", "tbl", "data.frame")
  )
)

set.seed(1234)
d <- sim_bnb(
  n = c(20, 30),
  mean1 = c(10, 15),
  ratio = c(1.3, 1.6),
  dispersion = c(10, 20)
)
d$data <- NULL

expect_equal(
  d,
  structure(
    list(
      n1 = structure(c(20, 30, 20, 30, 20, 30, 20, 30, 20, 30, 20, 30, 20, 30, 20, 30), label = "n1"),
      n2 = structure(c(20, 30, 20, 30, 20, 30, 20, 30, 20, 30, 20, 30, 20, 30, 20, 30), label = "n2"),
      mean1 = structure(c(10, 10, 10, 10, 15, 15, 15, 15, 10, 10, 10, 10, 15, 15, 15, 15), label = "Mean1"),
      mean2 = structure(c(13, 13, 16, 16, 19.5, 19.5, 24, 24, 13, 13, 16, 16, 19.5, 19.5, 24, 24), label = "Mean2"),
      ratio = structure(c(1.3, 1.3, 1.6, 1.6, 1.3, 1.3, 1.6, 1.6, 1.3, 1.3, 1.6, 1.6, 1.3, 1.3, 1.6, 1.6), label = "Ratio"),
      dispersion1 = structure(c(10, 10, 10, 10, 10, 10, 10, 10, 20, 20, 20, 20, 20, 20, 20, 20), label = "Dispersion"),
      nsims = structure(c(1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L), label = "N Simulations"),
      distribution = structure(rep("Dependent two-sample BNB", 16), label = "Distribution")
    ),
    row.names = c(NA, -16L),
    class = c("depower", "bnb", "tbl_df", "tbl", "data.frame")
  )
)

#-------------------------------------------------------------------------------
# Tests
#-------------------------------------------------------------------------------
# distributed as expected
set.seed(1234)
d <- sim_bnb(
  n = 300,
  mean1 = 10,
  ratio = 1.3,
  dispersion = 3,
  nsims = 1
)
expect_equal(mean(d$value1), 10, scale = 1, tolerance = 1)
expect_equal(mean(d$value2), 13, scale = 1, tolerance = 1)
expect_true(var(d$value1) > 10)
expect_true(var(d$value2) > 13)

set.seed(1234)
d <- sim_bnb(
  n = 300,
  mean1 = 10,
  ratio = 1.3,
  dispersion = 1000,
  nsims = 1
)
expect_equal(var(d$value1), 10, scale = 1, tolerance = 0.5)
expect_equal(var(d$value2), 13, scale = 1, tolerance = 0.5)

# zeros are excluded
set.seed(1234)
d <- sim_bnb(
  n = 5,
  mean1 = 10,
  ratio = 1.3,
  dispersion = 0.001,
  nsims = 2,
  max_zeros = 0.99
)
expect_equal(length(d$data[[1]]), 0L)
