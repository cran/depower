library(tinytest)

#-------------------------------------------------------------------------------
# Structure
#-------------------------------------------------------------------------------
set.seed(1234)
d <- sim_nb(
  n1 = 30,
  n2 = 35,
  mean1 = 10,
  ratio = 1.3,
  dispersion1 = 3,
  dispersion2 = 10,
  nsims = 1
)
expect_equal(length(d), 2L)
expect_true(inherits(d, what = "list"))

set.seed(1234)
d <- sim_nb(
  n1 = 30,
  n2 = 35,
  mean1 = 10,
  ratio = 1.3,
  dispersion1 = 3,
  dispersion2 = 10,
  nsims = 1,
  return_type = "data.frame"
)
expect_equal(length(d), 3L)
expect_true(is.data.frame(d))

set.seed(1234)
d <- sim_nb(
  n1 = 30,
  n2 = 35,
  mean1 = 10,
  ratio = 1.3,
  dispersion1 = 3,
  dispersion2 = 10,
  nsims = 2
)
d$data <- NULL

expect_equal(
  d,
  structure(
    list(
      n1 = structure(30, label = "n1"),
      n2 = structure(35, label = "n2"),
      mean1 = structure(10, label = "Mean1"),
      mean2 = structure(13, label = "Mean2"),
      ratio = structure(1.3, label = "Ratio"),
      dispersion1 = structure(3, label = "Dispersion1"),
      dispersion2 = structure(10, label = "Dispersion2"),
      nsims = structure(2, label = "N Simulations"),
      distribution = structure("Independent two-sample NB", label = "Distribution")),
    row.names = c(NA, -1L),
    class = c("depower", "nb", "tbl_df", "tbl", "data.frame")
  )
)

set.seed(1234)
d <- sim_nb(
  n1 = c(30, 40),
  n2 = c(35, 45),
  mean1 = c(10, 15),
  ratio = c(1.3, 1.6),
  dispersion1 = c(3, 5),
  dispersion2 = c(10, 15)
)
d$data <- NULL

expect_equal(
  d,
  structure(
    list(
      n1 = structure(c(30, 40, 30, 40, 30, 40, 30, 40,
                       30, 40, 30, 40, 30, 40, 30, 40, 30, 40, 30, 40, 30, 40, 30, 40,
                       30, 40, 30, 40, 30, 40, 30, 40, 30, 40, 30, 40, 30, 40, 30, 40,
                       30, 40, 30, 40, 30, 40, 30, 40, 30, 40, 30, 40, 30, 40, 30, 40,
                       30, 40, 30, 40, 30, 40, 30, 40), label = "n1"),
      n2 = structure(c(35, 35, 45, 45, 35, 35, 45, 45, 35, 35, 45, 45, 35, 35, 45, 45, 35,
                       35, 45, 45, 35, 35, 45, 45, 35, 35, 45, 45, 35, 35, 45, 45, 35,
                       35, 45, 45, 35, 35, 45, 45, 35, 35, 45, 45, 35, 35, 45, 45, 35,
                       35, 45, 45, 35, 35, 45, 45, 35, 35, 45, 45, 35, 35, 45, 45), label = "n2"),
      mean1 = structure(c(10, 10, 10, 10, 10, 10, 10, 10, 15, 15,
                          15, 15, 15, 15, 15, 15, 10, 10, 10, 10, 10, 10, 10, 10, 15,
                          15, 15, 15, 15, 15, 15, 15, 10, 10, 10, 10, 10, 10, 10, 10,
                          15, 15, 15, 15, 15, 15, 15, 15, 10, 10, 10, 10, 10, 10, 10,
                          10, 15, 15, 15, 15, 15, 15, 15, 15), label = "Mean1"),
      mean2 = structure(c(13, 13, 13, 13, 16, 16, 16, 16, 19.5, 19.5, 19.5, 19.5, 24, 24,
                          24, 24, 13, 13, 13, 13, 16, 16, 16, 16, 19.5, 19.5, 19.5,
                          19.5, 24, 24, 24, 24, 13, 13, 13, 13, 16, 16, 16, 16, 19.5,
                          19.5, 19.5, 19.5, 24, 24, 24, 24, 13, 13, 13, 13, 16, 16,
                          16, 16, 19.5, 19.5, 19.5, 19.5, 24, 24, 24, 24), label = "Mean2"),
      ratio = structure(c(1.3, 1.3, 1.3, 1.3, 1.6, 1.6, 1.6, 1.6,
                          1.3, 1.3, 1.3, 1.3, 1.6, 1.6, 1.6, 1.6, 1.3, 1.3, 1.3, 1.3,
                          1.6, 1.6, 1.6, 1.6, 1.3, 1.3, 1.3, 1.3, 1.6, 1.6, 1.6, 1.6,
                          1.3, 1.3, 1.3, 1.3, 1.6, 1.6, 1.6, 1.6, 1.3, 1.3, 1.3, 1.3,
                          1.6, 1.6, 1.6, 1.6, 1.3, 1.3, 1.3, 1.3, 1.6, 1.6, 1.6, 1.6,
                          1.3, 1.3, 1.3, 1.3, 1.6, 1.6, 1.6, 1.6), label = "Ratio"),
      dispersion1 = structure(c(3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
                                3, 3, 3, 3, 3, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
                                5, 5, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 5,
                                5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5), label = "Dispersion1"),
      dispersion2 = structure(c(10, 10, 10, 10, 10, 10, 10, 10,
                                10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10,
                                10, 10, 10, 10, 10, 10, 10, 10, 10, 15, 15, 15, 15, 15, 15,
                                15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15,
                                15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15), label = "Dispersion2"),
      nsims = structure(c(1L, 1L, 2L, 1L, 1L, 2L, 2L, 2L, 2L, 1L,
                          2L, 1L, 1L, 2L, 2L, 1L, 2L, 1L, 2L, 1L, 2L, 1L, 2L, 2L, 2L,
                          2L, 2L, 2L, 2L, 1L, 2L, 2L, 2L, 2L, 1L, 2L, 2L, 1L, 2L, 2L,
                          2L, 1L, 1L, 2L, 2L, 1L, 2L, 1L, 2L, 2L, 2L, 2L, 2L, 2L, 1L,
                          1L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L), label = "N Simulations"),
      distribution = structure(rep("Independent two-sample NB", 64), label = "Distribution")),
    row.names = c(NA, -64L),
    class = c("depower", "nb", "tbl_df", "tbl", "data.frame")
  )
)

#-------------------------------------------------------------------------------
# Tests
#-------------------------------------------------------------------------------
# distributed as expected
set.seed(1234)
d <- sim_nb(
  n1 = 300,
  n2 = 300,
  mean1 = 10,
  ratio = 1.3,
  dispersion1 = 3,
  dispersion2 = 10,
  nsims = 1
)
expect_equal(mean(d$value1), 10, scale = 1, tolerance = 1)
expect_equal(mean(d$value2), 13, scale = 1, tolerance = 1)
expect_true(var(d$value1) > 10)
expect_true(var(d$value2) > 13)

set.seed(1234)
d <- sim_nb(
  n1 = 300,
  n2 = 300,
  mean1 = 10,
  ratio = 1.3,
  dispersion1 = 1000,
  dispersion2 = 1000,
  nsims = 1
)
expect_equal(var(d$value1), 10, scale = 1, tolerance = 2)
expect_equal(var(d$value2), 13, scale = 1, tolerance = 2)

# zeros are excluded
set.seed(1234)
d <- sim_nb(
  n1 = 10,
  n2 = 10,
  mean1 = 10,
  ratio = 1.3,
  dispersion1 = 0.001,
  dispersion2 = 10,
  nsims = 2,
  max_zeros = 0.99
)
expect_equal(length(d$data[[1]]), 0L)
