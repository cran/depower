library(tinytest)

#-------------------------------------------------------------------------------
# Structure
#-------------------------------------------------------------------------------
# Independent two sample
set.seed(1234)
d <- sim_log_lognormal(
  n1 = 30,
  n2 = 30,
  ratio = 1.3,
  cv1 = 0.35,
  cv2 = 0.35,
  nsims = 1
)
expect_equal(length(d), 2L)
expect_true(inherits(d, what = "list"))

# one-sample
set.seed(1234)
d <- sim_log_lognormal(
  n1 = 30,
  ratio = 1.3,
  cv1 = 0.35,
  nsims = 1
)
expect_equal(length(d), 1L)
expect_true(inherits(d, what = "list"))


# Independent two sample
set.seed(1234)
d <- sim_log_lognormal(
  n1 = 30,
  n2 = 30,
  ratio = 1.3,
  cv1 = 0.35,
  cv2 = 0.35,
  nsims = 1,
  return_type = "data.frame"
)
expect_equal(length(d), 3L)
expect_true(is.data.frame(d))

# one-sample
set.seed(1234)
d <- sim_log_lognormal(
  n1 = 30,
  ratio = 1.3,
  cv1 = 0.35,
  nsims = 1,
  return_type = "data.frame"
)
expect_equal(length(d), 2L)
expect_true(is.data.frame(d))

# two-sample
set.seed(1234)
d <- sim_log_lognormal(
  n1 = 30,
  n2 = 30,
  ratio = 1.3,
  cv1 = 0.35,
  cv2 = 0.35,
  nsims = 2
)
d$data <- NULL

expect_equal(
  d,
  structure(
    list(
      n1 = structure(30, label = "n1"),
      n2 = structure(30, label = "n2"),
      ratio = structure(1.3, label = "Ratio"),
      cv1 = structure(0.35, label = "CV1"),
      cv2 = structure(0.35, label = "CV2"),
      cor = structure(0, label = "Correlation"),
      nsims = structure(2, label = "N Simulations"),
      distribution = structure(
        "Independent two-sample log(lognormal)",
        label = "Distribution"
      )
    ),
    row.names = c(NA, -1L),
    class = c(
      "depower",
      "log_lognormal_independent_two_sample",
      "tbl_df",
      "tbl",
      "data.frame"
    )
  )
)

# one-sample
set.seed(1234)
d <- sim_log_lognormal(
  n1 = 30,
  ratio = 1.3,
  cv1 = 0.35,
  nsims = 2
)
d$data <- NULL

expect_equal(
  d,
  structure(
    list(
      n1 = structure(30, label = "n Pairs"),
      ratio = structure(1.3, label = "Ratio"),
      cv1 = structure(0.35, label = "CV"),
      nsims = structure(2, label = "N Simulations"),
      distribution = structure(
        "One-sample log(lognormal)",
        label = "Distribution"
      )
    ),
    row.names = c(NA, -1L),
    class = c(
      "depower",
      "log_lognormal_one_sample",
      "tbl_df",
      "tbl",
      "data.frame"
    )
  )
)

# two-sample
set.seed(1234)
d <- sim_log_lognormal(
  n1 = 30,
  n2 = c(30, 40),
  ratio = c(1.3, 2),
  cv1 = 0.35,
  cv2 = 0.35,
  cor = c(0, 0.5),
  messages = FALSE
)
d$data <- NULL

expect_equal(
  d,
  structure(
    list(
      n1 = structure(c(30, 30, 30, 30, 30, 30), label = "n1"),
      n2 = structure(c(30, 40, 30, 40, 30, 30), label = "n2"),
      ratio = structure(c(1.3, 1.3, 2, 2, 1.3, 2), label = "Ratio"),
      cv1 = structure(c(0.35, 0.35, 0.35, 0.35, 0.35, 0.35), label = "CV1"),
      cv2 = structure(c(0.35, 0.35, 0.35, 0.35, 0.35, 0.35), label = "CV2"),
      cor = structure(c(0, 0, 0, 0, 0.5, 0.5), label = "Correlation"),
      nsims = structure(c(1, 1, 1, 1, 1, 1), label = "N Simulations"),
      distribution = structure(
        c(
          rep("Independent two-sample log(lognormal)", 4),
          rep("Dependent two-sample log(lognormal)", 2)
        ),
        label = "Distribution"
      )
    ),
    row.names = c(NA, -6L),
    class = c(
      "depower",
      "log_lognormal_mixed_two_sample",
      "tbl_df",
      "tbl",
      "data.frame"
    )
  )
)

# one-sample
set.seed(1234)
d <- sim_log_lognormal(
  n1 = 30,
  ratio = c(1.3, 2),
  cv1 = c(0.35, 0.5)
)
d$data <- NULL

expect_equal(
  d,
  structure(
    list(
      n1 = structure(c(30, 30, 30, 30), label = "n Pairs"),
      ratio = structure(c(1.3, 2, 1.3, 2), label = "Ratio"),
      cv1 = structure(c(0.35, 0.35, 0.5, 0.5), label = "CV"),
      nsims = structure(c(1L, 1L, 1L, 1L), label = "N Simulations"),
      distribution = structure(
        c(
          "One-sample log(lognormal)",
          "One-sample log(lognormal)",
          "One-sample log(lognormal)",
          "One-sample log(lognormal)"
        ),
        label = "Distribution"
      )
    ),
    row.names = c(NA, -4L),
    class = c(
      "depower",
      "log_lognormal_one_sample",
      "tbl_df",
      "tbl",
      "data.frame"
    )
  )
)

#-------------------------------------------------------------------------------
# Tests
#-------------------------------------------------------------------------------
# distributed as expected
# two-sample
set.seed(1234)
d <- sim_log_lognormal(
  n1 = 300,
  n2 = 300,
  ratio = 1.3,
  cv1 = 0.35,
  cv2 = 0.35,
  nsims = 1
)
expect_equal(mean(d$value1), 0, scale = 1, tolerance = 0.05)
expect_equal(mean(d$value2), log(1.3), scale = 1, tolerance = 0.05)
expect_equal(sqrt(exp(sd(d$value1)^2) - 1), 0.35, scale = 1, tolerance = 0.05)
expect_equal(sqrt(exp(sd(d$value2)^2) - 1), 0.35, scale = 1, tolerance = 0.05)

# one-sample
set.seed(1234)
d <- sim_log_lognormal(
  n1 = 300,
  ratio = 1.3,
  cv1 = 0.35,
  nsims = 1
)
expect_equal(mean(d$value1), log(1.3), scale = 1, tolerance = 0.05)
expect_equal(sqrt(exp(sd(d$value1)^2) - 1), 0.35, scale = 1, tolerance = 0.05)

# Messages
expect_message(
  current = sim_log_lognormal(
    n1 = 30,
    n2 = c(30, 40),
    ratio = c(1.3, 2),
    cv1 = 0.35,
    cv2 = 0.35,
    cor = c(0, 0.5),
    messages = TRUE
  ),
  pattern = "Arguments 'n1' and 'n2' must be the same for correlated data.
2 rows were removed because they had different sample sizes.
6 rows (valid parameter combinations) remain.",
  fixed = TRUE
)
