library(tinytest)

#-------------------------------------------------------------------------------
# Structure
#-------------------------------------------------------------------------------
# Independent two-sample t-Test
set.seed(1234)
d <- sim_log_lognormal(
  n1 = 20,
  n2 = 20,
  ratio = c(1.2, 1.4),
  cv1 = 0.4,
  cv2 = 0.4,
  cor = 0,
  nsims = 1000
) |>
  power()

expect_equal(
  d,
  structure(list(n1 = structure(c(20, 20), label = "n1"), n2 = structure(c(20,
20), label = "n2"), ratio = structure(c(1.2, 1.4), label = "Ratio"),
    cv1 = structure(c(0.4, 0.4), label = "CV1"), cv2 = structure(c(0.4,
    0.4), label = "CV2"), cor = structure(c(0, 0), label = "Correlation"),
    distribution = structure(c("Independent two-sample log(lognormal)",
    "Independent two-sample log(lognormal)"), label = "Distribution"),
    nsims = structure(c(1000, 1000), label = "N Simulations"),
    test = structure(c(`t_test_welch()` = "Welch's t-Test", `t_test_welch()` = "Welch's t-Test"
    ), label = "Test"), alpha = structure(c(0.05, 0.05), label = "Alpha"),
    power = structure(c(0.312, 0.772), label = "Power")), row.names = c(NA,
-2L), class = c("depower", "log_lognormal_independent_two_sample",
"tbl_df", "tbl", "data.frame"))
)

# Dependent two-sample t-Test
set.seed(1234)
d <- sim_log_lognormal(
  n1 = 20,
  n2 = 20,
  ratio = c(1.2, 1.4),
  cv1 = 0.4,
  cv2 = 0.4,
  cor = 0.5,
  nsims = 1000
) |>
  power()

expect_equal(
  d,
  structure(list(n1 = structure(c(20, 20), label = "n1"), n2 = structure(c(20,
20), label = "n2"), ratio = structure(c(1.2, 1.4), label = "Ratio"),
    cv1 = structure(c(0.4, 0.4), label = "CV1"), cv2 = structure(c(0.4,
    0.4), label = "CV2"), cor = structure(c(0.5, 0.5), label = "Correlation"),
    distribution = structure(c("Dependent two-sample log(lognormal)",
    "Dependent two-sample log(lognormal)"), label = "Distribution"),
    nsims = structure(c(1000, 1000), label = "N Simulations"),
    test = structure(c(`t_test_paired()` = "Paired t-Test", `t_test_paired()` = "Paired t-Test"
    ), label = "Test"), alpha = structure(c(0.05, 0.05), label = "Alpha"),
    power = structure(c(0.535, 0.957), label = "Power")), row.names = c(NA,
-2L), class = c("depower", "log_lognormal_dependent_two_sample",
"tbl_df", "tbl", "data.frame"))
)

# Mixed-type two-sample t-Test
set.seed(1234)
d <- sim_log_lognormal(
  n1 = 20,
  n2 = 20,
  ratio = c(1.2, 1.4),
  cv1 = 0.4,
  cv2 = 0.4,
  cor = c(0, 0.5),
  nsims = 1000
) |>
  power()

expect_equal(
  d,
  structure(list(n1 = structure(c(20, 20, 20, 20), label = "n1"),
    n2 = structure(c(20, 20, 20, 20), label = "n2"), ratio = structure(c(1.2,
    1.4, 1.2, 1.4), label = "Ratio"), cv1 = structure(c(0.4,
    0.4, 0.4, 0.4), label = "CV1"), cv2 = structure(c(0.4, 0.4,
    0.4, 0.4), label = "CV2"), cor = structure(c(0, 0, 0.5, 0.5
    ), label = "Correlation"), distribution = structure(c("Independent two-sample log(lognormal)",
    "Independent two-sample log(lognormal)", "Dependent two-sample log(lognormal)",
    "Dependent two-sample log(lognormal)"), label = "Distribution"),
    nsims = structure(c(1000, 1000, 1000, 1000), label = "N Simulations"),
    test = structure(c(`t_test_welch()` = "Welch's t-Test", `t_test_welch()` = "Welch's t-Test",
    `t_test_paired()` = "Paired t-Test", `t_test_paired()` = "Paired t-Test"
    ), label = "Test"), alpha = structure(c(0.05, 0.05, 0.05,
    0.05), label = "Alpha"), power = structure(c(0.312, 0.772,
    0.54, 0.963), label = "Power")), row.names = c(NA, -4L), class = c("depower",
"log_lognormal_mixed_two_sample", "tbl_df", "tbl", "data.frame"))
)

# One-sample t-Test
set.seed(1234)
d <- sim_log_lognormal(
  n1 = 20,
  ratio = c(1.2, 1.4),
  cv1 = 0.4,
  nsims = 1000
) |>
  power()

expect_equal(
  d,
  structure(list(n1 = structure(c(20, 20), label = "n1"), ratio = structure(c(1.2,
1.4), label = "Ratio"), cv1 = structure(c(0.4, 0.4), label = "CV1"),
    distribution = structure(c("One-sample log(lognormal)", "One-sample log(lognormal)"
    ), label = "Distribution"), nsims = structure(c(1000, 1000
    ), label = "N Simulations"), test = structure(c(`t_test_paired()` = "One-sample t-Test",
    `t_test_paired()` = "One-sample t-Test"), label = "Test"),
    alpha = structure(c(0.05, 0.05), label = "Alpha"), power = structure(c(0.515,
    0.958), label = "Power")), row.names = c(NA, -2L), class = c("depower",
"log_lognormal_one_sample", "tbl_df", "tbl", "data.frame"))
)

# NB test
set.seed(1234)
d <- sim_nb(
  n1 = 10,
  mean1 = 10,
  ratio = c(1.6, 2),
  dispersion1 = 2,
  dispersion2 = 2,
  nsims = 200
) |>
  power()

expect_equal(
  d,
  structure(list(n1 = structure(c(10, 10), label = "n1"), n2 = structure(c(10,
10), label = "n2"), mean1 = structure(c(10, 10), label = "Mean1"),
    mean2 = structure(c(16, 20), label = "Mean2"), ratio = structure(c(1.6,
    2), label = "Ratio"), dispersion1 = structure(c(2,
    2), label = "Dispersion1"), dispersion2 = structure(c(2,
    2), label = "Dispersion2"), distribution = structure(c("Independent two-sample NB",
    "Independent two-sample NB"), label = "Distribution"), nsims = structure(c(200,
    200), label = "N Simulations"), test = structure(c(`wald_test_nb()` = "NB Wald test",
    `wald_test_nb()` = "NB Wald test"), label = "Test"), alpha = structure(c(0.05,
    0.05), label = "Alpha"), power = structure(c(0.315, 0.545
    ), label = "Power")), row.names = c(NA, -2L), class = c("depower",
"nb", "tbl_df", "tbl", "data.frame"))
)

# BNB test
set.seed(1234)
d <- sim_bnb(
  n = 10,
  mean1 = 10,
  ratio = c(1.4, 1.6),
  dispersion = 10,
  nsims = 200
) |>
  power()

expect_equal(
  d,
  structure(list(n1 = structure(c(10, 10), label = "n1"), n2 = structure(c(10,
10), label = "n2"), mean1 = structure(c(10, 10), label = "Mean1"),
    mean2 = structure(c(14, 16), label = "Mean2"), ratio = structure(c(1.4,
    1.6), label = "Ratio"), dispersion1 = structure(c(10,
    10), label = "Dispersion1"), distribution = structure(c("Dependent two-sample BNB",
    "Dependent two-sample BNB"), label = "Distribution"), nsims = structure(c(200,
    200), label = "N Simulations"), test = structure(c(`wald_test_bnb()` = "BNB Wald test",
    `wald_test_bnb()` = "BNB Wald test"), label = "Test"), alpha = structure(c(0.05,
    0.05), label = "Alpha"), power = structure(c(0.86, 0.97), label = "Power")), row.names = c(NA,
-2L), class = c("depower", "bnb", "tbl_df", "tbl", "data.frame"
))
)

#-------------------------------------------------------------------------------
# Validated t-test power
# from SAS
#-------------------------------------------------------------------------------
set.seed(1234)
# Independent two sample
expect_equal(
  sim_log_lognormal(
    n1 = 30,
    n2 = 30,
    ratio = 1.3,
    cv1 = 0.35,
    cv2 = 0.35,
    nsims = 10000
  ) |>
    power() |>
    getElement("power") |>
    as.numeric(),
  0.8363,
  tolerance = 0.03,
  scale = 1
)

# Independent two sample (inverted fc)
expect_equal(
  sim_log_lognormal(
    n1 = 30,
    n2 = 30,
    ratio = 1/1.3,
    cv1 = 0.35,
    cv2 = 0.35,
    nsims = 10000
  ) |>
    power() |>
    getElement("power") |>
    as.numeric(),
  0.8363,
  tolerance = 0.03,
  scale = 1
)

# Independent two sample (different sample sizes)
expect_equal(
  sim_log_lognormal(
    n1 = 30,
    n2 = 100,
    ratio = 1.3,
    cv1 = 0.35,
    cv2 = 0.35,
    nsims = 10000
  ) |>
    power() |>
    getElement("power") |>
    as.numeric(),
  0.9573,
  tolerance = 0.02,
  scale = 1
)

# Independent two sample (different CV)
expect_equal(
  sim_log_lognormal(
    n1 = 30,
    n2 = 30,
    ratio = 1.1,
    cv1 = 0.3,
    cv2 = 0.5,
    cor = 0.345,
    nsims = 10000
  ) |>
    power() |>
    getElement("power") |>
    as.numeric(),
  0.19762,
  tolerance = 0.04,
  scale = 1
)

# Paired two sample (Same CV)
expect_equal(
  sim_log_lognormal(
    n1 = 30,
    n2 = 30,
    ratio = 1.1,
    cv1 = 0.3,
    cv2 = 0.3,
    cor = 0.345,
    nsims = 10000
  ) |>
    power() |>
    getElement("power") |>
    as.numeric(),
  0.32798,
  tolerance = 0.04,
  scale = 1
)

# Paired two sample (Different CV)
expect_equal(
  sim_log_lognormal(
    n1 = 30,
    n2 = 30,
    ratio = 1.3,
    cv1 = 0.7,
    cv2 = 0.4,
    cor = 0.765,
    nsims = 10000
  ) |>
    power() |>
    getElement("power") |>
    as.numeric(),
  0.93598,
  tolerance = 0.02,
  scale = 1
)

# one-sample
expect_equal(
  sim_log_lognormal(
    n1 = 30,
    ratio = 1.3,
    cv1 = 0.35,
    nsims = 10000
  ) |>
    power() |>
    getElement("power") |>
    as.numeric(),
  0.983,
  tolerance = 0.02,
  scale = 1
)

#-------------------------------------------------------------------------------
# Validated NB power
# table 1 and 2 in Rettiganti (2012) doi:10.1080/10543406.2010.528105
#-------------------------------------------------------------------------------
expect_equal(
  sim_nb(
    n1 = 76,
    mean1 = 5.9,
    ratio = 0.5,
    dispersion1 = 0.49,
    nsims = 1000
  ) |>
    power(lrt_nb()) |>
    getElement("power") |>
    as.numeric(),
  0.8,
  tolerance = 0.04,
  scale = 1
)

expect_equal(
  sim_nb(
    n1 = 87,
    mean1 = 5.9,
    ratio = 0.5,
    dispersion1 = 0.49,
    dispersion2 = 0.3675,
    nsims = 1000
  ) |>
    power(lrt_nb()) |>
    getElement("power") |>
    as.numeric(),
  0.8,
  tolerance = 0.05,
  scale = 1
)

#-------------------------------------------------------------------------------
# Validated BNB power
# table 4 in Rettiganti (2012) doi:10.1080/10543406.2010.528105
#-------------------------------------------------------------------------------
expect_equal(
  sim_bnb(
    n = 10,
    mean1 = 5.9,
    ratio = 0.5,
    dispersion = 0.49,
    nsims = 1000
  ) |>
    power(lrt_bnb()) |>
    getElement("power") |>
    as.numeric(),
  0.82,
  tolerance = 0.04,
  scale = 1
)

#-------------------------------------------------------------------------------
# Data with majority zeros
#-------------------------------------------------------------------------------
nsims <- 200
set.seed(1234)
d <- sim_nb(
  n1 = c(10),
  n2 = c(10),
  mean1 = c(10),
  mean2 = c(15),
  dispersion1 = c(0.01),
  dispersion2 = c(0.01),
  nsims = nsims
) |>
  power()

expect_true(d$nsims < nsims*0.5)

nsims <- 200
set.seed(1234)
d <- sim_bnb(
  n = 10,
  mean1 = 10,
  ratio = 1.3,
  dispersion = 0.01,
  nsims = 200
) |>
  power() |>
  suppressWarnings()

expect_true(d$nsims < nsims*0.5)
