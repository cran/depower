library(tinytest)

#-------------------------------------------------------------------------------
# fun_args
#-------------------------------------------------------------------------------
current <- depower:::fun_args(alist(wald_test_nb()))
target <- list(
  list(
    distribution = NA,
    nsims = NA,
    ncores = NA,
    alternative = NA,
    equal_dispersion = NA,
    ratio_null = NA,
    mean_null = NA
  )
)
expect_equal(current, target)

current <- depower:::fun_args(alist(wald_test_nb(distribution = asymptotic())))
target <- list(
  list(
    distribution = "asymptotic",
    nsims = NA,
    ncores = NA,
    alternative = NA,
    equal_dispersion = NA,
    ratio_null = NA,
    mean_null = NA
  )
)
expect_equal(current, target)

current <- depower:::fun_args(alist(wald_test_nb(distribution = simulated())))
target <- list(
  list(
    distribution = "simulated",
    nsims = 1000L,
    ncores = 1L,
    alternative = NA,
    equal_dispersion = NA,
    ratio_null = NA,
    mean_null = NA
  )
)
expect_equal(current, target)

current <- depower:::fun_args(
  alist(
    wald_test_nb(
      equal_dispersion = FALSE,
      ratio_null = 1.2,
      distribution = simulated(method = "approximate", nsims = 250, ncores = 4)
    )
  )
)
target <- list(
  list(
    distribution = "simulated",
    nsims = 250,
    ncores = 4,
    alternative = NA,
    equal_dispersion = FALSE,
    ratio_null = 1.2,
    mean_null = NA
  )
)
expect_equal(current, target)

current <- depower:::fun_args(
  alist(
    t_test_welch(
      alternative = "greater",
      mean_null = log(1.2)
    )
  )
)
target <- list(
  list(
    distribution = NA,
    nsims = NA,
    ncores = NA,
    alternative = "greater",
    equal_dispersion = NA,
    ratio_null = NA,
    mean_null = quote(log(1.2))
  )
)
expect_equal(current, target)

#-------------------------------------------------------------------------------
# data_and_fun_args()
#-------------------------------------------------------------------------------
data <- sim_nb(
  n1 = 10,
  n2 = 15,
  mean1 = 10,
  mean2 = 15,
  dispersion1 = 10,
  dispersion2 = 11,
  nsims = 100
)

.funs <- alist(
  "Asymptotic NB Wald Test" = wald_test_nb(),
  "Simulated NB Wald Test" = wald_test_nb(
    equal_dispersion = FALSE,
    link = "sqrt",
    ratio_null = 1.1,
    distribution = simulated(nsims = 150)
  )
)

res <- depower:::data_and_fun_args(data = data, .funs = .funs)
expect_equal(
  res,
  structure(
    list(
      n1 = structure(c(10, 10), label = "n1"),
      n2 = structure(c(15, 15), label = "n2"),
      mean1 = structure(c(10, 10), label = "Mean1"),
      mean2 = structure(c(15, 15), label = "Mean2"),
      ratio = structure(c(1.5, 1.5), label = "Ratio"),
      dispersion1 = structure(c(10, 10), label = "Dispersion1"),
      dispersion2 = structure(c(11, 11), label = "Dispersion2"),
      nsims = structure(c(100, 100), label = "N Simulations"),
      distribution = structure(
        c("Independent two-sample NB", "Independent two-sample NB"),
        label = "Distribution"
      ),
      test = c(
        `wald_test_nb()` = "Asymptotic NB Wald Test",
        `wald_test_nb(equal_dispersion = FALSE, link = "sqrt", ratio_null = 1.1, distribution = simulated(nsims = 150))` = "Simulated NB Wald Test"
      ),
      distribution_test_stat = c(NA, "simulated"),
      nsims_test_stat = c(NA, 150),
      ncores_test_stat = c(NA, 1L),
      alternative = c(NA, NA),
      equal_dispersion = c(NA, FALSE),
      ratio_null = c(NA, 1.1),
      mean_null = c(NA, NA)
    ),
    row.names = c(NA, -2L),
    class = c("depower", "nb", "tbl_df", "tbl", "data.frame")
  )
)

#-------------------------------------------------------------------------------
# sim_null()
#-------------------------------------------------------------------------------
set.seed(1234)
data <- sim_nb(
  n1 = 10,
  n2 = 15,
  mean1 = 10,
  mean2 = 30,
  dispersion1 = 10,
  dispersion2 = 11,
  nsims = 2
)

.funs <- alist(
  "Asymptotic NB Wald Test" = wald_test_nb(),
  "Simulated NB Wald Test" = wald_test_nb(
    equal_dispersion = FALSE,
    link = "sqrt",
    ratio_null = 1.1,
    distribution = simulated(nsims = 2)
  )
)

expect_equal(
  depower:::sim_null(data = data, .funs = .funs),
  structure(
    list(
      n1 = structure(c(10, 10), label = "n1"),
      n2 = structure(c(15, 15), label = "n2"),
      mean1 = structure(c(10, 10), label = "Mean1"),
      mean2 = structure(c(30, 30), label = "Mean2"),
      ratio = structure(c(3, 3), label = "Ratio"),
      dispersion1 = structure(c(10, 10), label = "Dispersion1"),
      dispersion2 = structure(c(11, 11), label = "Dispersion2"),
      nsims = structure(c(2, 2), label = "N Simulations"),
      distribution = structure(
        c("Independent two-sample NB", "Independent two-sample NB"),
        label = "Distribution"
      ),
      test = c(
        `wald_test_nb()` = "Asymptotic NB Wald Test",
        `wald_test_nb(equal_dispersion = FALSE, link = "sqrt", ratio_null = 1.1, distribution = simulated(nsims = 2))` = "Simulated NB Wald Test"
      ),
      distribution_test_stat = c(NA, "simulated"),
      nsims_test_stat = c(NA, 2),
      ncores_test_stat = c(NA, 1L),
      alternative = c(NA, NA),
      equal_dispersion = c(NA, FALSE),
      ratio_null = c(NA, 1.1),
      mean_null = c(NA, NA),
      data_null = list(
        NULL,
        list(
          list(
            value1 = c(15, 13, 6, 14, 2, 7, 8, 10, 11, 9),
            value2 = c(12, 10, 6, 17, 10, 14, 7, 13, 8, 14, 18, 4, 16, 11, 12)
          ),
          list(
            value1 = c(16, 13, 8, 14, 19, 11, 5, 5, 12, 17),
            value2 = c(16, 9, 10, 9, 16, 8, 12, 10, 4, 11, 13, 15, 20, 15, 14)
          )
        )
      )
    ),
    class = c("tbl_df", "tbl", "data.frame"),
    row.names = c(NA, -2L)
  )
)
