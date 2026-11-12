# Near-equal test statistics count as ties in randomization and simulated
# p-values, for pchisq2() and for the exact randomization tests that use it.
library(depower)
library(tinytest)

#-------------------------------------------------------------------------------
# pchisq2()
#-------------------------------------------------------------------------------
# Upper tail. The values 9.995 and 12 are counted, and 9.98 and 3 are not.
expect_equal(
  depower:::pchisq2(q = 10, q_null = c(9.995, 9.98, 12, 3)),
  3 / 5
)

# Lower tail. The values 10.005 and 3 are counted, and 10.02 and 12 are not.
expect_equal(
  depower:::pchisq2(
    q = 10,
    lower.tail = TRUE,
    q_null = c(10.005, 10.02, 12, 3)
  ),
  3 / 5
)

# A statistic of zero has no tolerance, and every value is counted.
expect_equal(depower:::pchisq2(q = 0, q_null = c(0, 1)), 1)

# The asymptotic p-value does not use the tolerance.
expect_identical(
  depower:::pchisq2(q = 3.84),
  pchisq(3.84, df = 1L, lower.tail = FALSE)
)

#-------------------------------------------------------------------------------
# Exact randomization tests
#-------------------------------------------------------------------------------
# With the log link, the BNB Wald statistic and the LRT depend only on the two
# totals. A swap of all pairs then gives a tied statistic.
# The expected p-value uses exact totals for all 16 swap patterns.
#
# The optimizer error for these datasets is much smaller than the tolerance.
# The comparison of p-values is exact because a p-value is a count divided by 17.
exact_p <- function(x, y, statistic) {
  patterns <- as.matrix(expand.grid(rep(list(0:1), length(x))))
  totals <- cbind(
    as.vector((1L - patterns) %*% x + patterns %*% y),
    as.vector(patterns %*% x + (1L - patterns) %*% y)
  )
  values <- statistic(totals[, 1L], totals[, 2L])
  (sum(values >= statistic(sum(x), sum(y)) * (1 - 1e-8)) + 1L) /
    (nrow(patterns) + 1L)
}
wald_totals <- function(x, y) {
  log(y / x)^2 / (1 / x + 1 / y)
}
lrt_totals <- function(x, y) {
  s <- x + y
  2 * (x * log(2 * x / s) + y * log(2 * y / s))
}

data_ties <- list(
  concordant = list(c(3L, 8L, 14L, 30L), c(9L, 20L, 41L, 80L)),
  one_reversed = list(c(2L, 4L, 25L, 8L), c(5L, 9L, 22L, 39L)),
  moderate = list(c(12L, 7L, 20L, 15L), c(30L, 25L, 41L, 44L))
)
for (name in names(data_ties)) {
  d <- data_ties[[name]]
  expect_equal(
    wald_test_bnb(d, distribution = simulated(method = "exact"))$p,
    exact_p(d[[1L]], d[[2L]], wald_totals),
    info = paste("wald_test_bnb", name)
  )
  expect_equal(
    lrt_bnb(d, distribution = simulated(method = "exact"))$p,
    exact_p(d[[1L]], d[[2L]], lrt_totals),
    info = paste("lrt_bnb", name)
  )
}

# Concordant pairs have one tied pair of largest statistics, the smallest
# possible count for a statistic that a swap of all pairs does not change.
expect_equal(
  exact_p(
    data_ties$concordant[[1L]],
    data_ties$concordant[[2L]],
    wald_totals
  ),
  3 / 17
)
