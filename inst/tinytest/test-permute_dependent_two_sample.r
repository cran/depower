library(tinytest)

#-------------------------------------------------------------------------------
# Data
#-------------------------------------------------------------------------------
df1 <- list(value1 = c(1:3), value2 = c(11:13))
df2 <- list(value1 = c(1:3, NA), value2 = c(11:14))

#-------------------------------------------------------------------------------
# Tests
#-------------------------------------------------------------------------------
expect_equal(
  depower:::permute_dependent_two_sample(
    data = df1,
    distribution = simulated()
  ),
  list(
    list(value1 = 1:3, value2 = 11:13),
    list(value1 = c(11L, 2L, 3L), value2 = c(1L, 12L, 13L)),
    list(value1 = c(12L, 1L, 3L), value2 = c(2L, 11L, 13L)),
    list(value1 = c(11L, 12L, 3L), value2 = c(1L, 2L, 13L)),
    list(value1 = c(13L, 1L, 2L), value2 = c(3L, 11L, 12L)),
    list(value1 = c(11L, 13L, 2L), value2 = c(1L, 3L, 12L)),
    list(value1 = c(12L, 13L, 1L), value2 = c(2L, 3L, 11L)),
    list(value1 = 11:13, value2 = 1:3)
  )
)

expect_equal(
  depower:::permute_dependent_two_sample(
    data = df2,
    distribution = simulated()
  ),
  list(
    list(value1 = 1:3, value2 = 11:13),
    list(value1 = c(11L, 2L, 3L), value2 = c(1L, 12L, 13L)),
    list(value1 = c(12L, 1L, 3L), value2 = c(2L, 11L, 13L)),
    list(value1 = c(11L, 12L, 3L), value2 = c(1L, 2L, 13L)),
    list(value1 = c(13L, 1L, 2L), value2 = c(3L, 11L, 12L)),
    list(value1 = c(11L, 13L, 2L), value2 = c(1L, 3L, 12L)),
    list(value1 = c(12L, 13L, 1L), value2 = c(2L, 3L, 11L)),
    list(value1 = 11:13, value2 = 1:3)
  )
)

expect_error(
  depower:::permute_dependent_two_sample(
    data = list(1:3, 1:4),
    distribution = simulated()
  ),
  pattern = "Argument 'data' must have the same sample size for both samples."
)

expect_warning(
  depower:::permute_dependent_two_sample(
    data = list(1:21, 1:21),
    distribution = simulated()
  ),
  pattern = "Using approximate randomization with 1000 resamples. The exact randomization test is only used for sample size per group 20 or fewer (max 2^20=1048576 resamples).",
  fixed = TRUE
)
