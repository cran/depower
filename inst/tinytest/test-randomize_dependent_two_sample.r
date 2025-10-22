library(tinytest)

#-------------------------------------------------------------------------------
# Data
#-------------------------------------------------------------------------------
df1 <- list(value1 = c(1:3), value2 = c(11:13))
df2 <- list(value1 = c(1:3, NA), value2 = c(11:14))

#-------------------------------------------------------------------------------
# Tests
#-------------------------------------------------------------------------------
set.seed(123)
res <- depower:::randomize_dependent_two_sample(
  data = df1,
  distribution = simulated(nsims = 5)
)
expect_equal(
  res,
  list(
    list(value1 = 1:3, value2 = 11:13),
    list(value1 = 1:3, value2 = 11:13),
    list(value1 = c(2L, 11L, 13L), value2 = c(12L, 1L, 3L)),
    list(value1 = c(3L, 11L, 12L), value2 = c(13L, 1L, 2L)),
    list(value1 = c(1L, 12L, 13L), value2 = c(11L, 2L, 3L)),
    list(value1 = c(2L, 11L, 13L), value2 = c(12L, 1L, 3L))
  )
)

set.seed(123)
res2 <- depower:::randomize_dependent_two_sample(
  data = df2,
  distribution = simulated(nsims = 5)
)
expect_equal(
  res2,
  list(
    list(value1 = c(1L, 2L, 3L, NA), value2 = 11:14),
    list(value1 = 1:3, value2 = 11:13),
    list(value1 = c(2L, 11L, 13L), value2 = c(12L, 1L, 3L)),
    list(value1 = c(3L, 11L, 12L), value2 = c(13L, 1L, 2L)),
    list(value1 = c(1L, 12L, 13L), value2 = c(11L, 2L, 3L)),
    list(value1 = c(2L, 11L, 13L), value2 = c(12L, 1L, 3L))
  )
)

expect_error(
  depower:::randomize_dependent_two_sample(
    data = list(1:3, 1:4),
    distribution = simulated(nsims = 5)
  ),
  pattern = "Argument 'data' must have the same sample size for both samples."
)
