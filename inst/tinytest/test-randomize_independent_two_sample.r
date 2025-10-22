library(tinytest)

#-------------------------------------------------------------------------------
# Data
#-------------------------------------------------------------------------------
df1 <- list(value1 = c(1:3), value2 = c(11:13))
df2 <- list(value1 = c(1:3, NA), value2 = c(11:13))

#-------------------------------------------------------------------------------
# Tests
#-------------------------------------------------------------------------------
set.seed(123)
res <- depower:::randomize_independent_two_sample(
  data = df1,
  distribution = simulated(nsims = 5)
)
expect_equal(
  res,
  list(
    list(value1 = 1:3, value2 = 11:13),
    list(value1 = c(3L, 13L, 2L), value2 = c(11L, 12L, 1L)),
    list(value1 = c(12L, 11L, 2L), value2 = c(13L, 1L, 3L)),
    list(value1 = c(3L, 12L, 13L), value2 = c(11L, 1L, 2L)),
    list(value1 = c(1L, 13L, 12L), value2 = c(3L, 2L, 11L)),
    list(value1 = c(2L, 1L, 13L), value2 = c(3L, 11L, 12L))
  )
)

set.seed(123)
res <- depower:::randomize_independent_two_sample(
  data = df2,
  distribution = simulated(nsims = 5)
)
expect_equal(
  res,
  list(
    list(value1 = c(1L, 2L, 3L, NA), value2 = 11:13),
    list(value1 = c(13L, 3L, 12L, 2L), value2 = c(NA, 11L, 1L)),
    list(value1 = c(11L, NA, 1L, 2L), value2 = c(3L, 12L, 13L)),
    list(value1 = c(3L, 13L, 1L, NA), value2 = c(11L, 12L, 2L)),
    list(value1 = c(3L, 2L, 12L, 1L), value2 = c(11L, NA, 13L)),
    list(value1 = c(12L, 1L, 3L, 11L), value2 = c(13L, 2L, NA))
  )
)
