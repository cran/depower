library(tinytest)

#-------------------------------------------------------------------------------
# Data
#-------------------------------------------------------------------------------
df1 <- list(value1 = c(1:3), value2 = c(11:13))
df2 <- list(value1 = c(1:2), value2 = c(11:13))

#-------------------------------------------------------------------------------
# Tests
#-------------------------------------------------------------------------------
expect_equal(
  depower:::permute_independent_two_sample(
    data = df1,
    distribution = simulated()
  ),
  list(
    list(value1 = 1:3, value2 = 11:13),
    list(value1 = c(1L, 2L, 11L), value2 = c(3L, 12L, 13L)),
    list(value1 = c(1L, 2L, 12L), value2 = c(3L, 11L, 13L)),
    list(value1 = c(1L, 2L, 13L), value2 = c(3L, 11L, 12L)),
    list(value1 = c(1L, 3L, 11L), value2 = c(2L, 12L, 13L)),
    list(value1 = c(1L, 3L, 12L), value2 = c(2L, 11L, 13L)),
    list(value1 = c(1L, 3L, 13L), value2 = c(2L, 11L, 12L)),
    list(value1 = c(1L, 11L, 12L), value2 = c(2L, 3L, 13L)),
    list(value1 = c(1L, 11L, 13L), value2 = c(2L, 3L, 12L)),
    list(value1 = c(1L, 12L, 13L), value2 = c(2L, 3L, 11L)),
    list(value1 = c(2L, 3L, 11L), value2 = c(1L, 12L, 13L)),
    list(value1 = c(2L, 3L, 12L), value2 = c(1L, 11L, 13L)),
    list(value1 = c(2L, 3L, 13L), value2 = c(1L, 11L, 12L)),
    list(value1 = c(2L, 11L, 12L), value2 = c(1L, 3L, 13L)),
    list(value1 = c(2L, 11L, 13L), value2 = c(1L, 3L, 12L)),
    list(value1 = c(2L, 12L, 13L), value2 = c(1L, 3L, 11L)),
    list(value1 = c(3L, 11L, 12L), value2 = c(1L, 2L, 13L)),
    list(value1 = c(3L, 11L, 13L), value2 = c(1L, 2L, 12L)),
    list(value1 = c(3L, 12L, 13L), value2 = c(1L, 2L, 11L)),
    list(value1 = 11:13, value2 = 1:3)
  )
)

expect_equal(
  depower:::permute_independent_two_sample(
    data = df2,
    distribution = simulated()
  ),
  list(
    list(value1 = 1:2, value2 = 11:13),
    list(value1 = c(1L, 11L), value2 = c(2L, 12L, 13L)),
    list(value1 = c(1L, 12L), value2 = c(2L, 11L, 13L)),
    list(value1 = c(1L, 13L), value2 = c(2L, 11L, 12L)),
    list(value1 = c(2L, 11L), value2 = c(1L, 12L, 13L)),
    list(value1 = c(2L, 12L), value2 = c(1L, 11L, 13L)),
    list(value1 = c(2L, 13L), value2 = c(1L, 11L, 12L)),
    list(value1 = 11:12, value2 = c(1L, 2L, 13L)),
    list(value1 = c(11L, 13L), value2 = c(1L, 2L, 12L)),
    list(value1 = 12:13, value2 = c(1L, 2L, 11L))
  )
)

expect_warning(
  depower:::permute_independent_two_sample(
    data = list(1:20, 1:20),
    distribution = simulated()
  ),
  pattern = "Using approximate randomization with 1000 resamples. The exact randomization test is only used when choose(20+20,20)=1.38e+11 < 1e6",
  fixed = TRUE
)
