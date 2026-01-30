library(depower)
library(tinytest)

# Access internal function
qbetabinom <- depower:::qbetabinom

#-------------------------------------------------------------------------------
# Scalar input: dispatches to qbetabinom_scalar
#-------------------------------------------------------------------------------
expect_equal(
  qbetabinom(p = 0.5, size = 100, shape1 = 2, shape2 = 3),
  depower:::qbetabinom_scalar(p = 0.5, size = 100, shape1 = 2, shape2 = 3),
  info = "Scalar input matches qbetabinom_scalar"
)

#-------------------------------------------------------------------------------
# Return type
#-------------------------------------------------------------------------------
expect_true(
  is.integer(qbetabinom(p = 0.5, size = 100, shape1 = 2, shape2 = 3)),
  info = "Scalar input returns integer"
)

expect_true(
  is.integer(qbetabinom(p = c(0.25, 0.75), size = 100, shape1 = 2, shape2 = 3)),
  info = "Vector input returns integer vector"
)

#-------------------------------------------------------------------------------
# Vectorized over p
#-------------------------------------------------------------------------------
p_vals <- c(0.1, 0.25, 0.5, 0.75, 0.9)
result <- qbetabinom(p = p_vals, size = 100, shape1 = 5, shape2 = 5)

expect_equal(
  length(result),
  length(p_vals),
  info = "Vectorized over p: correct length"
)

# Verify each element matches scalar call
for (i in seq_along(p_vals)) {
  expect_equal(
    result[i],
    depower:::qbetabinom_scalar(p_vals[i], 100, 5, 5),
    info = sprintf("Vectorized over p: element %d matches scalar", i)
  )
}

# Monotonicity
expect_true(
  all(diff(result) >= 0),
  info = "Vectorized over p: monotonically non-decreasing"
)

#-------------------------------------------------------------------------------
# Vectorized over size
#-------------------------------------------------------------------------------
sizes <- c(50L, 100L, 200L)
result <- qbetabinom(p = 0.5, size = sizes, shape1 = 5, shape2 = 5)

expect_equal(
  length(result),
  length(sizes),
  info = "Vectorized over size: correct length"
)

for (i in seq_along(sizes)) {
  expect_equal(
    result[i],
    depower:::qbetabinom_scalar(0.5, sizes[i], 5, 5),
    info = sprintf("Vectorized over size: element %d matches scalar", i)
  )
}

#-------------------------------------------------------------------------------
# Vectorized over shape1
#-------------------------------------------------------------------------------
shape1_vals <- c(1, 5, 10)
result <- qbetabinom(p = 0.5, size = 100, shape1 = shape1_vals, shape2 = 5)

expect_equal(
  length(result),
  length(shape1_vals),
  info = "Vectorized over shape1: correct length"
)

for (i in seq_along(shape1_vals)) {
  expect_equal(
    result[i],
    depower:::qbetabinom_scalar(0.5, 100, shape1_vals[i], 5),
    info = sprintf("Vectorized over shape1: element %d matches scalar", i)
  )
}

# Larger shape1 should give larger quantile (distribution shifts right)
expect_true(
  all(diff(result) >= 0),
  info = "Vectorized over shape1: larger shape1 gives larger quantile"
)

#-------------------------------------------------------------------------------
# Vectorized over shape2
#-------------------------------------------------------------------------------
shape2_vals <- c(1, 5, 10)
result <- qbetabinom(p = 0.5, size = 100, shape1 = 5, shape2 = shape2_vals)

expect_equal(
  length(result),
  length(shape2_vals),
  info = "Vectorized over shape2: correct length"
)

for (i in seq_along(shape2_vals)) {
  expect_equal(
    result[i],
    depower:::qbetabinom_scalar(0.5, 100, 5, shape2_vals[i]),
    info = sprintf("Vectorized over shape2: element %d matches scalar", i)
  )
}

# Larger shape2 should give smaller quantile (distribution shifts left)
expect_true(
  all(diff(result) <= 0),
  info = "Vectorized over shape2: larger shape2 gives smaller quantile"
)

#-------------------------------------------------------------------------------
# Vectorized over all arguments (same length)
#-------------------------------------------------------------------------------
p_vec <- c(0.25, 0.5, 0.75)
size_vec <- c(50L, 100L, 150L)
shape1_vec <- c(2, 5, 10)
shape2_vec <- c(3, 5, 2)

result <- qbetabinom(
  p = p_vec,
  size = size_vec,
  shape1 = shape1_vec,
  shape2 = shape2_vec
)

expect_equal(
  length(result),
  3L,
  info = "Vectorized over all: correct length"
)

for (i in 1:3) {
  expect_equal(
    result[i],
    depower:::qbetabinom_scalar(
      p_vec[i],
      size_vec[i],
      shape1_vec[i],
      shape2_vec[i]
    ),
    info = sprintf("Vectorized over all: element %d matches scalar", i)
  )
}

#-------------------------------------------------------------------------------
# Recycling: length-1 recycled to match longer vector
#-------------------------------------------------------------------------------
p_vec <- c(0.25, 0.5, 0.75)
result <- qbetabinom(p = p_vec, size = 100, shape1 = 5, shape2 = 5)

expect_equal(
  length(result),
  3L,
  info = "Recycling length-1: correct output length"
)

# Equivalent to explicit recycling
result_explicit <- qbetabinom(
  p = p_vec,
  size = rep(100L, 3),
  shape1 = rep(5, 3),
  shape2 = rep(5, 3)
)
expect_equal(
  result,
  result_explicit,
  info = "Recycling length-1: matches explicit recycling"
)

#-------------------------------------------------------------------------------
# Recycling: length-2 recycled to length-4
#-------------------------------------------------------------------------------
p_vec <- c(0.25, 0.5, 0.75, 0.9)
shape1_vec <- c(2, 8)

result <- qbetabinom(p = p_vec, size = 100, shape1 = shape1_vec, shape2 = 5)

expect_equal(
  length(result),
  4L,
  info = "Recycling length-2 to length-4: correct output length"
)

# shape1 recycled as c(2, 8, 2, 8)
expected <- vapply(
  1:4,
  function(i) {
    depower:::qbetabinom_scalar(
      p_vec[i],
      100,
      shape1_vec[((i - 1) %% 2) + 1],
      5
    )
  },
  integer(1)
)
expect_equal(
  result,
  expected,
  info = "Recycling length-2 to length-4: correct values"
)

#-------------------------------------------------------------------------------
# Error on incompatible lengths
#-------------------------------------------------------------------------------
expect_error(
  qbetabinom(p = c(0.25, 0.5, 0.75), size = c(50, 100), shape1 = 5, shape2 = 5),
  pattern = "Lengths.*not compatible",
  info = "Error on incompatible lengths: 3 and 2"
)

expect_error(
  qbetabinom(p = 0.5, size = c(50, 100, 150), shape1 = c(2, 4), shape2 = 5),
  pattern = "Lengths.*not compatible",
  info = "Error on incompatible lengths: 3 and 2"
)

expect_error(
  qbetabinom(
    p = c(0.1, 0.2, 0.3, 0.4, 0.5),
    size = 100,
    shape1 = c(1, 2, 3),
    shape2 = 5
  ),
  pattern = "Lengths.*not compatible",
  info = "Error on incompatible lengths: 5 and 3"
)

#-------------------------------------------------------------------------------
# Edge cases with vectorization
#-------------------------------------------------------------------------------
# Mix of edge case p values
p_edge <- c(0, 0.5, 1)
result <- qbetabinom(p = p_edge, size = 100, shape1 = 5, shape2 = 5)

expect_equal(
  result[1],
  0L,
  info = "Vectorized edge case: p=0 returns 0"
)
expect_equal(
  result[3],
  100L,
  info = "Vectorized edge case: p=1 returns size"
)

#-------------------------------------------------------------------------------
# Large vectors
#-------------------------------------------------------------------------------
n_large <- 100L
p_large <- seq(0.01, 0.99, length.out = n_large)

result <- qbetabinom(p = p_large, size = 200, shape1 = 10, shape2 = 10)

expect_equal(
  length(result),
  n_large,
  info = "Large vector: correct length"
)

expect_true(
  all(is.integer(result)),
  info = "Large vector: all integers"
)

expect_true(
  all(result >= 0L & result <= 200L),
  info = "Large vector: all within valid range"
)

expect_true(
  all(diff(result) >= 0),
  info = "Large vector: monotonically non-decreasing"
)

#-------------------------------------------------------------------------------
# Consistency with examples from documentation
#-------------------------------------------------------------------------------
# Single quantile example
result_single <- qbetabinom(0.5, size = 100, shape1 = 2, shape2 = 2)
expect_true(
  is.integer(result_single) && length(result_single) == 1L,
  info = "Documentation example: single quantile"
)

# Multiple quantiles example
result_multi <- qbetabinom(c(0.025, 0.975), size = 100, shape1 = 10, shape2 = 5)
expect_equal(
  length(result_multi),
  2L,
  info = "Documentation example: multiple quantiles length"
)
expect_true(
  result_multi[1] < result_multi[2],
  info = "Documentation example: lower quantile < upper quantile"
)

# Fully vectorized example
result_full <- qbetabinom(
  p = c(0.025, 0.975),
  size = c(100, 200),
  shape1 = c(10, 20),
  shape2 = c(5, 10)
)
expect_equal(
  length(result_full),
  2L,
  info = "Documentation example: fully vectorized length"
)
