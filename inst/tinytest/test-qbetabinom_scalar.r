library(depower)
library(tinytest)

# Access internal function
qbetabinom_scalar <- depower:::qbetabinom_scalar

#-------------------------------------------------------------------------------
# Helper: compute Beta-Binomial CDF manually for verification
#-------------------------------------------------------------------------------
pbetabinom <- function(k, size, shape1, shape2) {
  if (k < 0) {
    return(0)
  }
  if (k >= size) {
    return(1)
  }
  x <- 0:k
  pmf <- choose(size, x) *
    exp(lbeta(x + shape1, size - x + shape2) - lbeta(shape1, shape2))
  sum(pmf)
}

#-------------------------------------------------------------------------------
# Edge cases: p <= 0 and p >= 1
#-------------------------------------------------------------------------------
expect_equal(
  qbetabinom_scalar(p = 0, size = 100, shape1 = 2, shape2 = 3),
  0L,
  info = "p = 0 returns 0"
)

expect_equal(
  qbetabinom_scalar(p = -0.5, size = 100, shape1 = 2, shape2 = 3),
  0L,
  info = "p < 0 returns 0"
)

expect_equal(
  qbetabinom_scalar(p = 1, size = 100, shape1 = 2, shape2 = 3),
  100L,
  info = "p = 1 returns size"
)

expect_equal(
  qbetabinom_scalar(p = 1.5, size = 100, shape1 = 2, shape2 = 3),
  100L,
  info = "p > 1 returns size"
)

#-------------------------------------------------------------------------------
# Edge case: size = 0
#-------------------------------------------------------------------------------
expect_equal(
  qbetabinom_scalar(p = 0.5, size = 0, shape1 = 2, shape2 = 3),
  0L,
  info = "size = 0 returns 0 for any p in (0,1)"
)

#-------------------------------------------------------------------------------
# Return type
#-------------------------------------------------------------------------------
expect_true(
  is.integer(qbetabinom_scalar(p = 0.5, size = 100, shape1 = 2, shape2 = 3)),
  info = "Return type is integer"
)

#-------------------------------------------------------------------------------
# Small size: verify against manual CDF computation
#-------------------------------------------------------------------------------
# For size = 5, alpha = 2, beta = 2, compute all CDF values
size <- 5L
alpha <- 2
beta <- 2
cdf_vals <- vapply(
  0:size,
  pbetabinom,
  numeric(1),
  size = size,
  shape1 = alpha,
  shape2 = beta
)

# Quantile is smallest k such that CDF(k) >= p
# Test p values that are not exactly at CDF boundaries to avoid floating-point issues
# CDF values are: 0.107, 0.286, 0.500, 0.714, 0.893, 1.000
# Use p values slightly below actual CDF values
for (p in c(0.1, 0.28, 0.49, 0.71, 0.89)) {
  expected <- min(which(cdf_vals >= p)) - 1L # -1 because cdf_vals is 0-indexed via 0:size
  result <- qbetabinom_scalar(p = p, size = size, shape1 = alpha, shape2 = beta)
  expect_equal(
    result,
    expected,
    info = sprintf("size=5, alpha=2, beta=2, p=%.2f", p)
  )
}

#-------------------------------------------------------------------------------
# Symmetric case: alpha = beta
#-------------------------------------------------------------------------------
# For symmetric Beta-Binomial, median should be near size/2
size <- 100L
result <- qbetabinom_scalar(p = 0.5, size = size, shape1 = 5, shape2 = 5)
expect_true(
  result >= 45L && result <= 55L,
  info = "Symmetric case: median near size/2"
)

#-------------------------------------------------------------------------------
# Skewed cases
#-------------------------------------------------------------------------------
# Large alpha relative to beta: distribution skewed toward size
result_high <- qbetabinom_scalar(p = 0.5, size = 100, shape1 = 10, shape2 = 1)
expect_true(
  result_high > 50L,
  info = "Large alpha/beta ratio: median > size/2"
)

# Small alpha relative to beta: distribution skewed toward 0
result_low <- qbetabinom_scalar(p = 0.5, size = 100, shape1 = 1, shape2 = 10)
expect_true(
  result_low < 50L,
  info = "Small alpha/beta ratio: median < size/2"
)

#-------------------------------------------------------------------------------
# Monotonicity: larger p gives larger or equal quantile
#-------------------------------------------------------------------------------
size <- 50L
alpha <- 3
beta <- 4
p_vals <- c(0.1, 0.25, 0.5, 0.75, 0.9)
q_vals <- vapply(
  p_vals,
  qbetabinom_scalar,
  integer(1),
  size = size,
  shape1 = alpha,
  shape2 = beta
)

expect_true(
  all(diff(q_vals) >= 0),
  info = "Quantiles are monotonically non-decreasing in p"
)

#-------------------------------------------------------------------------------
# Consistency: CDF at quantile should be >= p
#-------------------------------------------------------------------------------
size <- 30L
alpha <- 2.5
beta <- 3.5

for (p in c(0.05, 0.25, 0.5, 0.75, 0.95)) {
  q <- qbetabinom_scalar(p = p, size = size, shape1 = alpha, shape2 = beta)
  cdf_at_q <- pbetabinom(q, size = size, shape1 = alpha, shape2 = beta)
  expect_true(
    cdf_at_q >= p,
    info = sprintf("CDF at quantile >= p for p=%.2f", p)
  )
  # Also check that q-1 has CDF < p (unless q = 0)
  if (q > 0L) {
    cdf_at_q_minus_1 <- pbetabinom(
      q - 1L,
      size = size,
      shape1 = alpha,
      shape2 = beta
    )
    expect_true(
      cdf_at_q_minus_1 < p,
      info = sprintf("CDF at quantile-1 < p for p=%.2f", p)
    )
  }
}

#-------------------------------------------------------------------------------
# Large size: tests the underflow branch
#-------------------------------------------------------------------------------
# When size is large and alpha is small, P(Y=0) can underflow
size <- 10000L
alpha <- 0.5
beta <- 0.5

# Should not error and should return valid integer
result <- qbetabinom_scalar(p = 0.5, size = size, shape1 = alpha, shape2 = beta)
expect_true(
  is.integer(result) && !is.na(result),
  info = "Large size with small shape parameters: returns valid integer"
)
expect_true(
  result >= 0L && result <= size,
  info = "Large size: quantile within valid range"
)

# Verify monotonicity still holds for large size
q_low <- qbetabinom_scalar(p = 0.25, size = size, shape1 = alpha, shape2 = beta)
q_high <- qbetabinom_scalar(
  p = 0.75,
  size = size,
  shape1 = alpha,
  shape2 = beta
)
expect_true(
  q_low <= q_high,
  info = "Large size: monotonicity preserved"
)

#-------------------------------------------------------------------------------
# Large size with larger shape parameters (tests sequential branch)
#-------------------------------------------------------------------------------
size <- 500L
alpha <- 10
beta <- 10

result <- qbetabinom_scalar(p = 0.5, size = size, shape1 = alpha, shape2 = beta)
expect_true(
  is.integer(result) && !is.na(result),
  info = "Moderate size with larger shapes: returns valid integer"
)
# Symmetric case should have median near size/2
expect_true(
  result >= 200L && result <= 300L,
  info = "Moderate size symmetric case: median near size/2"
)

#-------------------------------------------------------------------------------
# Extreme probabilities near boundaries
#-------------------------------------------------------------------------------
size <- 100L
alpha <- 2
beta <- 2

# Very small p (but > 0)
result_small_p <- qbetabinom_scalar(
  p = 1e-10,
  size = size,
  shape1 = alpha,
  shape2 = beta
)
expect_true(
  result_small_p >= 0L,
  info = "Very small p returns non-negative quantile"
)

# p very close to 1 (but < 1)
result_large_p <- qbetabinom_scalar(
  p = 1 - 1e-10,
  size = size,
  shape1 = alpha,
  shape2 = beta
)
expect_true(
  result_large_p <= size,
  info = "p near 1 returns quantile <= size"
)

#-------------------------------------------------------------------------------
# Small shape parameters (Jeffreys prior)
#-------------------------------------------------------------------------------
size <- 50L
alpha <- 0.5
beta <- 0.5

result <- qbetabinom_scalar(p = 0.5, size = size, shape1 = alpha, shape2 = beta)
expect_true(
  is.integer(result) && result >= 0L && result <= size,
  info = "Jeffreys prior parameters: valid result"
)

#-------------------------------------------------------------------------------
# Non-integer shape parameters
#-------------------------------------------------------------------------------
size <- 40L
alpha <- 1.7
beta <- 2.3

result <- qbetabinom_scalar(p = 0.5, size = size, shape1 = alpha, shape2 = beta)
expect_true(
  is.integer(result) && result >= 0L && result <= size,
  info = "Non-integer shape parameters: valid result"
)

# Verify against manual CDF
cdf_at_result <- pbetabinom(result, size = size, shape1 = alpha, shape2 = beta)
expect_true(
  cdf_at_result >= 0.5,
  info = "Non-integer shapes: CDF at quantile >= p"
)
