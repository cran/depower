library(tinytest)

#-------------------------------------------------------------------------------
# Data
#-------------------------------------------------------------------------------
set.seed(20231201)
df <- MASS::mvrnorm(
  n = 30,
  mu = c(2, 3),
  Sigma = matrix(c(2, 0.5, 0.5, 2), 2L, 2L)
)

#-------------------------------------------------------------------------------
# Structure
#-------------------------------------------------------------------------------
# Dependent two sample
expect_equal(
  t_test_paired(
    data = list(value1 = df[, 1], value2 = df[, 2]),
    alternative = "two.sided",
    ci_level = 0.95
  ),
  list(
    t = 2.73630930779187,
    df = 29L,
    p = 0.0104953703353502,
    mean_diff = list(
      estimate = 1.03397837894942,
      lower = 0.261140591960272,
      upper = 1.80681616593858
    ),
    n = 30L,
    method = "One-sample t-test for 'two-sided' alternative",
    alternative = "two.sided",
    ci_level = 0.95,
    mean_null = 0
  ),
  tolerance = 0.001,
  scale = 1
)

# Dependent one sample
expect_equal(
  t_test_paired(
    data = list(value1 = df[, 2] - df[, 1]),
    alternative = "two.sided",
    ci_level = 0.95
  ),
  list(
    t = 2.73630930779187,
    df = 29L,
    p = 0.0104953703353502,
    mean_diff = list(
      estimate = 1.03397837894942,
      lower = 0.261140591960272,
      upper = 1.80681616593858
    ),
    n = 30L,
    method = "One-sample t-test for 'two-sided' alternative",
    alternative = "two.sided",
    ci_level = 0.95,
    mean_null = 0
  ),
  tolerance = 0.001,
  scale = 1
)

#-------------------------------------------------------------------------------
# Test expected p-value
#-------------------------------------------------------------------------------
# Dependent two sample
expect_equal(
  t_test_paired(
    data = list(value1 = df[, 1], value2 = df[, 2]),
    alternative = "two.sided"
  )$p,
  t.test(
    x = df[, 2],
    y = df[, 1],
    alternative = "two.sided",
    paired = TRUE
  )$p.value,
  tolerance = 0.0001,
  scale = 1
)

expect_equal(
  t_test_paired(
    data = list(value1 = df[, 1], value2 = df[, 2]),
    alternative = "greater"
  )$p,
  t.test(
    x = df[, 2],
    y = df[, 1],
    alternative = "greater",
    paired = TRUE
  )$p.value,
  tolerance = 0.0001,
  scale = 1
)

expect_equal(
  t_test_paired(
    data = list(value1 = df[, 1], value2 = df[, 2]),
    alternative = "less"
  )$p,
  t.test(
    x = df[, 2],
    y = df[, 1],
    alternative = "less",
    paired = TRUE
  )$p.value,
  tolerance = 0.0001,
  scale = 1
)

# Dependent one sample
expect_equal(
  t_test_paired(
    data = list(value1 = df[, 2] - df[, 1]),
    alternative = "two.sided"
  )$p,
  t.test(
    x = df[, 2] - df[, 1],
    alternative = "two.sided"
  )$p.value,
  tolerance = 0.0001,
  scale = 1
)

expect_equal(
  t_test_paired(
    data = list(value1 = df[, 2] - df[, 1]),
    alternative = "greater"
  )$p,
  t.test(
    x = df[, 2] - df[, 1],
    alternative = "greater"
  )$p.value,
  tolerance = 0.0001,
  scale = 1
)

expect_equal(
  t_test_paired(
    data = list(value1 = df[, 2] - df[, 1]),
    alternative = "less"
  )$p,
  t.test(
    x = df[, 2] - df[, 1],
    alternative = "less"
  )$p.value,
  tolerance = 0.0001,
  scale = 1
)

# Dependent two sample vs. Dependent one sample
expect_equal(
  t_test_paired(
    data = list(value1 = df[, 1], value2 = df[, 2]),
    alternative = "two.sided"
  ),
  t_test_paired(
    data = list(value1 = df[, 2] - df[, 1]),
    alternative = "two.sided"
  ),
  tolerance = 0.0001,
  scale = 1
)

expect_equal(
  t_test_paired(
    data = list(value1 = df[, 1], value2 = df[, 2]),
    alternative = "greater"
  ),
  t_test_paired(
    data = list(value1 = df[, 2] - df[, 1]),
    alternative = "greater"
  ),
  tolerance = 0.0001,
  scale = 1
)

expect_equal(
  t_test_paired(
    data = list(value1 = df[, 1], value2 = df[, 2]),
    alternative = "less"
  ),
  t_test_paired(
    data = list(value1 = df[, 2] - df[, 1]),
    alternative = "less"
  ),
  tolerance = 0.0001,
  scale = 1
)

#-------------------------------------------------------------------------------
# Test expected CI
#-------------------------------------------------------------------------------
# Dependent two sample
expect_equal(
  t_test_paired(
    data = list(value1 = df[, 1], value2 = df[, 2]),
    alternative = "two.sided",
    ci_level = 0.95
  )$mean_diff$lower,
  t.test(
    x = df[, 2],
    y = df[, 1],
    alternative = "two.sided",
    paired = TRUE
  )$conf.int[1],
  tolerance = 0.0001,
  scale = 1
)

expect_equal(
  t_test_paired(
    data = list(value1 = df[, 1], value2 = df[, 2]),
    alternative = "two.sided",
    ci_level = 0.95
  )$mean_diff$upper,
  t.test(
    x = df[, 2],
    y = df[, 1],
    alternative = "two.sided",
    paired = TRUE
  )$conf.int[2],
  tolerance = 0.0001,
  scale = 1
)

expect_equal(
  t_test_paired(
    data = list(value1 = df[, 1], value2 = df[, 2]),
    alternative = "greater",
    ci_level = 0.95
  )$mean_diff$lower,
  t.test(
    x = df[, 2],
    y = df[, 1],
    alternative = "greater",
    paired = TRUE
  )$conf.int[1],
  tolerance = 0.0001,
  scale = 1
)

expect_equal(
  t_test_paired(
    data = list(value1 = df[, 1], value2 = df[, 2]),
    alternative = "greater",
    ci_level = 0.95
  )$mean_diff$upper,
  t.test(
    x = df[, 2],
    y = df[, 1],
    alternative = "greater",
    paired = TRUE
  )$conf.int[2],
  tolerance = 0.0001,
  scale = 1
)

expect_equal(
  t_test_paired(
    data = list(value1 = df[, 1], value2 = df[, 2]),
    alternative = "less",
    ci_level = 0.95
  )$mean_diff$lower,
  t.test(
    x = df[, 2],
    y = df[, 1],
    alternative = "less",
    paired = TRUE
  )$conf.int[1],
  tolerance = 0.0001,
  scale = 1
)

expect_equal(
  t_test_paired(
    data = list(value1 = df[, 1], value2 = df[, 2]),
    alternative = "less",
    ci_level = 0.95
  )$mean_diff$upper,
  t.test(
    x = df[, 2],
    y = df[, 1],
    alternative = "less",
    paired = TRUE
  )$conf.int[2],
  tolerance = 0.0001,
  scale = 1
)

# Dependent one sample
expect_equal(
  t_test_paired(
    data = list(value1 = df[, 2] - df[, 1]),
    alternative = "two.sided",
    ci_level = 0.95
  )$mean_diff$lower,
  t.test(
    x = df[, 2] - df[, 1],
    alternative = "two.sided"
  )$conf.int[1],
  tolerance = 0.0001,
  scale = 1
)

expect_equal(
  t_test_paired(
    data = list(value1 = df[, 2] - df[, 1]),
    alternative = "two.sided",
    ci_level = 0.95
  )$mean_diff$upper,
  t.test(
    x = df[, 2] - df[, 1],
    alternative = "two.sided"
  )$conf.int[2],
  tolerance = 0.0001,
  scale = 1
)

expect_equal(
  t_test_paired(
    data = list(value1 = df[, 2] - df[, 1]),
    alternative = "greater",
    ci_level = 0.95
  )$mean_diff$lower,
  t.test(
    x = df[, 2] - df[, 1],
    alternative = "greater"
  )$conf.int[1],
  tolerance = 0.0001,
  scale = 1
)

expect_equal(
  t_test_paired(
    data = list(value1 = df[, 2] - df[, 1]),
    alternative = "greater",
    ci_level = 0.95
  )$mean_diff$upper,
  t.test(
    x = df[, 2] - df[, 1],
    alternative = "greater"
  )$conf.int[2],
  tolerance = 0.0001,
  scale = 1
)

expect_equal(
  t_test_paired(
    data = list(value1 = df[, 2] - df[, 1]),
    alternative = "less",
    ci_level = 0.95
  )$mean_diff$lower,
  t.test(
    x = df[, 2] - df[, 1],
    alternative = "less"
  )$conf.int[1],
  tolerance = 0.0001,
  scale = 1
)

expect_equal(
  t_test_paired(
    data = list(value1 = df[, 2] - df[, 1]),
    alternative = "less",
    ci_level = 0.95
  )$mean_diff$upper,
  t.test(
    x = df[, 2] - df[, 1],
    alternative = "less"
  )$conf.int[2],
  tolerance = 0.0001,
  scale = 1
)

#-------------------------------------------------------------------------------
# Test absence of CI
#-------------------------------------------------------------------------------
# Dependent two sample
expect_equal(
  t_test_paired(
    data = list(value1 = df[, 1], value2 = df[, 2]),
    alternative = "less"
  )$mean_diff[2:3],
  list(NA_real_, NA_real_),
  check.attributes = FALSE
)

# Dependent one sample
expect_equal(
  t_test_paired(
    data = list(value1 = df[, 2] - df[, 1]),
    alternative = "two.sided"
  )$mean_diff[2:3],
  list(NA_real_, NA_real_),
  check.attributes = FALSE
)

#-------------------------------------------------------------------------------
# CI level
#-------------------------------------------------------------------------------
# Dependent two sample
expect_equal(
  t_test_paired(
    data = list(value1 = df[, 1], value2 = df[, 2]),
    alternative = "two.sided",
    ci_level = 0.8
  )$mean_diff[2:3] |>
    unlist(),
  t.test(
    x = df[, 2],
    y = df[, 1],
    alternative = "two.sided",
    paired = TRUE,
    conf.level = 0.8
  )$conf.int,
  check.attributes = FALSE,
  tolerance = 0.0001,
  scale = 1
)

# Dependent one sample
expect_equal(
  t_test_paired(
    data = list(value1 = df[, 2] - df[, 1]),
    alternative = "two.sided",
    ci_level = 0.8
  )$mean_diff[2:3] |>
    unlist(),
  t.test(
    x = df[, 2] - df[, 1],
    alternative = "two.sided",
    conf.level = 0.8
  )$conf.int,
  check.attributes = FALSE,
  tolerance = 0.0001,
  scale = 1
)

#-------------------------------------------------------------------------------
# True mean/difference
#-------------------------------------------------------------------------------
# Dependent two sample
expect_equal(
  t_test_paired(
    data = list(value1 = df[, 1], value2 = df[, 2]),
    alternative = "two.sided",
    ci_level = 0.95,
    mean_null = 0.5
  )$p,
  t.test(
    x = df[, 2],
    y = df[, 1],
    alternative = "two.sided",
    paired = TRUE,
    mu = 0.5
  )$p.value,
  tolerance = 0.0001,
  scale = 1
)

expect_equal(
  t_test_paired(
    data = list(value1 = df[, 1], value2 = df[, 2]),
    alternative = "greater",
    ci_level = 0.95,
    mean_null = 0.5
  )$p,
  t.test(
    x = df[, 2],
    y = df[, 1],
    alternative = "greater",
    paired = TRUE,
    mu = 0.5
  )$p.value,
  tolerance = 0.0001,
  scale = 1
)

expect_equal(
  t_test_paired(
    data = list(value1 = df[, 1], value2 = df[, 2]),
    alternative = "less",
    ci_level = 0.95,
    mean_null = 0.5
  )$p,
  t.test(
    x = df[, 2],
    y = df[, 1],
    alternative = "less",
    paired = TRUE,
    mu = 0.5
  )$p.value,
  tolerance = 0.0001,
  scale = 1
)
