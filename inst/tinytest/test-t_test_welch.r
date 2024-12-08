library(tinytest)

#-------------------------------------------------------------------------------
# Data
#-------------------------------------------------------------------------------
set.seed(20231201)
df <- MASS::mvrnorm(n = 30, mu = c(2, 3), Sigma = matrix(c(2, 0, 0, 2), 2L, 2L))

#-------------------------------------------------------------------------------
# Structure
#-------------------------------------------------------------------------------
expect_equal(
  t_test_welch(
    data = list(value1 = df[,1], value2 = df[,2]),
    alternative = "two.sided",
    ci_level = 0.95
  ),
  list(
    t = 2.22469325299248,
    df = 49.9197931267991,
    p = 0.0306517980644275,
    diff_mean = list(
      estimate = 0.819712860095389,
      lower = 0.0796075540992306,
      upper = 1.55981816609155),
    mean1 = 1.97225676976233,
    mean2 = 2.79196962985772,
    n1 = 30L,
    n2 = 30L,
    method = "Welch's two-sample t-test for 'two-sided' alternative",
    alternative = "two.sided",
    ci_level = 0.95,
    mean_null = 0
  ),
  tolerance = 0.0001,
  scale = 1
)

#-------------------------------------------------------------------------------
# Test expected p-value
#-------------------------------------------------------------------------------
# Independent two sample
expect_equal(
  t_test_welch(
    data = list(value1 = df[,1], value2 = df[,2]),
    alternative = "two.sided"
  )$p,
  t.test(
    x = df[,1],
    y = df[,2],
    alternative = "two.sided",
    paired = FALSE,
    var.equal = FALSE
  )$p.value
)

expect_equal(
  t_test_welch(
    data = list(value1 = df[,2], value2 = df[,1]),
    alternative = "greater"
  )$p,
  t.test(
    x = df[,1],
    y = df[,2],
    alternative = "greater",
    paired = FALSE,
    var.equal = FALSE
  )$p.value
)

expect_equal(
  t_test_welch(
    data = list(value1 = df[,2], value2 = df[,1]),
    alternative = "less"
  )$p,
  t.test(
    x = df[,1],
    y = df[,2],
    alternative = "less",
    paired = FALSE,
    var.equal = FALSE
  )$p.value
)

#-------------------------------------------------------------------------------
# Test expected CI
#-------------------------------------------------------------------------------
# Independent two sample
expect_equal(
  t_test_welch(
    data = list(value1 = df[,1], value2 = df[,2]),
    alternative = "two.sided",
    ci_level = 0.95
  )$diff_mean$lower,
  t.test(
    x = df[,2],
    y = df[,1],
    alternative = "two.sided",
    paired = FALSE,
    var.equal = FALSE
  )$conf.int[1]
)

expect_equal(
  t_test_welch(
    data = list(value1 = df[,1], value2 = df[,2]),
    alternative = "two.sided",
    ci_level = 0.95
  )$diff_mean$upper,
  t.test(
    x = df[,2],
    y = df[,1],
    alternative = "two.sided",
    paired = FALSE,
    var.equal = FALSE
  )$conf.int[2]
)

expect_equal(
  t_test_welch(
    data = list(value1 = df[,1], value2 = df[,2]),
    alternative = "greater",
    ci_level = 0.95
  )$diff_mean$lower,
  t.test(
    x = df[,2],
    y = df[,1],
    alternative = "greater",
    paired = FALSE,
    var.equal = FALSE
  )$conf.int[1]
)

expect_equal(
  t_test_welch(
    data = list(value1 = df[,1], value2 = df[,2]),
    alternative = "greater",
    ci_level = 0.95
  )$diff_mean$upper,
  t.test(
    x = df[,2],
    y = df[,1],
    alternative = "greater",
    paired = FALSE,
    var.equal = FALSE
  )$conf.int[2]
)

expect_equal(
  t_test_welch(
    data = list(value1 = df[,1], value2 = df[,2]),
    alternative = "less",
    ci_level = 0.95
  )$diff_mean$lower,
  t.test(
    x = df[,2],
    y = df[,1],
    alternative = "less",
    paired = FALSE,
    var.equal = FALSE
  )$conf.int[1]
)

expect_equal(
  t_test_welch(
    data = list(value1 = df[,1], value2 = df[,2]),
    alternative = "less",
    ci_level = 0.95
  )$diff_mean$upper,
  t.test(
    x = df[,2],
    y = df[,1],
    alternative = "less",
    paired = FALSE,
    var.equal = FALSE
  )$conf.int[2]
)

#-------------------------------------------------------------------------------
# Test absence of CI
#-------------------------------------------------------------------------------
# Independent two sample
expect_equal(
  t_test_welch(
    data = list(value1 = df[,1], value2 = df[,2]),
    alternative = "less"
  )$diff_mean[2:3],
  list(NA_real_, NA_real_),
  check.attributes = FALSE
)

#-------------------------------------------------------------------------------
# CI level
#-------------------------------------------------------------------------------
# Independent two sample
expect_equal(
  t_test_welch(
    data = list(value1 = df[,1], value2 = df[,2]),
    alternative = "two.sided",
    ci_level = 0.8
  )$diff_mean[2:3] |> unlist(),
  t.test(
    x = df[,2],
    y = df[,1],
    alternative = "two.sided",
    var.equal = FALSE,
    conf.level = 0.8
  )$conf.int,
  check.attributes = FALSE
)

#-------------------------------------------------------------------------------
# True mean/difference
#-------------------------------------------------------------------------------
# Independent two sample
expect_equal(
  t_test_welch(
    data = list(value1 = df[,1], value2 = df[,2]),
    alternative = "two.sided",
    ci_level = 0.95,
    mean_null = 0.5
  )$p,
  t.test(
    x = df[,2],
    y = df[,1],
    alternative = "two.sided",
    var.equal = FALSE,
    mu = 0.5
  )$p.value
)

expect_equal(
  t_test_welch(
    data = list(value1 = df[,1], value2 = df[,2]),
    alternative = "greater",
    ci_level = 0.95,
    mean_null = 0.5
  )$p,
  t.test(
    x = df[,2],
    y = df[,1],
    alternative = "greater",
    var.equal = FALSE,
    mu = 0.5
  )$p.value
)

expect_equal(
  t_test_welch(
    data = list(value1 = df[,1], value2 = df[,2]),
    alternative = "less",
    ci_level = 0.95,
    mean_null = 0.5
  )$p,
  t.test(
    x = df[,2],
    y = df[,1],
    alternative = "less",
    var.equal = FALSE,
    mu = 0.5
  )$p.value
)
