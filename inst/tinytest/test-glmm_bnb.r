library(tinytest)

#-------------------------------------------------------------------------------
# Data
#-------------------------------------------------------------------------------
set.seed(1234)
d <- sim_bnb(
  n = 40,
  mean1 = 10,
  ratio = 1.5,
  dispersion = 2
)

lrt <- glmm_bnb(d, test = "lrt", ci_level = 0.95)
wald <- glmm_bnb(d, test = "wald", ci_level = 0.95)

#-------------------------------------------------------------------------------
# Structure
#-------------------------------------------------------------------------------
expect_equal(
  lrt,
  list(chisq = 17.5773165609801, df = 1L, p = 2.7585931727773e-05,
    ratio = list(estimate = 1.45762368911064, lower = 1.24458390253481,
        upper = 1.71115203153423), mean1 = list(estimate = 6.92417304255506,
        lower = 5.0925254629797, upper = 9.27811111143692), mean2 = list(
        estimate = 10.0928386543296, lower = 6.33807521447318,
        upper = 15.8762586771356), dispersion = list(estimate = 43.900202301203,
        lower = 5.26076300503589, upper = 366.339970122528),
    item_sd = list(estimate = 0.84413064875369, lower = 0.647631172737355,
        upper = 1.10025054716491), n1 = 40L, n2 = 40L,
    method = "GLMM for bivariate negative binomial ratio of means",
    test = "lrt", alternative = "two.sided", ci_level = 0.95,
    hessian = "Hessian appears to be positive definite.",
    convergence = "relative convergence (4)"),
  tolerance = 0.0001,
  scale = 1
)

expect_equal(
  wald,
  list(chisq = 23.1477529148312, df = 1L, p = 1.5001841082508e-06,
    ratio = list(estimate = 1.45762368911064, lower = 1.25020281771652,
        upper = 1.69945771113938), mean1 = list(estimate = 6.92417304255506,
        lower = 5.1653376861414, upper = 9.28190473429849), mean2 = list(
        estimate = 10.0928386543296, lower = 6.45771972967132,
        upper = 15.7742045747647), dispersion = list(estimate = 43.900202301203,
        lower = 5.26076300503589, upper = 366.339970122528),
    item_sd = list(estimate = 0.84413064875369, lower = 0.647631172737355,
        upper = 1.10025054716491), n1 = 40L, n2 = 40L,
    method = "GLMM for bivariate negative binomial ratio of means",
    test = "wald", alternative = "two.sided", ci_level = 0.95,
    hessian = "Hessian appears to be positive definite.",
    convergence = "relative convergence (4)"),
  tolerance = 0.0001,
  scale = 1
)
