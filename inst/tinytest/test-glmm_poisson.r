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

lrt <- glmm_poisson(d, test = "lrt", ci_level = 0.95)
wald <- glmm_poisson(d, test = "wald", ci_level = 0.95)

#-------------------------------------------------------------------------------
# Structure
#-------------------------------------------------------------------------------
expect_equal(
  lrt,
  list(
    chisq = 32.525582707165,
    df = 1L,
    p = 1.17633805049949e-08,
    ratio = list(
      estimate = 1.45336781526169,
      lower = 1.27717807844155,
      upper = 1.65542865749175
    ),
    mean1 = list(
      estimate = 6.89552391371006,
      lower = 5.08790784218065,
      upper = 9.20682992879782
    ),
    mean2 = list(
      estimate = 10.0217325255535,
      lower = 7.543392,
      upper = 13.31432
    ),
    item_sd = list(
      estimate = 0.850882520374344,
      lower = 0.655998178980542,
      upper = 1.10366322144939
    ),
    n1 = 40L,
    n2 = 40L,
    method = "GLMM for two dependent Poisson ratio of means",
    test = "lrt",
    alternative = "two.sided",
    ci_level = 0.95,
    hessian = "Hessian appears to be positive definite.",
    convergence = "relative convergence (4)"
  ),
  tolerance = 0.001,
  scale = 1
)

expect_equal(
  wald,
  list(
    chisq = 31.9648555673617,
    df = 1L,
    p = 1.5698720626369e-08,
    ratio = list(
      estimate = 1.45336781526169,
      lower = 1.27668983158746,
      upper = 1.65449583303417
    ),
    mean1 = list(
      estimate = 6.89552391371006,
      lower = 5.16226787476225,
      upper = 9.2107289273006
    ),
    mean2 = list(
      estimate = 10.0217325255535,
      lower = 7.543392,
      upper = 13.31432
    ),
    item_sd = list(
      estimate = 0.850882520374344,
      lower = 0.655998178980542,
      upper = 1.10366322144939
    ),
    n1 = 40L,
    n2 = 40L,
    method = "GLMM for two dependent Poisson ratio of means",
    test = "wald",
    alternative = "two.sided",
    ci_level = 0.95,
    hessian = "Hessian appears to be positive definite.",
    convergence = "relative convergence (4)"
  ),
  tolerance = 0.001,
  scale = 1
)
