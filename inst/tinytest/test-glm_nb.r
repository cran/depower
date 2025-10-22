library(tinytest)

#-------------------------------------------------------------------------------
# Data
#-------------------------------------------------------------------------------
set.seed(1234)
d <- sim_nb(
  n1 = 60,
  n2 = 40,
  mean1 = 10,
  ratio = 1.5,
  dispersion1 = 2,
  dispersion2 = 8
)

lrt <- glm_nb(d, equal_dispersion = FALSE, test = "lrt", ci_level = 0.95)
wald <- glm_nb(d, equal_dispersion = FALSE, test = "wald", ci_level = 0.95)

#-------------------------------------------------------------------------------
# Structure
#-------------------------------------------------------------------------------
expect_equal(
  lrt,
  list(
    chisq = 10.0246446613575,
    df = 1L,
    p = 0.00154459471894568,
    ratio = list(
      estimate = 1.54293426717091,
      lower = 1.18982361092849,
      upper = 1.9827851367772
    ),
    mean1 = list(
      estimate = 9.31667005453285,
      lower = 7.49376459638661,
      upper = 11.7239913902618
    ),
    mean2 = list(
      estimate = 14.3750094830638,
      lower = 12.70001,
      upper = 16.27092
    ),
    dispersion1 = list(
      estimate = 1.5454223071454,
      lower = 1.00531585918218,
      upper = 2.35860005631339
    ),
    dispersion2 = list(
      estimate = 11.0799848708247,
      lower = 4.980081,
      upper = 24.65142
    ),
    n1 = 60L,
    n2 = 40L,
    method = "GLM for independent negative binomial ratio of means",
    test = "lrt",
    alternative = "two.sided",
    equal_dispersion = FALSE,
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
    chisq = 11.3515984221669,
    df = 1L,
    p = 0.000753830618317084,
    ratio = list(
      estimate = 1.54293426717091,
      lower = 1.19889335038147,
      upper = 1.9857030252547
    ),
    mean1 = list(
      estimate = 9.31667005453285,
      lower = 7.4784965450777,
      upper = 11.6066565494585
    ),
    mean2 = list(
      estimate = 14.3750094830638,
      lower = 12.70001,
      upper = 16.27092
    ),
    dispersion1 = list(
      estimate = 1.5454223071454,
      lower = 1.01028646457562,
      upper = 2.36401277376892
    ),
    dispersion2 = list(
      estimate = 11.0799848708247,
      lower = 4.980081,
      upper = 24.65142
    ),
    n1 = 60L,
    n2 = 40L,
    method = "GLM for independent negative binomial ratio of means",
    test = "wald",
    alternative = "two.sided",
    equal_dispersion = FALSE,
    ci_level = 0.95,
    hessian = "Hessian appears to be positive definite.",
    convergence = "relative convergence (4)"
  ),
  tolerance = 0.001,
  scale = 1
)
