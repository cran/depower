#-------------------------------------------------------------------------------
# Overhead of power()? There, but acceptable.
#-------------------------------------------------------------------------------
set.seed(1234)
data <- sim_log_lognormal(
  n1 = 20,
  n2 = 20,
  ratio = c(1.2),
  cv1 = 0.4,
  cv2 = 0.4,
  cor = 0,
  nsims = 100000
)

microbenchmark::microbenchmark(
  lapply = lapply(data$data[[1]], function(x) t_test_welch(x)),
  power = power(data),
  times = 20
)

# Unit: seconds
# expr   min      lq       mean     median   uq       max      neval cld
# lapply 1.284657 1.442017 1.467022 1.488866 1.523565 1.536142 20    a
# power  1.613477 1.715695 1.779597 1.803967 1.838617 1.861551 20    b
