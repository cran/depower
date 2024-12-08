#-------------------------------------------------------------------------------
# Two sample dependent t_test
#-------------------------------------------------------------------------------
grid <- expand.grid(
  ln1 = 1:3,
  ln2 = 1:3,
  lratio = 1:3,
  lcv1 = 1:3,
  lcv2 = 1:3,
  lcor = 1:3
)

grid <- grid |>
  dplyr::rowwise() |>
  dplyr::mutate(
    params = list(list(
      n1 = sample(10:30, ln1, replace = FALSE),
      n2 = sample(10:30, ln2, replace = FALSE),
      ratio = runif(n = lratio, min = 1.2, max = 1.8),
      cv1 = runif(n = lcv1, min = 0.2, max = 0.7),
      cv2 = runif(n = lcv2, min = 0.2, max = 0.7),
      cor = runif(n = lcor, min = 0.1, max = 0.5)
    ))
  )

test_plot <- function(fun, grid) {
  for(i in sample(seq_len(nrow(grid)), size = 10, replace = FALSE)) {#seq_len(nrow(grid))) {
    res <- fun(
      n1 = grid$params[[i]]$n1,
      n2 = grid$params[[i]]$n1,
      ratio = grid$params[[i]]$ratio,
      cv1 = grid$params[[i]]$cv1,
      cv2 = grid$params[[i]]$cv2,
      cor = grid$params[[i]]$cor,
      nsims = 300
    )

    print(res |> power() |> plot())
  }
}

test_plot(fun = sim_log_lognormal, grid = grid)

sim_log_lognormal(
  n1 = 10:20,
  n2 = 10:20,
  ratio = c(1.7, 2, 2.5),
  cv1 = 0.5,
  cv2 = 0.5,
  cor = c(0.3, 0.5),
  nsims = 300
) |>
  power() |>
  plot()

sim_log_lognormal(
  n1 = 10:20,
  n2 = 10:20,
  ratio = c(1.7, 2, 2.5),
  cv1 = 0.5,
  cv2 = 0.5,
  cor = c(0, 0.3, 0.5),
  nsims = 300
) |>
  power() |>
  plot(facet_row = "ratio")

#-------------------------------------------------------------------------------
# Two sample independent t_test
#-------------------------------------------------------------------------------
grid <- expand.grid(
  ln1 = 1:3,
  ln2 = 1:3,
  lratio = 1:3,
  lcv1 = 1:3,
  lcv2 = 1:3
)

grid <- grid |>
  dplyr::rowwise() |>
  dplyr::mutate(
    params = list(list(
      n1 = sample(10:30, ln1, replace = FALSE),
      n2 = sample(10:30, ln2, replace = FALSE),
      ratio = runif(n = lratio, min = 1.2, max = 1.8),
      cv1 = runif(n = lcv1, min = 0.2, max = 0.7),
      cv2 = runif(n = lcv2, min = 0.2, max = 0.7),
      cor = 0
    ))
  )

test_plot <- function(fun, grid) {
  for(i in sample(seq_len(nrow(grid)), size = 10, replace = FALSE)) {#seq_len(nrow(grid))) {
    res <- fun(
      n1 = grid$params[[i]]$n1,
      n2 = grid$params[[i]]$n2,
      ratio = grid$params[[i]]$ratio,
      cv1 = grid$params[[i]]$cv1,
      cv2 = grid$params[[i]]$cv2,
      cor = grid$params[[i]]$cor,
      nsims = 300
    )

    print(res |> power() |> plot())
  }
}

test_plot(fun = sim_log_lognormal, grid = grid)


sim_log_lognormal(
  n1 = 10:20,
  n2 = 10:20,
  ratio = c(1.7, 2, 2.5),
  cv1 = 0.5,
  cv2 = 0.5,
  cor = c(0),
  nsims = 300
) |>
  power() |>
  plot()

sim_log_lognormal(
  n1 = 10:20,
  n2 = 10:20,
  ratio = c(1.7, 2, 2.5),
  cv1 = 0.5,
  cv2 = 0.5,
  cor = c(0),
  nsims = 300
) |>
  power() |>
  plot(facet_row = "ratio")

#-------------------------------------------------------------------------------
# One sample t_test
#-------------------------------------------------------------------------------
grid <- expand.grid(
  ln1 = 1:3,
  lratio = 1:3,
  lcv1 = 1:3
)

grid <- grid |>
  dplyr::rowwise() |>
  dplyr::mutate(
    params = list(list(
      n1 = sample(10:30, ln1, replace = FALSE),
      ratio = runif(n = lratio, min = 1.2, max = 1.8),
      cv1 = runif(n = lcv1, min = 0.2, max = 0.7)
    ))
  )

test_plot <- function(fun, grid) {
  for(i in sample(seq_len(nrow(grid)), size = 10, replace = FALSE)) {#seq_len(nrow(grid))) {
    res <- fun(
      n1 = grid$params[[i]]$n1,
      ratio = grid$params[[i]]$ratio,
      cv1 = grid$params[[i]]$cv1,
      nsims = 300
    )

    print(res |> power() |> plot())
  }
}

test_plot(fun = sim_log_lognormal, grid = grid)


sim_log_lognormal(
  n1 = 10:20,
  ratio = c(1.7, 2, 2.5),
  cv1 = 0.5,
  nsims = 300
) |>
  power() |>
  plot()


sim_log_lognormal(
  n1 = 10:20,
  ratio = c(1.7, 2, 2.5),
  cv1 = 0.5,
  nsims = 300
) |>
  power() |>
  plot(facet_row = "ratio")

#-------------------------------------------------------------------------------
# Two sample independent nb test
#-------------------------------------------------------------------------------
grid <- expand.grid(
  ln1 = 1:3,
  ln2 = 1:3,
  lmean1 = 1:3,
  lratio = 1:3,
  ldispersion1 = 1:3,
  ldispersion2 = 1:3,
  ltest = 1:2
)

grid <- grid |>
  dplyr::rowwise() |>
  dplyr::mutate(
    params = list(list(
      n1 = sample(10:30, ln1, replace = FALSE),
      n2 = sample(10:30, ln2, replace = FALSE),
      mean1 = rnorm(n = lmean1, mean = 20, sd = 5),
      ratio = runif(n = lratio, min = 1.2, max = 1.8),
      dispersion1 = runif(n = ldispersion1, min = 1, max = 8),
      dispersion2 = runif(n = ldispersion2, min = 1, max = 8),
      test = sample(c(`Wald test` = "wald_test_nb", `LRT` = "lrt_nb"), ltest, replace = FALSE)
    ))
  )

test_plot <- function(fun, grid) {
  for(i in sample(seq_len(nrow(grid)), size = 5, replace = FALSE)) {#seq_len(nrow(grid))) {
    res <- fun(
      n1 = grid$params[[i]]$n1,
      n2 = grid$params[[i]]$n2,
      mean1 = grid$params[[i]]$mean1,
      ratio = grid$params[[i]]$ratio,
      dispersion1 = grid$params[[i]]$dispersion1,
      dispersion2 = grid$params[[i]]$dispersion2,
      nsims = 75
    )

    print(res |> power() |> plot())
  }
}

test_plot(fun = sim_nb, grid = grid)

sim_nb(
  n1 = c(15, 20),
  n2 = c(15, 20),
  mean1 = c(10, 20),
  ratio = c(1.7, 2),
  dispersion1 = 2,
  dispersion2 = 2,
  nsims = 100
) |>
  power() |>
  plot()

sim_nb(
  n1 = c(15, 20),
  n2 = c(15, 20),
  mean1 = c(10, 20),
  ratio = c(1.7, 2),
  dispersion1 = 2,
  dispersion2 = 2,
  nsims = 100
) |>
  power() |>
  plot(facet_row = "ratio")

#-------------------------------------------------------------------------------
# dependent bnb test
#-------------------------------------------------------------------------------
grid <- expand.grid(
  ln = 1:3,
  lmean1 = 1:3,
  lratio = 1:3,
  ldispersion = 1:3,
  ltest = 1:2
)

grid <- grid |>
  dplyr::rowwise() |>
  dplyr::mutate(
    params = list(list(
      n = sample(10:30, ln, replace = FALSE),
      mean1 = rnorm(n = lmean1, mean = 20, sd = 5),
      ratio = runif(n = lratio, min = 1.2, max = 1.8),
      dispersion = runif(n = ldispersion, min = 1, max = 8),
      test = sample(c(`Wald test` = "wald_test_bnb", `LRT` = "lrt_bnb"), ltest, replace = FALSE)
    ))
  )

test_plot <- function(fun, grid) {
  for(i in sample(seq_len(nrow(grid)), size = 5, replace = FALSE)) {#seq_len(nrow(grid))) {
    res <- fun(
      n = grid$params[[i]]$n,
      mean1 = grid$params[[i]]$mean1,
      ratio = grid$params[[i]]$ratio,
      dispersion = grid$params[[i]]$dispersion,
      nsims = 75
    )

    print(res |> power() |> plot())
  }
}

test_plot(fun = sim_bnb, grid = grid)

sim_bnb(
  n = c(15, 20),
  mean1 = c(10, 20),
  ratio = c(1.2, 1.4),
  dispersion = 10,
  nsims = 500
) |>
  power() |>
  plot()

sim_bnb(
  n = c(15, 20),
  mean1 = c(10, 20),
  ratio = c(1.2, 1.4),
  dispersion = 10,
  nsims = 500
) |>
  power() |>
  plot(facet_row = "ratio")

#-------------------------------------------------------------------------------
# When nsims varies, make sure all plot points are shown and caption is correct.
#-------------------------------------------------------------------------------
set.seed(1234)
sim_bnb(
  n = c(10),
  mean1 = 5.9,
  ratio = c(0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2),
  dispersion = c(0.49, 2),
  nsims = 100,
  return_type = "list"
) |>
  power("Wald Test" = wald_test_bnb(link = "squared"), alpha = 0.05) |>
  plot(hline = 0.05, x_axis = "ratio")

#-------------------------------------------------------------------------------
# If you jump the queue in axis sorting, are those who got pushed back still
# plotted?
# n1 and n2 should be on plot.
#-------------------------------------------------------------------------------
set.seed(1234)
sim_nb(
  n1 = c(15, 20),
  n2 = c(15, 20),
  mean1 = c(10, 20),
  ratio = c(1.7, 2),
  dispersion1 = 2,
  dispersion2 = 2,
  nsims = 100
) |>
  power() |>
  plot(facet_row = "ratio", x_axis = "n1")
