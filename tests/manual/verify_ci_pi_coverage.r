#-------------------------------------------------------------------------------
# Results
#-------------------------------------------------------------------------------
# =======================================================================
#   Coverage Probability Verification
# =======================================================================
#   Replicates: 600
# Simulations per replicate: 200
# Simulations for true power: 30000
# Nominal CI level: 0.95
# Nominal PI level: 0.95
# =======================================================================
#
#   SCENARIO 1: Welch's t-test (log-lognormal data)
# -----------------------------------------------------------------------
# Parameters: n1 = 20 , n2 = 20 , ratio = 1.3 , cv1 = 0.4 , cv2 = 0.4
#
# Estimating true power with 30000 simulations...
# True power estimate: 0.5538
#
# Running 600 replicate studies...
#   Completed 100 of 600 replicates
#   Completed 200 of 600 replicates
#   Completed 300 of 600 replicates
#   Completed 400 of 600 replicates
#   Completed 500 of 600 replicates
#   Completed 600 of 600 replicates
#
# Results for Welch's t-test:
#   CI coverage (Wilson):  0.9633  (nominal: 0.95 )
# CI coverage (Exact):   0.9633  (nominal: 0.95 )
# PI coverage:           0.9683  (nominal: 0.95 )
#
# SCENARIO 2: Wald test (negative binomial data)
# -----------------------------------------------------------------------
#   Parameters: n1 = 15 , n2 = 15 , mean1 = 10 , ratio = 1.5 , dispersion = 2
#
# Estimating true power with 30000 simulations...
# Using partial cluster of size 1
# True power estimate: 0.3373
#
# Running 600 replicate studies...
# Completed 100 of 600 replicates
# Completed 200 of 600 replicates
# Completed 300 of 600 replicates
# Completed 400 of 600 replicates
# Completed 500 of 600 replicates
# Completed 600 of 600 replicates
#
# Results for Wald test (NB):
#   CI coverage (Wilson):  0.95  (nominal: 0.95 )
# CI coverage (Exact):   0.955  (nominal: 0.95 )
# PI coverage:           0.97  (nominal: 0.95 )
#
# SCENARIO 3: Approximate parametric Wald test (negative binomial data)
# -----------------------------------------------------------------------
#   Parameters: n1 = 15 , n2 = 15 , mean1 = 10 , ratio = 1.5 , dispersion = 2
#
# Estimating true power with 30000 simulations...
# Using partial cluster of size 1
# Using partial cluster of size 1
# True power estimate: 0.2827
#
# Running 600 replicate studies...
# Completed 100 of 600 replicates
# Completed 200 of 600 replicates
# Completed 300 of 600 replicates
# Completed 400 of 600 replicates
# Completed 500 of 600 replicates
# Completed 600 of 600 replicates
#
# Results for Wald test (NB):
#   CI coverage (Wilson):  0.9417  (nominal: 0.95 )
# CI coverage (Exact):   0.9533  (nominal: 0.95 )
# PI coverage:           0.945  (nominal: 0.95 )
#
# =======================================================================
#   SUMMARY
# =======================================================================
#
#   Nominal coverage level: 0.95
#
# Confidence Interval Coverage (true power in CI):
#   Welch t-test (Wilson): 0.9633
# Welch t-test (Exact):  0.9633
# Wald NB test (Wilson): 0.95
# Wald NB test (Exact):  0.955
# Approx Parametric Wald NB test (Wilson): 0.9417
# Approx Parametric Wald NB test (Exact):  0.9533
#
# Prediction Interval Coverage (future estimate in PI):
#   Welch t-test:          0.9683
# Wald NB test:          0.97
# Approx Parametric Wald NB test:          0.945
#
# Interpretation:
#   - CI coverage should be close to 0.95
# - PI coverage should be close to 0.95
# - Exact CI tends to be conservative (coverage > nominal)
# - Wilson CI provides coverage closer to nominal
#
# Monte Carlo SE of coverage estimates: 0.0089
# =======================================================================

#-------------------------------------------------------------------------------
# Verify Coverage Probability for Confidence Intervals and Prediction Intervals
#-------------------------------------------------------------------------------
# This script verifies that:
# 1. Confidence intervals achieve nominal coverage of the true power
# 2. Prediction intervals achieve nominal coverage of future power estimates
#
# For each scenario, we:
# 1. Estimate "true" power using a very large number of simulations
# 2. Repeatedly:
#    a. Run a smaller simulation study to get a power estimate
#    b. Compute CI and PI
#    c. Check if true power falls in CI (for CI coverage)
#    d. Run another small simulation and check if it falls in PI (for PI coverage)
# 3. Report empirical coverage rates
#-------------------------------------------------------------------------------
library(depower)
set.seed(42)

#-------------------------------------------------------------------------------
# Configuration
#-------------------------------------------------------------------------------
n_replicates <- 600 # Number of replicate studies
nsims_small <- 200 # Simulations per replicate study
nsims_large <- 30000 # Simulations for "true" power estimate
ci_level <- 0.95
pi_level <- 0.95

cat("=======================================================================\n")
cat("Coverage Probability Verification\n")
cat("=======================================================================\n")
cat("Replicates:", n_replicates, "\n")
cat("Simulations per replicate:", nsims_small, "\n")
cat("Simulations for true power:", nsims_large, "\n")
cat("Nominal CI level:", ci_level, "\n")
cat("Nominal PI level:", pi_level, "\n")
cat(
  "=======================================================================\n\n"
)

#-------------------------------------------------------------------------------
# Scenario 1: Welch's t-test on log-lognormal data
#-------------------------------------------------------------------------------
cat("SCENARIO 1: Welch's t-test (log-lognormal data)\n")
cat("-----------------------------------------------------------------------\n")

# Parameters
n1_welch <- 20
n2_welch <- 20
ratio_welch <- 1.3
cv1_welch <- 0.4
cv2_welch <- 0.4

cat(
  "Parameters: n1 =",
  n1_welch,
  ", n2 =",
  n2_welch,
  ", ratio =",
  ratio_welch,
  ", cv1 =",
  cv1_welch,
  ", cv2 =",
  cv2_welch,
  "\n\n"
)

# Step 1: Estimate "true" power with large simulation
cat("Estimating true power with", nsims_large, "simulations...\n")
true_power_welch <- sim_log_lognormal(
  n1 = n1_welch,
  n2 = n2_welch,
  ratio = ratio_welch,
  cv1 = cv1_welch,
  cv2 = cv2_welch,
  cor = 0,
  nsims = nsims_large
) |>
  power(t_test_welch())

true_power_welch_value <- true_power_welch$power
cat("True power estimate:", round(true_power_welch_value, 4), "\n\n")

# Step 2: Run replicate studies
cat("Running", n_replicates, "replicate studies...\n")

ci_covers_true_wilson <- logical(n_replicates)
ci_covers_true_exact <- logical(n_replicates)
pi_covers_future <- logical(n_replicates)

for (r in seq_len(n_replicates)) {
  # Run a small simulation study
  power_result <- sim_log_lognormal(
    n1 = n1_welch,
    n2 = n2_welch,
    ratio = ratio_welch,
    cv1 = cv1_welch,
    cv2 = cv2_welch,
    cor = 0,
    nsims = nsims_small
  ) |>
    power(t_test_welch())

  # Compute CI
  ci_wilson <- add_power_ci(
    power_result,
    ci_level = ci_level,
    method = "wilson"
  )
  ci_exact <- add_power_ci(power_result, ci_level = ci_level, method = "exact")

  # Check CI coverage (does CI contain true power?)
  ci_covers_true_wilson[r] <- (true_power_welch_value >=
    ci_wilson$power_ci_lower) &
    (true_power_welch_value <= ci_wilson$power_ci_upper)
  ci_covers_true_exact[r] <- (true_power_welch_value >=
    ci_exact$power_ci_lower) &
    (true_power_welch_value <= ci_exact$power_ci_upper)

  # Compute PI
  pi_result <- add_power_pi(
    power_result,
    future_nsims = nsims_small,
    pi_level = pi_level
  )

  # Run another independent study to get future power estimate
  future_power <- sim_log_lognormal(
    n1 = n1_welch,
    n2 = n2_welch,
    ratio = ratio_welch,
    cv1 = cv1_welch,
    cv2 = cv2_welch,
    cor = 0,
    nsims = nsims_small
  ) |>
    power(t_test_welch())

  # Check PI coverage (does PI contain future estimate?)
  pi_covers_future[r] <- (future_power$power >= pi_result$power_pi_lower) &
    (future_power$power <= pi_result$power_pi_upper)

  if (r %% 100 == 0) {
    cat("  Completed", r, "of", n_replicates, "replicates\n")
  }
}

# Report results
cat("\nResults for Welch's t-test:\n")
cat(
  "  CI coverage (Wilson): ",
  round(mean(ci_covers_true_wilson), 4),
  " (nominal:",
  ci_level,
  ")\n"
)
cat(
  "  CI coverage (Exact):  ",
  round(mean(ci_covers_true_exact), 4),
  " (nominal:",
  ci_level,
  ")\n"
)
cat(
  "  PI coverage:          ",
  round(mean(pi_covers_future), 4),
  " (nominal:",
  pi_level,
  ")\n\n"
)

#-------------------------------------------------------------------------------
# Scenario 2: Wald test on negative binomial data
#-------------------------------------------------------------------------------
cat("SCENARIO 2: Wald test (negative binomial data)\n")
cat("-----------------------------------------------------------------------\n")

# Parameters
n1_nb <- 15
n2_nb <- 15
mean1_nb <- 10
ratio_nb <- 1.5
disp1_nb <- 2
disp2_nb <- 2

cat(
  "Parameters: n1 =",
  n1_nb,
  ", n2 =",
  n2_nb,
  ", mean1 =",
  mean1_nb,
  ", ratio =",
  ratio_nb,
  ", dispersion =",
  disp1_nb,
  "\n\n"
)

# Step 1: Estimate "true" power with large simulation
cat("Estimating true power with", nsims_large, "simulations...\n")
true_power_nb <- sim_nb(
  n1 = n1_nb,
  n2 = n2_nb,
  mean1 = mean1_nb,
  ratio = ratio_nb,
  dispersion1 = disp1_nb,
  dispersion2 = disp2_nb,
  nsims = nsims_large
) |>
  power(wald_test_nb(), ncores = 4)

true_power_nb_value <- true_power_nb$power
cat("True power estimate:", round(true_power_nb_value, 4), "\n\n")

# Step 2: Run replicate studies
cat("Running", n_replicates, "replicate studies...\n")

ci_covers_true_wilson_nb <- logical(n_replicates)
ci_covers_true_exact_nb <- logical(n_replicates)
pi_covers_future_nb <- logical(n_replicates)

for (r in seq_len(n_replicates)) {
  # Run a small simulation study
  power_result <- sim_nb(
    n1 = n1_nb,
    n2 = n2_nb,
    mean1 = mean1_nb,
    ratio = ratio_nb,
    dispersion1 = disp1_nb,
    dispersion2 = disp2_nb,
    nsims = nsims_small
  ) |>
    power(wald_test_nb())

  # Compute CI
  ci_wilson <- add_power_ci(
    power_result,
    ci_level = ci_level,
    method = "wilson"
  )
  ci_exact <- add_power_ci(power_result, ci_level = ci_level, method = "exact")

  # Check CI coverage (does CI contain true power?)
  ci_covers_true_wilson_nb[r] <- (true_power_nb_value >=
    ci_wilson$power_ci_lower) &
    (true_power_nb_value <= ci_wilson$power_ci_upper)
  ci_covers_true_exact_nb[r] <- (true_power_nb_value >=
    ci_exact$power_ci_lower) &
    (true_power_nb_value <= ci_exact$power_ci_upper)

  # Compute PI
  pi_result <- add_power_pi(
    power_result,
    future_nsims = nsims_small,
    pi_level = pi_level
  )

  # Run another independent study to get future power estimate
  future_power <- sim_nb(
    n1 = n1_nb,
    n2 = n2_nb,
    mean1 = mean1_nb,
    ratio = ratio_nb,
    dispersion1 = disp1_nb,
    dispersion2 = disp2_nb,
    nsims = nsims_small
  ) |>
    power(wald_test_nb())

  # Check PI coverage (does PI contain future estimate?)
  pi_covers_future_nb[r] <- (future_power$power >= pi_result$power_pi_lower) &
    (future_power$power <= pi_result$power_pi_upper)

  if (r %% 100 == 0) {
    cat("  Completed", r, "of", n_replicates, "replicates\n")
  }
}

# Report results
cat("\nResults for Wald test (NB):\n")
cat(
  "  CI coverage (Wilson): ",
  round(mean(ci_covers_true_wilson_nb), 4),
  " (nominal:",
  ci_level,
  ")\n"
)
cat(
  "  CI coverage (Exact):  ",
  round(mean(ci_covers_true_exact_nb), 4),
  " (nominal:",
  ci_level,
  ")\n"
)
cat(
  "  PI coverage:          ",
  round(mean(pi_covers_future_nb), 4),
  " (nominal:",
  pi_level,
  ")\n\n"
)

#-------------------------------------------------------------------------------
# Scenario 3: Approximate parametric Wald test on negative binomial data
#-------------------------------------------------------------------------------
cat("SCENARIO 3: Approximate parametric Wald test (negative binomial data)\n")
cat("-----------------------------------------------------------------------\n")

# Parameters
n1_nb2 <- 15
n2_nb2 <- 15
mean1_nb2 <- 10
ratio_nb2 <- 1.5
disp1_nb2 <- 2
disp2_nb2 <- 2

cat(
  "Parameters: n1 =",
  n1_nb2,
  ", n2 =",
  n2_nb2,
  ", mean1 =",
  mean1_nb2,
  ", ratio =",
  ratio_nb2,
  ", dispersion =",
  disp1_nb2,
  "\n\n"
)

# Step 1: Estimate "true" power with large simulation
cat("Estimating true power with", nsims_large, "simulations...\n")
true_power_nb2 <- sim_nb(
  n1 = n1_nb2,
  n2 = n2_nb2,
  mean1 = mean1_nb2,
  ratio = ratio_nb2,
  dispersion1 = disp1_nb2,
  dispersion2 = disp2_nb2,
  nsims = nsims_large
) |>
  power(wald_test_nb(distribution = simulated(nsims = 20000)), ncores = 4)

true_power_nb2_value <- true_power_nb2$power
cat("True power estimate:", round(true_power_nb2_value, 4), "\n\n")

# Step 2: Run replicate studies
cat("Running", n_replicates, "replicate studies...\n")

ci_covers_true_wilson_nb2 <- logical(n_replicates)
ci_covers_true_exact_nb2 <- logical(n_replicates)
pi_covers_future_nb2 <- logical(n_replicates)

for (r in seq_len(n_replicates)) {
  # Run a small simulation study
  power_result <- sim_nb(
    n1 = n1_nb2,
    n2 = n2_nb2,
    mean1 = mean1_nb2,
    ratio = ratio_nb2,
    dispersion1 = disp1_nb2,
    dispersion2 = disp2_nb2,
    nsims = nsims_small
  ) |>
    power(wald_test_nb(distribution = simulated(nsims = 4000)))

  # Compute CI
  ci_wilson <- add_power_ci(
    power_result,
    ci_level = ci_level,
    method = "wilson"
  )
  ci_exact <- add_power_ci(power_result, ci_level = ci_level, method = "exact")

  # Check CI coverage (does CI contain true power?)
  ci_covers_true_wilson_nb2[r] <- (true_power_nb2_value >=
    ci_wilson$power_ci_lower) &
    (true_power_nb2_value <= ci_wilson$power_ci_upper)
  ci_covers_true_exact_nb2[r] <- (true_power_nb2_value >=
    ci_exact$power_ci_lower) &
    (true_power_nb2_value <= ci_exact$power_ci_upper)

  # Compute PI
  pi_result <- add_power_pi(
    power_result,
    future_nsims = nsims_small,
    pi_level = pi_level
  )

  # Run another independent study to get future power estimate
  future_power <- sim_nb(
    n1 = n1_nb2,
    n2 = n2_nb2,
    mean1 = mean1_nb2,
    ratio = ratio_nb2,
    dispersion1 = disp1_nb2,
    dispersion2 = disp2_nb2,
    nsims = nsims_small
  ) |>
    power(wald_test_nb(distribution = simulated(nsims = 4000)))

  # Check PI coverage (does PI contain future estimate?)
  pi_covers_future_nb2[r] <- (future_power$power >= pi_result$power_pi_lower) &
    (future_power$power <= pi_result$power_pi_upper)

  if (r %% 100 == 0) {
    cat("  Completed", r, "of", n_replicates, "replicates\n")
  }
}

# Report results
cat("\nResults for Wald test (NB):\n")
cat(
  "  CI coverage (Wilson): ",
  round(mean(ci_covers_true_wilson_nb2), 4),
  " (nominal:",
  ci_level,
  ")\n"
)
cat(
  "  CI coverage (Exact):  ",
  round(mean(ci_covers_true_exact_nb2), 4),
  " (nominal:",
  ci_level,
  ")\n"
)
cat(
  "  PI coverage:          ",
  round(mean(pi_covers_future_nb2), 4),
  " (nominal:",
  pi_level,
  ")\n\n"
)

#-------------------------------------------------------------------------------
# Summary
#-------------------------------------------------------------------------------
cat("=======================================================================\n")
cat("SUMMARY\n")
cat("=======================================================================\n")
cat("\nNominal coverage level:", ci_level, "\n\n")

cat("Confidence Interval Coverage (true power in CI):\n")
cat("  Welch t-test (Wilson):", round(mean(ci_covers_true_wilson), 4), "\n")
cat("  Welch t-test (Exact): ", round(mean(ci_covers_true_exact), 4), "\n")
cat("  Wald NB test (Wilson):", round(mean(ci_covers_true_wilson_nb), 4), "\n")
cat("  Wald NB test (Exact): ", round(mean(ci_covers_true_exact_nb), 4), "\n")
cat(
  "  Approx Parametric Wald NB test (Wilson):",
  round(mean(ci_covers_true_wilson_nb2), 4),
  "\n"
)
cat(
  "  Approx Parametric Wald NB test (Exact): ",
  round(mean(ci_covers_true_exact_nb2), 4),
  "\n"
)

cat("\nPrediction Interval Coverage (future estimate in PI):\n")
cat("  Welch t-test:         ", round(mean(pi_covers_future), 4), "\n")
cat("  Wald NB test:         ", round(mean(pi_covers_future_nb), 4), "\n")
cat(
  "  Approx Parametric Wald NB test:         ",
  round(mean(pi_covers_future_nb2), 4),
  "\n"
)

cat("\nInterpretation:\n")
cat("  - CI coverage should be close to", ci_level, "\n")
cat("  - PI coverage should be close to", pi_level, "\n")
cat("  - Exact CI tends to be conservative (coverage > nominal)\n")
cat("  - Wilson CI provides coverage closer to nominal\n")

cat(
  "\nMonte Carlo SE of coverage estimates:",
  round(sqrt(0.95 * 0.05 / n_replicates), 4),
  "\n"
)
cat("=======================================================================\n")
