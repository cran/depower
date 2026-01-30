# depower 2026.1.30

### New Features

- Add functions `add_power_ci()` and `add_power_pi()` for confidence intervals and prediction intervals of the power estimate, respectively.
  - They add uncertainty intervals for power to the object return by `power()`.
  - Function `plot.depower()` will add shaded regions for the intervals.
- Add functions `eval_power_ci()` and `eval_power_pi()` for confidence intervals and prediction intervals of the power estimate, respectively.
  - In contrast to the `add_power_*()` variants, these are for planning power simulation studies.
  - Assess uncertainty in hypothetical simulation-based power estimates before committing computational time to a simulation run.

### Other Updates

- General improvements in documentation.

# depower 2025.10.21

### New Features

- Add functions `asymptotic()` and `simulated()` for specification of type of
  test.
- Add argument `distribution` to functions:
    - `wald_test_nb()`
    - `wald_test_bnb()`
    - `lrt_nb()`
    - `lrt_bnb()`
    
  for control of running an asymptotic, randomization, or parametric simulation test.

### Bug Fixes

- Fix bug where `equal_dispersion` argument in function `lrt_nb()` was not passed downstream and instead fixed to `FALSE`.
- Fix calculation of confidence limits for `mean2` and `dispersion2` in functions `glm_nb()`, `glmm_bnb()`, and `glmm_poisson()`.

### Other Updates

- Add and update tolerances in tests.
- General improvements in documentation.
- General improvements in code.

# depower 2025.1.20

- Fix typos in documentation.
- Fix test error resulting from ATLAS numerical differences.

# depower 2024.12.4

First release to CRAN.

### New Features

- Power
    - `power()`
- Plot
    - `plot.depower()`
- Simulate
    - `sim_log_lognormal()`
    - `sim_nb()`
    - `sim_bnb()`
- Test
    - `t_test_welch()`
    - `t_test_paired()`
    - `wald_test_nb()`
    - `wald_test_bnb()`
    - `lrt_nb()`
    - `lrt_bnb()`
    - `glm_nb()`
    - `glmm_bnb()`
    - `glmm_poisson()`
- MLE
    - `mle_nb_null()` and `mle_nb_alt()`
    - `mle_bnb_null()` and `mle_bnb_alt()`
- Log-likelihood
    - `nll_nb_null()` and `nll_nb_alt()`
    - `nll_bnb_null()` and `nll_bnb_alt()`

# depower 2024.8.1

Proof of concept.
