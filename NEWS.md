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
