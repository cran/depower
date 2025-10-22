#' Simulate data from a normal distribution
#'
#' Simulate data from the log transformed lognormal distribution (i.e. a normal
#' distribution). This function handles all three cases:
#' 1. One-sample data
#' 2. Dependent two-sample data
#' 3. Independent two-sample data
#'
#' Based on assumed characteristics of the original lognormal distribution, data
#' is simulated from the corresponding log-transformed (normal) distribution.
#' This simulated data is suitable for assessing power of a hypothesis for the
#' geometric mean or ratio of geometric means from the original lognormal data.
#'
#' This method can also be useful for other population distributions which are
#' positive and where it makes sense to describe the ratio of geometric means.
#' However, the lognormal distribution is theoretically correct in the sense
#' that you can log transform to a normal distribution, compute the summary
#' statistic, then apply the inverse transformation to summarize on the original
#' lognormal scale.
#'
#' Let \eqn{GM(\cdot)} be the geometric mean and \eqn{AM(\cdot)} be the
#' arithmetic mean. For independent lognormal samples \eqn{X_1} and \eqn{X_2}
#'
#' \deqn{\text{Fold Change} = \frac{GM(X_2)}{GM(X_1)}}
#'
#' For dependent lognormal samples \eqn{X_1} and \eqn{X_2}
#'
#' \deqn{\text{Fold Change} = GM\left( \frac{X_2}{X_1} \right)}
#'
#' Unlike ratios and the arithmetic mean, for equal sample sizes of
#' \eqn{X_1} and \eqn{X_2} it follows that
#' \eqn{\frac{GM(X_2)}{GM(X_1)} = GM \left( \frac{X_2}{X_1} \right) =
#' e^{AM(\ln X_2) - AM(\ln X_1)} = e^{AM(\ln X_2 - \ln X_1)}}.
#'
#' The coefficient of variation (CV) for \eqn{X} is defined as
#'
#' \deqn{CV = \frac{SD(X)}{AM(X)}}
#'
#' The relationship between sample statistics for the original lognormal data
#' (\eqn{X}) and the natural logged data (\eqn{\ln{X}}) are
#'
#' \deqn{
#' \begin{aligned}
#' AM(X)  &= e^{AM(\ln{X}) + \frac{Var(\ln{X})}{2}} \\
#' GM(X)  &= e^{AM(\ln{X})} \\
#' Var(X) &= AM(X)^2 \left( e^{Var(\ln{X})} - 1 \right) \\
#' CV(X)  &= \frac{\sqrt{AM(X)^2 \left( e^{Var(\ln{X})} - 1 \right)}}{AM(X)} \\
#'        &= \sqrt{e^{Var(\ln{X})} - 1}
#' \end{aligned}
#' }
#'
#' and
#'
#' \deqn{
#' \begin{aligned}
#' AM(\ln{X})  &= \ln \left( \frac{AM(X)}{\sqrt{CV(X)^2 + 1}} \right) \\
#' Var(\ln{X}) &= \ln(CV(X)^2 + 1) \\
#' Cor(\ln{X_1}, \ln{X_2}) &= \frac{\ln \left( Cor(X_1, X_2)CV(X_1)CV(X_2) + 1 \right)}{SD(\ln{X_1})SD(\ln{X_2})}
#' \end{aligned}
#' }
#'
#' Based on the properties of correlation and variance,
#'
#' \deqn{
#' \begin{aligned}
#' Var(X_2 - X_1) &= Var(X_1) + Var(X_2) - 2Cov(X_1, X_2) \\
#'            &= Var(X_1) + Var(X_2) - 2Cor(X_1, X_2)SD(X_1)SD(X_2) \\
#' SD(X_2 - X_1)  &= \sqrt{Var(X_2 - X_1)}
#' \end{aligned}
#' }
#'
#' The standard deviation of the differences gets smaller the more positive
#' the correlation and conversely gets larger the more negative the correlation.
#' For the special case where the two samples are uncorrelated and each has the
#' same variance, it follows that
#'
#' \deqn{
#' \begin{aligned}
#' Var(X_2 - X_1) &= \sigma^2 + \sigma^2 \\
#' SD(X_2 - X_1) &= \sqrt{2}\sigma
#' \end{aligned}
#' }
#'
#' @references
#' \insertRef{julious_2004}{depower}
#'
#' \insertRef{hauschke_1992}{depower}
#'
#' \insertRef{johnson_1994}{depower}
#'
#' @param n1 (integer: `[2, Inf)`)\cr
#'        The sample size(s) of sample 1.
#' @param n2 (integer: `NULL`; `[2, Inf)`)\cr
#'        The sample size(s) of sample 2. Set as `NULL` if you want to simulate
#'        for the one-sample case.
#' @param ratio (numeric: `(0, Inf)`)\cr
#'        The assumed population fold change(s) of sample 2 with respect to
#'        sample 1.
#' - For one-sample data, `ratio` is defined as the geometric mean (GM) of the
#'   original lognormal population distribution.
#' - For dependent two-sample data, `ratio` is defined by
#'   GM(sample 2 / sample 1) of the original lognormal population distributions.
#'     - e.g. `ratio = 2` assumes that the geometric mean of all paired
#'       ratios (sample 2 / sample 1) is 2.
#' - For independent two-sample data, the `ratio` is defined by
#'   GM(group 2) / GM(group 1) of the original lognormal population distributions.
#'     - e.g. `ratio = 2` assumes that the geometric mean of sample 2 is 2
#'       times larger than the geometric mean of sample 1.
#'
#' See 'Details' for additional information.
#' @param cv1 (numeric: `(0, Inf)`)\cr
#'        The coefficient of variation(s) of sample 1 in the original lognormal
#'        data.
#' @param cv2 (numeric: `NULL`; `(0, Inf)`)\cr
#'        The coefficient of variation(s) of sample 2 in the original lognormal
#'        data. Set as `NULL` if you want to simulate for the one-sample case.
#' @param cor (numeric: `0`; `[-1, 1]`)\cr
#'        The correlation(s) between sample 1 and sample 2 in the original
#'        lognormal data. Not used for the one-sample case.
#' @param nsims (Scalar integer: `1L`; `[1,Inf)`)\cr
#'        The number of simulated data sets. If `nsims > 1`, the data is
#'        returned in a list-column of a depower simulation data frame.
#' @param return_type (string: `"list"`; `c("list", "data.frame")`)\cr
#'        The data structure of the simulated data. If `"list"` (default), a
#'        list object is returned. If `"data.frame"` a data frame in tall format
#'        is returned. The list object provides computational efficiency and the
#'        data frame object is convenient for formulas. See 'Value'.
#' @param messages (Scalar logical: `TRUE`)\cr
#'        Whether or not to display messages for pathological simulation cases.
#'
#' @return If `nsims = 1` and the number of unique parameter combinations is
#' one, the following objects are returned:
#' - If one-sample data with `return_type = "list"`, a list:
#' \tabular{lll}{
#'   Slot \tab Name \tab Description \cr
#'   1 \tab \tab One sample of simulated normal values.
#' }
#' - If one-sample data with `return_type = "data.frame"`, a data frame:
#' \tabular{lll}{
#'   Column \tab Name \tab Description \cr
#'   1 \tab `item` \tab Pair/subject/item indicator. \cr
#'   2 \tab `value` \tab Simulated normal values.
#' }
#' - If two-sample data with `return_type = "list"`, a list:
#' \tabular{lll}{
#'   Slot \tab Name \tab Description \cr
#'   1 \tab \tab Simulated normal values from sample 1. \cr
#'   2 \tab \tab Simulated normal values from sample 2.
#' }
#' - If two-sample data with `return_type = "data.frame"`, a data frame:
#' \tabular{lll}{
#'   Column \tab Name \tab Description \cr
#'   1 \tab `item` \tab Pair/subject/item indicator. \cr
#'   2 \tab `condition` \tab Time/group/condition indicator. \cr
#'   3 \tab `value` \tab Simulated normal values.
#' }
#'
#' If `nsims > 1` or the number of unique parameter combinations is greater than
#' one, each object described above is returned in data frame, located in a
#' list-column named `data`.
#' - If one-sample data, a data frame:
#' \tabular{lll}{
#'   Column \tab Name \tab Description \cr
#'   1 \tab `n1` \tab The sample size. \cr
#'   2 \tab `ratio` \tab Geometric mean \[GM(sample 1)\]. \cr
#'   3 \tab `cv1` \tab Coefficient of variation for sample 1. \cr
#'   4 \tab `nsims` \tab Number of data simulations. \cr
#'   5 \tab `distribution` \tab Distribution sampled from. \cr
#'   6 \tab `data` \tab List-column of simulated data.
#' }
#' - If two-sample data, a data frame:
#' \tabular{lll}{
#'   Column \tab Name \tab Description \cr
#'   1 \tab `n1` \tab Sample size of sample 1. \cr
#'   2 \tab `n2` \tab Sample size of sample 2. \cr
#'   3 \tab `ratio` \tab Ratio of geometric means \[GM(sample 2) / GM(sample 1)\]
#'                       or geometric mean ratio \[GM(sample 2 / sample 1)\]. \cr
#'   4 \tab `cv1` \tab Coefficient of variation for sample 1. \cr
#'   5 \tab `cv2` \tab Coefficient of variation for sample 2. \cr
#'   6 \tab `cor` \tab Correlation between samples. \cr
#'   7 \tab `nsims` \tab Number of data simulations. \cr
#'   8 \tab `distribution` \tab Distribution sampled from. \cr
#'   9 \tab `data` \tab List-column of simulated data.
#' }
#'
#' @examples
#' #----------------------------------------------------------------------------
#' # sim_log_lognormal() examples
#' #----------------------------------------------------------------------------
#' library(depower)
#'
#' # Independent two-sample data returned in a data frame
#' sim_log_lognormal(
#'   n1 = 10,
#'   n2 = 10,
#'   ratio = 1.3,
#'   cv1 = 0.35,
#'   cv2 = 0.35,
#'   cor = 0,
#'   nsims = 1,
#'   return_type = "data.frame"
#' )
#'
#' # Independent two-sample data returned in a list
#' sim_log_lognormal(
#'   n1 = 10,
#'   n2 = 10,
#'   ratio = 1.3,
#'   cv1 = 0.35,
#'   cv2 = 0.35,
#'   cor = 0,
#'   nsims = 1,
#'   return_type = "list"
#' )
#'
#' # Dependent two-sample data returned in a data frame
#' sim_log_lognormal(
#'   n1 = 10,
#'   n2 = 10,
#'   ratio = 1.3,
#'   cv1 = 0.35,
#'   cv2 = 0.35,
#'   cor = 0.4,
#'   nsims = 1,
#'   return_type = "data.frame"
#' )
#'
#' # Dependent two-sample data returned in a list
#' sim_log_lognormal(
#'   n1 = 10,
#'   n2 = 10,
#'   ratio = 1.3,
#'   cv1 = 0.35,
#'   cv2 = 0.35,
#'   cor = 0.4,
#'   nsims = 1,
#'   return_type = "list"
#' )
#'
#' # One-sample data returned in a data frame
#' sim_log_lognormal(
#'   n1 = 10,
#'   ratio = 1.3,
#'   cv1 = 0.35,
#'   nsims = 1,
#'   return_type = "data.frame"
#' )
#'
#' # One-sample data returned in a list
#' sim_log_lognormal(
#'   n1 = 10,
#'   ratio = 1.3,
#'   cv1 = 0.35,
#'   nsims = 1,
#'   return_type = "list"
#' )
#'
#' # Independent two-sample data: two simulations for four parameter combinations.
#' # Returned as a list-column of lists within a data frame
#' sim_log_lognormal(
#'   n1 = c(10, 20),
#'   n2 = c(10, 20),
#'   ratio = 1.3,
#'   cv1 = 0.35,
#'   cv2 = 0.35,
#'   cor = 0,
#'   nsims = 2,
#'   return_type = "list"
#' )
#'
#' # Dependent two-sample data: two simulations for two parameter combinations.
#' # Returned as a list-column of lists within a data frame
#' sim_log_lognormal(
#'   n1 = c(10, 20),
#'   n2 = c(10, 20),
#'   ratio = 1.3,
#'   cv1 = 0.35,
#'   cv2 = 0.35,
#'   cor = 0.4,
#'   nsims = 2,
#'   return_type = "list"
#' )
#'
#' # One-sample data: two simulations for two parameter combinations
#' # Returned as a list-column of lists within a data frame
#' sim_log_lognormal(
#'   n1 = c(10, 20),
#'   ratio = 1.3,
#'   cv1 = 0.35,
#'   nsims = 2,
#'   return_type = "list"
#' )
#'
#' @importFrom mvnfast rmvn
#' @importFrom stats rnorm
#'
#' @export
sim_log_lognormal <- function(
  n1,
  n2 = NULL,
  ratio,
  cv1,
  cv2 = NULL,
  cor = 0,
  nsims = 1L,
  return_type = "list",
  messages = TRUE
) {
  #-----------------------------------------------------------------------------
  # Check arguments
  #-----------------------------------------------------------------------------
  if (!is.numeric(n1) || any(n1 < 2L)) {
    stop("Argument 'n1' must be an integer vector from [2, Inf).")
  }
  if (!is.null(n2) && (!is.numeric(n2) || any(n2 < 2L))) {
    stop("Argument 'n2' must be an integer vector from [2, Inf).")
  }
  if (!is.numeric(ratio) || any(ratio <= 0)) {
    stop("Argument 'ratio' must be a positive numeric vector.")
  }
  if (!is.numeric(cv1) || any(cv1 <= 0)) {
    stop("Argument 'cv1' must be a positive numeric vector.")
  }
  if (!is.null(cv2) && (!is.numeric(cv2) || any(cv2 <= 0))) {
    stop("Argument 'cv2' must be a positive numeric vector.")
  }
  if (!is.numeric(cor) || any(cor < -1) || any(cor > 1)) {
    stop("Argument 'cor' must be a numeric vector from [-1, 1].")
  }
  if (length(return_type) != 1L) {
    stop("Argument 'return_type' must be one of 'list' or 'data.frame'.")
  }
  if (!is.numeric(nsims) || length(nsims) != 1L || nsims < 1L) {
    stop("Argument 'nsims' must be a scalar integer from [1, Inf).")
  }
  data.frame <- switch(
    return_type,
    "list" = FALSE,
    "data.frame" = TRUE,
    stop("Argument 'return_type' must be one of 'list' or 'data.frame'.")
  )

  if ((is.null(n2) && !is.null(cv2)) || (!is.null(n2) && is.null(cv2))) {
    stop(
      "Arguments 'n2' and 'cv2' must both be NULL (one sample) or both numeric (two samples)."
    )
  }

  needs_grid <- any(c(nsims, lengths(list(n1, n2, ratio, cv1, cv2, cor))) > 1L)

  #-----------------------------------------------------------------------------
  # Simulate data
  #-----------------------------------------------------------------------------
  res <- if (is.null(n2)) {
    # n2 = NULL implies one sample case
    if (needs_grid) {
      grid_log_lognormal_one_sample(
        n1 = n1,
        ratio = ratio,
        cv1 = cv1,
        nsims = nsims,
        data.frame = data.frame
      )
    } else {
      sim_log_lognormal_one_sample(
        n1 = n1,
        ratio = ratio,
        cv1 = cv1,
        nsims = nsims,
        data.frame = data.frame
      )
    }
  } else {
    # two sample case
    if (needs_grid) {
      grid_log_lognormal_two_sample(
        n1 = n1,
        n2 = n2,
        ratio = ratio,
        cv1 = cv1,
        cv2 = cv2,
        cor = cor,
        nsims = nsims,
        data.frame = data.frame,
        messages = messages
      )
    } else {
      sim_log_lognormal_two_sample(
        n1 = n1,
        n2 = n2,
        ratio = ratio,
        cv1 = cv1,
        cv2 = cv2,
        cor = cor,
        nsims = nsims,
        data.frame = data.frame
      )
    }
  }

  #-----------------------------------------------------------------------------
  # Return
  #-----------------------------------------------------------------------------
  res
}

grid_log_lognormal_one_sample <- function(
  n1,
  ratio,
  cv1,
  nsims,
  data.frame
) {
  #-----------------------------------------------------------------------------
  # Unique combinations for simulating data.
  #-----------------------------------------------------------------------------
  grid_sim <- expand.grid(
    n1 = n1,
    ratio = ratio,
    cv1 = cv1,
    nsims = nsims,
    distribution = "One-sample log(lognormal)",
    stringsAsFactors = FALSE
  )

  #-----------------------------------------------------------------------------
  # Simulate data
  #-----------------------------------------------------------------------------
  res <- grid_sim |>
    dplyr::rowwise() |>
    dplyr::mutate(
      data = list(
        sim_log_lognormal_one_sample(
          n1 = n1,
          ratio = ratio,
          cv1 = cv1,
          nsims = nsims,
          data.frame = data.frame
        )
      )
    ) |>
    dplyr::ungroup()

  #-----------------------------------------------------------------------------
  # Prepare return
  #-----------------------------------------------------------------------------
  # Column labels
  vars <- c(
    "n Pairs" = "n1",
    "Ratio" = "ratio",
    "CV" = "cv1",
    "N Simulations" = "nsims",
    "Distribution" = "distribution",
    "Data" = "data"
  )
  idx <- match(names(res), vars)
  if (anyNA(idx)) {
    stop("Unknown variable found while labeling data frame.")
  }
  for (i in seq_len(ncol(res))) {
    attr(res[[i]], "label") <- names(vars)[idx][i]
  }

  # Class
  class(res) <- c("depower", "log_lognormal_one_sample", class(res))

  #-----------------------------------------------------------------------------
  # Return
  #-----------------------------------------------------------------------------
  res
}

grid_log_lognormal_two_sample <- function(
  n1,
  n2,
  ratio,
  cv1,
  cv2,
  cor,
  nsims,
  data.frame,
  messages
) {
  #-----------------------------------------------------------------------------
  # Unique combinations for simulating data
  #-----------------------------------------------------------------------------
  grid_sim <- expand.grid(
    n1 = n1,
    n2 = n2,
    ratio = ratio,
    cv1 = cv1,
    cv2 = cv2,
    cor = cor,
    nsims = nsims,
    stringsAsFactors = FALSE
  )

  # Handle case of differing sample sizes for correlated data.
  bad_rows <- which(
    grid_sim[["cor"]] != 0 & (grid_sim[["n1"]] != grid_sim[["n2"]])
  )
  if (length(bad_rows) > 0L) {
    # Message isn't necessary if you only have correlated data
    if (any(grid_sim[["cor"]] == 0) && any(grid_sim[["cor"]] != 0)) {
      message <- paste0(
        "Message from depower::sim_log_lognormal():\n",
        "Arguments 'n1' and 'n2' must be the same for correlated data.\n",
        length(bad_rows),
        " rows were removed because they had different sample sizes.\n",
        nrow(grid_sim) - length(bad_rows),
        " rows (valid parameter combinations) remain."
      )
      if (messages) {
        message(message)
      }
    }
    # This must be below/after/under the message
    grid_sim <- grid_sim |>
      dplyr::filter(!(cor != 0 & (n1 != n2)))
  }

  grid_sim[["distribution"]] <- ifelse(
    grid_sim[["cor"]] == 0,
    "Independent two-sample log(lognormal)",
    "Dependent two-sample log(lognormal)"
  )

  #-----------------------------------------------------------------------------
  # Simulate data
  #-----------------------------------------------------------------------------
  res <- grid_sim |>
    dplyr::rowwise() |>
    dplyr::mutate(
      data = list(
        sim_log_lognormal_two_sample(
          n1 = n1,
          n2 = n2,
          ratio = ratio,
          cv1 = cv1,
          cv2 = cv2,
          cor = cor,
          nsims = nsims,
          data.frame = data.frame
        )
      )
    ) |>
    dplyr::ungroup()

  #-----------------------------------------------------------------------------
  # Prepare return
  #-----------------------------------------------------------------------------
  # Column labels
  vars <- c(
    "n1" = "n1",
    "n2" = "n2",
    "Ratio" = "ratio",
    "CV1" = "cv1",
    "CV2" = "cv2",
    "Correlation" = "cor",
    "N Simulations" = "nsims",
    "Distribution" = "distribution",
    "Data" = "data"
  )
  idx <- match(names(res), vars)
  if (anyNA(idx)) {
    stop("Unknown variable found while labeling data frame.")
  }
  for (i in seq_len(ncol(res))) {
    attr(res[[i]], "label") <- names(vars)[idx][i]
  }

  # Class
  cls <- if (all(res$cor == 0L)) {
    "log_lognormal_independent_two_sample"
  } else if (all(res$cor != 0)) {
    "log_lognormal_dependent_two_sample"
  } else {
    "log_lognormal_mixed_two_sample"
  }
  class(res) <- c("depower", cls, class(res))

  #-----------------------------------------------------------------------------
  # Return
  #-----------------------------------------------------------------------------
  res
}

sim_log_lognormal_one_sample <- function(
  n1,
  ratio,
  cv1,
  nsims = 1L,
  data.frame = FALSE
) {
  #-----------------------------------------------------------------------------
  # Simulate
  #-----------------------------------------------------------------------------
  log_ratio <- log(ratio)
  log_sigma <- sqrt(log(cv1^2 + 1))

  if (nsims > 1L) {
    res <- lapply(
      seq_len(nsims),
      function(x) {
        list(value1 = rnorm(n = n1, mean = log_ratio, sd = log_sigma))
      }
    )
    if (data.frame) {
      res <- lapply(res, list_to_df)
    }
  } else {
    res <- list(value1 = rnorm(n = n1, mean = log_ratio, sd = log_sigma))
    if (data.frame) {
      res <- list_to_df(res)
    }
  }

  #-----------------------------------------------------------------------------
  # Return
  #-----------------------------------------------------------------------------
  res
}

sim_log_lognormal_two_sample <- function(
  n1,
  n2 = NULL,
  ratio,
  cv1,
  cv2 = NULL,
  cor = 0,
  nsims = 1L,
  data.frame = FALSE
) {
  #-----------------------------------------------------------------------------
  # Simulate
  #-----------------------------------------------------------------------------
  log_ratio <- log(ratio)

  is_paired <- cor != 0L
  is_unequal_ss <- n1 != n2
  if (is_unequal_ss & is_paired) {
    stop("Arguments 'n1' and 'n2' must be the same for dependent data.")
  }
  maxn <- max(n1, n2)
  is_n1_smaller <- n1 < n2

  # cor * cv1 * cv2 + 1 must be > 0
  lower_cor <- -1 / (cv1 * cv2)
  if (cor < lower_cor) {
    stop(sprintf(
      "Infeasible 'cor' for given CVs: require cor >= %.4f but got %.4f.",
      lower_cor,
      cor
    ))
  }

  log_sigma1 <- sqrt(log(cv1^2 + 1))
  log_sigma2 <- sqrt(log(cv2^2 + 1))
  log_cor <- log(cor * cv1 * cv2 + 1) / (log_sigma1 * log_sigma2)
  log_cormat <- matrix(c(1, log_cor, log_cor, 1), nrow = 2L, ncol = 2L)
  log_varmat <- matrix(
    c(
      log_sigma1^2,
      log_sigma1 * log_sigma2,
      log_sigma1 * log_sigma2,
      log_sigma2^2
    ),
    nrow = 2L,
    ncol = 2L
  )
  log_covmat <- log_cormat * log_varmat

  # use pre-allocated matrix for mvnfast::rmvn()
  A <- matrix(data = NA_real_, nrow = maxn, ncol = 2L)

  if (nsims > 1) {
    res <- lapply(
      seq_len(nsims),
      function(x) {
        rnorm_two_sample(
          maxn = maxn,
          log_ratio = log_ratio,
          log_covmat = log_covmat,
          is_unequal_ss = is_unequal_ss,
          is_n1_smaller = is_n1_smaller,
          n1 = n1,
          n2 = n2,
          A = A
        )
      }
    )
    if (data.frame) {
      res <- lapply(res, list_to_df)
    }
  } else {
    # nsims == 1
    res <- rnorm_two_sample(
      maxn = maxn,
      log_ratio = log_ratio,
      log_covmat = log_covmat,
      is_unequal_ss = is_unequal_ss,
      is_n1_smaller = is_n1_smaller,
      n1 = n1,
      n2 = n2,
      A = A
    )
    if (data.frame) {
      res <- list_to_df(res)
    }
  }

  #-----------------------------------------------------------------------------
  # Return
  #-----------------------------------------------------------------------------
  res
}

rnorm_two_sample <- function(
  maxn,
  log_ratio,
  log_covmat,
  is_unequal_ss,
  is_n1_smaller,
  n1,
  n2,
  A
) {
  rmvn(
    n = maxn,
    mu = c(0, log_ratio),
    sigma = log_covmat,
    A = A
  )

  if (is_unequal_ss) {
    # two sample data with unequal sample sizes
    if (is_n1_smaller) {
      # Return
      list(value1 = A[, 1L][seq_len(n1)], value2 = A[, 2L])
    } else {
      # Return
      list(value1 = A[, 1L], value2 = A[, 2L][seq_len(n2)])
    }
  } else {
    # two sample data with equal sample size
    # Return
    list(value1 = A[, 1L], value2 = A[, 2L])
  }
}
