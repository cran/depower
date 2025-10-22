#' Simulate data from a BNB distribution
#'
#' Simulate data from the bivariate negative binomial (BNB) distribution. The
#' BNB distribution is used to simulate count data where the event counts are
#' jointly dependent (correlated). For independent data, see
#' [depower::sim_nb()].
#'
#' The negative binomial distribution may be defined using a gamma-Poisson
#' mixture distribution. In this case, the Poisson parameter \eqn{\lambda}
#' is a random variable with gamma distribution. Equivalence between different
#' parameterizations are demonstrated below:
#'
#' ```{r, gamma_poisson_mixture, dev='png', dpi = 160, fig.width=4, fig.height=4, fig.show='hide'}
#' # Define constants and their relationships
#' n <- 10000
#' dispersion <- 8
#' mu <- 4
#' p <- dispersion / (dispersion + mu)
#' q <- mu / (mu + dispersion)
#' variance <- mu + (mu^2 / dispersion)
#' rate <- p / (1 - p)
#' scale <- (1 - p) / p
#'
#' # alternative formula for mu
#' mu_alt <- (dispersion * (1 - p)) / p
#' stopifnot(isTRUE(all.equal(mu, mu_alt)))
#'
#' set.seed(20240321)
#'
#' # Using built-in rnbinom with dispersion and mean
#' w <- rnbinom(n = n, size = dispersion, mu = mu)
#'
#' # Using gamma-Poisson mixture with gamma rate parameter
#' x <- rpois(
#'   n = n,
#'   lambda = rgamma(n = n, shape = dispersion, rate = rate)
#' )
#'
#' # Using gamma-Poisson mixture with gamma scale parameter
#' y <- rpois(
#'   n = n,
#'   lambda = rgamma(n = n, shape = dispersion, scale = scale)
#' )
#'
#' # Using gamma-Poisson mixture with multiplicative mean and
#' # gamma scale parameter
#' z <- rpois(
#'   n = n,
#'   lambda = mu * rgamma(n = n, shape = dispersion, scale = 1/dispersion)
#' )
#'
#' # Compare CDFs
#' par(mar=c(4,4,1,1))
#' plot(
#'   x = sort(w),
#'   y = (1:n)/n,
#'   xlim = range(c(w,x,y,z)),
#'   ylim = c(0,1),
#'   col = 'green',
#'   lwd = 4,
#'   type = 'l',
#'   main = 'CDF'
#' )
#' lines(x = sort(x), y = (1:n)/n, col = 'red', lwd = 2)
#' lines(x = sort(y), y = (1:n)/n, col = 'yellow', lwd = 1.5)
#' lines(x = sort(z), y = (1:n)/n, col = 'black')
#' ```
#'
#' \if{html}{\out{<div style="display: flex; justify-content: center; padding-top: 10px; padding-bottom: 10px;">
#'   <img style="max-width: 400px; width: 100\%; height: auto;" src="figures/gamma_poisson_mixture-1.png" alt="Gamma-Poisson Mixture CDF" />
#' </div>
#' <br>}}
#' \if{latex}{
#'   \out{\begin{center}}
#'   \figure{gamma_poisson_mixture-1.png}{options: width=3in}
#'   \out{\end{center}}
#' }
#'
#' The BNB distribution is implemented by compounding two conditionally
#' independent Poisson random variables \eqn{X_1 \mid G = g \sim \text{Poisson}(\mu g)}
#' and \eqn{X_2 \mid G = g \sim \text{Poisson}(r \mu g)} with a gamma random
#' variable \eqn{G \sim \text{Gamma}(\theta, \theta^{-1})}. The probability mass
#' function for the joint distribution of \eqn{X_1, X_2} is
#'
#' \deqn{
#' P(X_1 = x_1, X_2 = x_2) = \frac{\Gamma(x_1 + x_2 + \theta)}{(\mu + r \mu + \theta)^{x_1 + x_2 + \theta}}
#'   \frac{\mu^{x_1}}{x_1!} \frac{(r \mu)^{x_2}}{x_2!}
#'   \frac{\theta^{\theta}}{\Gamma(\theta)}
#' }
#'
#' where \eqn{x_1,x_2 \in \mathbb{N}^{\geq 0}} are specific values of the count
#' outcomes, \eqn{\theta \in \mathbb{R}^{> 0}} is the `dispersion` parameter
#' which controls the dispersion and level of correlation between the two
#' samples (otherwise known as the shape parameter of the gamma distribution),
#' \eqn{\mu \in \mathbb{R}^{> 0}} is the mean parameter, and
#' \eqn{r = \frac{\mu_2}{\mu_1} \in \mathbb{R}^{> 0}} is the `ratio` parameter
#' representing the multiplicative change in the mean of the second sample
#' relative to the first sample. \eqn{G} denotes the random subject effect and
#' the gamma distribution scale parameter is assumed to be the inverse of the
#' dispersion parameter (\eqn{\theta^{-1}}) for identifiability.
#'
#' Correlation decreases from 1 to 0 as the dispersion parameter increases from
#' 0 to infinity. For a given dispersion, increasing means also increases the
#' correlation. See 'Examples' for a demonstration.
#'
#' See 'Details' in [depower::sim_nb()] for additional information on the
#' negative binomial distribution.
#'
#' @references
#' \insertRef{yu_2020}{depower}
#'
#' \insertRef{rettiganti_2012}{depower}
#'
#' \insertRef{aban_2009}{depower}
#'
#' @param n (integer: `[2, Inf)`)\cr
#'        The number(s) of paired observations.
#' @param mean1 (numeric: `(0, Inf)`)\cr
#'        The mean(s) of sample 1 \eqn{(\mu_1)}.
#' @param mean2,ratio (numeric: `(0, Inf)`)\cr
#'        Only specify one of these arguments.
#' - `mean2`: The mean(s) of sample 2 \eqn{(\mu_2)}.
#' - `ratio`: The ratio(s) of means for sample 2 with respect to sample 1
#'        \eqn{\left( r = \frac{\mu_2}{\mu_1} \right)}.
#'
#' `mean2 = ratio * mean1`
#' @param dispersion (numeric: `(0, Inf)`)\cr
#'        The gamma distribution shape parameter(s) \eqn{(\theta)} which control
#'        the dispersion and the correlation between sample 1 and 2.
#'        See 'Details' and 'Examples'.
#' @param nsims (Scalar integer: `1L`; `[1, Inf)`)\cr
#'        The expected number of simulated data sets. If `nsims > 1`, the data
#'        is returned in a list-column of a depower simulation data frame.
#'        `nsims` may be reduced depending on `max_zeros`.
#' @param return_type (string: `"list"`; `c("list", "data.frame")`)\cr
#'        The data structure of the simulated data. If `"list"` (default), a
#'        list object is returned. If `"data.frame"` a data frame in tall format
#'        is returned. The list object provides computational efficiency and the
#'        data frame object is convenient for formulas. See 'Value'.
#' @param max_zeros (Scalar numeric: `0.99`; `[0, 1]`)\cr
#'        The maximum proportion of zeros each group in a simulated dataset is
#'        allowed to have. If the proportion of zeros is greater than this
#'        value, the corresponding data is excluded from the set of simulations.
#'        This is most likely to occur when the sample size is small and the
#'        dispersion parameter is small.
#'
#' @return If `nsims = 1` and the number of unique parameter combinations is
#' one, the following objects are returned:
#' - If `return_type = "list"`, a list:
#' \tabular{lll}{
#'   Slot \tab Name \tab Description \cr
#'   1 \tab \tab Simulated counts from sample 1. \cr
#'   2 \tab \tab Simulated counts from sample 2.
#' }
#' - If `return_type = "data.frame"`, a data frame:
#' \tabular{lll}{
#'   Column \tab Name \tab Description \cr
#'   1 \tab `item` \tab Subject/item indicator. \cr
#'   2 \tab `condition` \tab Sample/condition indicator. \cr
#'   3 \tab `value` \tab Simulated counts.
#' }
#'
#' If `nsims > 1` or the number of unique parameter combinations is greater than
#' one, each object described above is returned in a list-column named `data` in
#' a depower simulation data frame:
#' \tabular{lll}{
#'   Column \tab Name \tab Description \cr
#'   1 \tab `n1` \tab Sample size of sample 1. \cr
#'   2 \tab `n2` \tab Sample size of sample 2. \cr
#'   3 \tab `mean1` \tab Mean for sample 1. \cr
#'   4 \tab `mean2` \tab Mean for sample 2. \cr
#'   5 \tab `ratio` \tab Ratio of means (sample 2 / sample 1). \cr
#'   6 \tab `dispersion` \tab Gamma distribution shape parameter (dispersion). \cr
#'   7 \tab `nsims` \tab Number of valid simulation iterations. \cr
#'   8 \tab `distribution` \tab Distribution sampled from. \cr
#'   9 \tab `data` \tab List-column of simulated data.
#' }
#'
#' @seealso [depower::sim_nb()]
#'
#' @examples
#' #----------------------------------------------------------------------------
#' # sim_bnb() examples
#' #----------------------------------------------------------------------------
#' library(depower)
#'
#' # Paired two-sample data returned in a data frame
#' sim_bnb(
#'   n = 10,
#'   mean1 = 10,
#'   ratio = 1.6,
#'   dispersion = 3,
#'   nsims = 1,
#'   return_type = "data.frame"
#' )
#'
#' # Paired two-sample data returned in a list
#' sim_bnb(
#'   n = 10,
#'   mean1 = 10,
#'   ratio = 1.6,
#'   dispersion = 3,
#'   nsims = 1,
#'   return_type = "list"
#' )
#'
#' # Two simulations of paired two-sample data
#' # returned as a list of data frames
#' sim_bnb(
#'   n = 10,
#'   mean1 = 10,
#'   ratio = 1.6,
#'   dispersion = 3,
#'   nsims = 2,
#'   return_type = "data.frame"
#' )
#'
#' # Two simulations of Paired two-sample data
#' # returned as a list of lists
#' sim_bnb(
#'   n = 10,
#'   mean1 = 10,
#'   ratio = 1.6,
#'   dispersion = 3,
#'   nsims = 2,
#'   return_type = "list"
#' )
#'
#' #----------------------------------------------------------------------------
#' # Visualization of the BNB distribution as dispersion varies.
#' #----------------------------------------------------------------------------
#' set.seed(1234)
#' data <- lapply(
#'   X = c(1, 10, 100, 1000),
#'   FUN = function(x) {
#'     d <- sim_bnb(
#'       n = 1000,
#'       mean1 = 10,
#'       ratio = 1.5,
#'       dispersion = x,
#'       nsims = 1,
#'       return_type = "data.frame"
#'     )
#'     cor <- cor(
#'       x = d[d$condition == "1", ]$value,
#'       y = d[d$condition == "2", ]$value
#'     )
#'     cbind(dispersion = x, correlation = cor, d)
#'   }
#' )
#'
#' data <- do.call(what = "rbind", args = data)
#'
#' ggplot2::ggplot(
#'   data = data,
#'   mapping = ggplot2::aes(x = value, fill = condition)
#' ) +
#'   ggplot2::facet_wrap(
#'     facets = ggplot2::vars(.data$dispersion),
#'     ncol = 2,
#'     labeller = ggplot2::labeller(.rows = ggplot2::label_both)
#'   ) +
#'   ggplot2::geom_density(alpha = 0.3) +
#'   ggplot2::coord_cartesian(xlim = c(0, 60)) +
#'   ggplot2::geom_text(
#'     mapping = ggplot2::aes(
#'       x = 30,
#'       y = 0.12,
#'       label = paste0("Correlation: ", round(correlation, 2))
#'     ),
#'     check_overlap = TRUE
#'   ) +
#'   ggplot2::labs(
#'     x = "Value",
#'     y = "Density",
#'     fill = "Condition",
#'     caption = "Mean1=10, Mean2=15, ratio=1.5"
#'   )
#'
#' @importFrom stats rgamma rpois
#'
#' @export
sim_bnb <- function(
  n,
  mean1,
  mean2,
  ratio,
  dispersion,
  nsims = 1L,
  return_type = "list",
  max_zeros = 0.99
) {
  #-----------------------------------------------------------------------------
  # Check arguments
  #-----------------------------------------------------------------------------
  if (!is.numeric(n) || any(n < 2L)) {
    stop("Argument 'n' must be an integer vector from [2, Inf).")
  }
  if (!is.numeric(mean1) || any(mean1 <= 0)) {
    stop("Argument 'mean1' must be a positive numeric vector.")
  }

  missing_mean2 <- missing(mean2)
  missing_ratio <- missing(ratio)
  if (missing_mean2 && missing_ratio) {
    stop("You must specify one of the arguments: 'mean2' or 'ratio'.")
  }
  if (!missing_mean2 && !missing_ratio) {
    stop(
      "Arguments 'mean2' and 'ratio' were both specified. You may only specify one."
    )
  }
  if (missing_ratio) {
    if (!is.numeric(mean2) || any(mean2 <= 0)) {
      stop("Argument 'mean2' must be a positive numeric vector.")
    }
    ratio_keep <- NULL # For filtering step below
    ratio <- NULL # For checking length(ratio)
  }
  if (missing_mean2) {
    if (!is.numeric(ratio) || any(ratio <= 0)) {
      stop("Argument 'ratio' must be a positive numeric vector.")
    }
    mean2 <- as.numeric(tcrossprod(ratio, mean1))
    # For filtering step below
    # Numerical accuracy should be fine to 5 decimals...
    ratio_keep <- round(ratio, 5)
  }

  if (!is.numeric(dispersion) || any(dispersion <= 0)) {
    stop("Argument 'dispersion' must be a positive numeric vector.")
  }
  if (!is.numeric(nsims) || length(nsims) != 1L || nsims < 1L) {
    stop("Argument 'nsims' must be a positive scalar integer.")
  }
  if (length(return_type) != 1L) {
    stop("Argument 'return_type' must be one of 'list' or 'data.frame'.")
  }
  data.frame <- switch(
    return_type,
    "list" = FALSE,
    "data.frame" = TRUE,
    stop("Argument 'return_type' must be one of 'list' or 'data.frame'.")
  )
  if (
    !is.numeric(max_zeros) ||
      length(max_zeros) != 1L ||
      max_zeros < 0 ||
      max_zeros > 1
  ) {
    stop("Argument 'max_zeros' must be a scalar numeric from [0,1].")
  }

  needs_grid <- any(
    c(nsims, lengths(list(n, mean1, mean2, ratio, dispersion))) > 1L
  )

  #-----------------------------------------------------------------------------
  # Simulate data
  #-----------------------------------------------------------------------------
  res <- if (needs_grid) {
    grid_bnb(
      n1 = n,
      mean1 = mean1,
      mean2 = mean2,
      ratio = ratio,
      dispersion1 = dispersion,
      nsims = nsims,
      data.frame = data.frame,
      max_zeros = max_zeros,
      ratio_keep = ratio_keep
    )
  } else {
    sim_bnb_two_sample(
      n1 = n,
      mean1 = mean1,
      mean2 = mean2,
      dispersion1 = dispersion,
      nsims = nsims,
      data.frame = data.frame
    )
  }

  #-----------------------------------------------------------------------------
  # Return
  #-----------------------------------------------------------------------------
  res
}

grid_bnb <- function(
  n1,
  mean1,
  mean2,
  ratio,
  dispersion1,
  nsims,
  data.frame,
  max_zeros,
  ratio_keep
) {
  #-----------------------------------------------------------------------------
  # Unique combinations for simulating data.
  # Compute ratio in second step to filter unwanted rows for the case when
  # mean1 and ratio are provided in function.
  #-----------------------------------------------------------------------------
  grid_sim <- expand.grid(
    n1 = n1,
    mean1 = mean1,
    mean2 = mean2,
    dispersion1 = dispersion1,
    nsims = nsims,
    distribution = "Dependent two-sample BNB",
    stringsAsFactors = FALSE
  ) |>
    dplyr::mutate(ratio = round(mean2 / mean1, 5), .after = "mean2") |>
    {
      \(.) {
        if (!is.null(ratio_keep)) {
          dplyr::filter(.data = ., ratio %in% ratio_keep)
        } else {
          .
        }
      }
    }()

  #-----------------------------------------------------------------------------
  # Simulate data
  #-----------------------------------------------------------------------------
  # Simulate data
  res <- grid_sim |>
    dplyr::rowwise() |>
    dplyr::mutate(
      data = list(
        sim_bnb_two_sample(
          n1 = n1,
          mean1 = mean1,
          mean2 = mean2,
          dispersion1 = dispersion1,
          nsims = nsims,
          data.frame = data.frame
        )
      )
    ) |>
    dplyr::ungroup() |>
    dplyr::mutate(n2 = n1, .after = "n1")

  # Check if simulated data is all zeros or if the vast majority of data is all
  # zeros. This primarily occurs if you've selected a small dispersion (<0.1)
  # and have a small sample size.
  if (any_zeros(res[["data"]], max_zeros)) {
    res <- res |>
      dplyr::rowwise() |>
      dplyr::mutate(
        data = list(.data$data[not_zeros(.data$data, max_zeros)]),
        nsims = length(.data$data)
      ) |>
      dplyr::ungroup()
  }

  #-----------------------------------------------------------------------------
  # Prepare return
  #-----------------------------------------------------------------------------
  # Column labels
  vars <- c(
    "n1" = "n1",
    "n2" = "n2",
    "Mean1" = "mean1",
    "Mean2" = "mean2",
    "Ratio" = "ratio",
    "Dispersion" = "dispersion1",
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
  class(res) <- c("depower", "bnb", class(res))

  #-----------------------------------------------------------------------------
  # Return
  #-----------------------------------------------------------------------------
  res
}

sim_bnb_two_sample <- function(
  n1,
  mean1,
  mean2,
  dispersion1,
  nsims,
  data.frame
) {
  #-----------------------------------------------------------------------------
  # Simulate
  #-----------------------------------------------------------------------------
  if (nsims > 1L) {
    res <- lapply(seq_len(nsims), function(x) {
      g <- rgamma(n = n1, shape = dispersion1, scale = 1 / dispersion1)
      mean1 <- mean1 * g
      mean2 <- mean2 * g

      list(
        value1 = rpois(n = n1, lambda = mean1),
        value2 = rpois(n = n1, lambda = mean2)
      )
    })
    if (data.frame) {
      res <- lapply(res, list_to_df)
    }
  } else {
    # nsims == 1
    g <- rgamma(n = n1, shape = dispersion1, scale = 1 / dispersion1)
    mean1 <- mean1 * g
    mean2 <- mean2 * g

    res <- list(
      value1 = rpois(n = n1, lambda = mean1),
      value2 = rpois(n = n1, lambda = mean2)
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
