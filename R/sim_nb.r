#' Simulate data from a NB distribution
#'
#' Simulate data from two independent negative binomial (NB) distributions. For
#' paired data, see [depower::sim_bnb()].
#'
#' The negative binomial distribution has many parameterizations. In the
#' regression modeling context, it is common to specify the distribution in
#' terms of its mean and dispersion. We use the following probability
#' mass function:
#'
#' \deqn{
#' \begin{aligned}
#' P(X = x) &= \dbinom{x + \theta - 1}{x}
#'             \left( \frac{\theta}{\theta + \mu} \right)^{\theta}
#'             \left( \frac{\mu}{\mu + \theta} \right)^x \\
#'          &= \frac{\Gamma(x + \theta)}{x! \Gamma(\theta)}
#'             \left( \frac{\theta}{\theta + \mu} \right)^{\theta}
#'             \left( \frac{\mu}{\mu + \theta} \right)^{x} \\
#'          &= \frac{\Gamma(x + \theta)}{(\theta + \mu)^{\theta + x}}
#'             \frac{\theta^{\theta}}{\Gamma(\theta)} \frac{\mu^{x}}{x!}
#' \end{aligned}
#' }
#'
#' where \eqn{x \in \mathbb{N}^{\geq 0}}, \eqn{\theta \in \mathbb{R}^{> 0}}
#' is the dispersion parameter, and \eqn{\mu \in \mathbb{R}^{> 0}} is the mean.
#' This is analogous to the typical formulation where \eqn{X} is counting
#' \eqn{x} failures given \eqn{\theta} successes and
#' \eqn{p = \frac{\theta}{\theta + \mu}} is the probability of success on each
#' trial. It follows that \eqn{E(X) = \mu} and
#' \eqn{Var(X) = \mu + \frac{\mu^2}{\theta}}. The \eqn{\theta} parameter
#' describes the 'dispersion' among observations. Smaller values of \eqn{\theta}
#' lead to overdispersion and larger values of \eqn{\theta} decrease the
#' overdispersion, eventually converging to the Poisson distribution.
#'
#' Described above is the 'indirect quadratic parameterization' of the negative
#' binomial distribution, which is commonly found in the R ecosystem. However, it
#' is somewhat counterintuitive because the smaller \eqn{\theta} gets, the larger
#' the overdispersion. The 'direct quadratic parameterization' of the negative
#' binomial distribution may be found in some R packages and other languages
#' such as SAS and Stata. The direct parameterization is defined by substituting
#' \eqn{\alpha = \frac{1}{\theta}} (\eqn{\alpha > 0}) which results in
#' \eqn{Var(X) = \mu + \alpha\mu^2}. In this case, the larger \eqn{\alpha} gets
#' the larger the overdispersion, and the Poisson distribution is a special case
#' of the negative binomial distribution where \eqn{\alpha = 0}.
#'
#' A general class of negative binomial models may be defined with mean
#' \eqn{\mu} and variance \eqn{\mu + \alpha\mu^{p}}. The 'linear
#' parameterization' is then found by setting \eqn{p=1}, resulting in
#' \eqn{Var(X) = \mu + \alpha\mu}. It's common to label the linear
#' parameterization as 'NB1' and the direct quadratic parameterization as 'NB2'.
#'
#' See 'Details' in [depower::sim_bnb()] for additional information on the
#' gamma-Poisson mixture formulation of the negative binomial distribution.
#'
#' @references
#' \insertRef{yu_2017}{depower}
#'
#' \insertRef{rettiganti_2012}{depower}
#'
#' \insertRef{aban_2009}{depower}
#'
#' \insertRef{hilbe_2011}{depower}
#'
#' \insertRef{hilbe_2014}{depower}
#'
#' \insertRef{cameron_2013}{depower}
#'
#' @param n1 (integer: `[2, Inf)`)\cr
#'        The sample size(s) of group 1.
#' @param n2 (integer: `n1`; `[2, Inf)`)\cr
#'        The sample size(s) of group 2.
#' @param mean1 (numeric: `(0, Inf)`)\cr
#'        The mean(s) of group 1 \eqn{(\mu_1)}.
#' @param mean2,ratio (numeric: `(0, Inf)`)\cr
#'        Only specify one of these arguments.
#' - `mean2`: The mean(s) of group 2 \eqn{(\mu_2)}.
#' - `ratio`: The ratio(s) of means for group 2 with respect to group 1
#'        \eqn{\left( r = \frac{\mu_2}{\mu_1} \right)}.
#'
#' `mean2 = ratio * mean1`
#' @param dispersion1 (numeric: `(0, Inf)`)\cr
#'        The dispersion parameter(s) of group 1 \eqn{(\theta_1)}. See 'Details'
#'        and 'Examples'.
#' @param dispersion2 (numeric: `dispersion1`; `(0, Inf)`)\cr
#'        The dispersion parameter(s) of group 2 \eqn{(\theta_2)}. See 'Details'
#'        and 'Examples'.
#' @param nsims (Scalar integer: `1L`; `[1,Inf)`)\cr
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
#' @param ncores (Scalar integer: `1L`; `[1,Inf)`)\cr
#'        The number of cores (number of worker processes) to use. Do not set
#'        greater than the value returned by [parallel::detectCores()]. May be
#'        helpful when the number of parameter combinations is large and `nsims`
#'        is large.
#'
#' @return If `nsims = 1` and the number of unique parameter combinations is
#' one, the following objects are returned:
#' - If `return_type = "list"`, a list:
#' \tabular{lll}{
#'   Slot \tab Name \tab Description \cr
#'   1 \tab \tab Simulated counts from group 1. \cr
#'   2 \tab \tab Simulated counts from group 2.
#' }
#' - If `return_type = "data.frame"`, a data frame:
#' \tabular{lll}{
#'   Column \tab Name \tab Description \cr
#'   1 \tab `item` \tab Subject/item indicator. \cr
#'   2 \tab `condition` \tab Group/condition indicator. \cr
#'   3 \tab `value` \tab Simulated counts.
#' }
#'
#' If `nsims > 1` or the number of unique parameter combinations is greater than
#' one, each object described above is returned in a list-column named `data` in
#' a depower simulation data frame:
#' \tabular{lll}{
#'   Column \tab Name \tab Description \cr
#'   1 \tab `n1` \tab Sample size of group 1. \cr
#'   2 \tab `n2` \tab Sample size of group 2. \cr
#'   3 \tab `mean1` \tab Mean for group 1. \cr
#'   4 \tab `mean2` \tab Mean for group 2. \cr
#'   5 \tab `ratio` \tab Ratio of means (group 2 / group 1). \cr
#'   6 \tab `dispersion1` \tab Dispersion parameter for group 1. \cr
#'   7 \tab `dispersion2` \tab Dispersion parameter for group 2. \cr
#'   8 \tab `nsims` \tab Number of valid simulation iterations. \cr
#'   9 \tab `distribution` \tab Distribution sampled from. \cr
#'   10 \tab `data` \tab List-column of simulated data.
#' }
#'
#' @seealso [depower::sim_bnb()], [stats::rnbinom()]
#'
#' @examples
#' #----------------------------------------------------------------------------
#' # sim_nb() examples
#' #----------------------------------------------------------------------------
#' library(depower)
#'
#' # Independent two-sample NB data returned in a data frame
#' sim_nb(
#'   n1 = 10,
#'   mean1 = 5,
#'   ratio = 1.6,
#'   dispersion1 = 0.5,
#'   dispersion2 = 0.5,
#'   nsims = 1,
#'   return_type = "data.frame"
#' )
#'
#' # Independent two-sample NB data returned in a list
#' sim_nb(
#'   n1 = 10,
#'   mean1 = 5,
#'   ratio = 1.6,
#'   dispersion1 = 0.5,
#'   dispersion2 = 0.5,
#'   nsims = 1,
#'   return_type = "list"
#' )
#'
#' # Two simulations of independent two-sample data
#' # returned as a list of data frames
#' sim_nb(
#'   n1 = 10,
#'   mean1 = 5,
#'   ratio = 1.6,
#'   dispersion1 = 0.5,
#'   dispersion2 = 0.5,
#'   nsims = 2,
#'   return_type = "data.frame"
#' )
#'
#' # Two simulations of independent two-sample data
#' # returned as a list of lists
#' sim_nb(
#'   n1 = 10,
#'   mean1 = 5,
#'   ratio = 1.6,
#'   dispersion1 = 0.5,
#'   dispersion2 = 0.5,
#'   nsims = 2,
#'   return_type = "list"
#' )
#'
#' #----------------------------------------------------------------------------
#' # Visualization of the NB distribution as dispersion varies between groups.
#' #----------------------------------------------------------------------------
#' disp <- expand.grid(c(1, 10, 100), c(1, 10, 100))
#' set.seed(1234)
#' data <- mapply(
#'   FUN = function(disp1, disp2) {
#'     d <- sim_nb(
#'       n1 = 1000,
#'       mean1 = 10,
#'       ratio = 1.5,
#'       dispersion1 = disp1,
#'       dispersion2 = disp2,
#'       nsims = 1,
#'       return_type = "data.frame"
#'     )
#'     cbind(dispersion1 = disp1, dispersion2 = disp2, d)
#'   },
#'   disp1 = disp[[1]],
#'   disp2 = disp[[2]],
#'   SIMPLIFY = FALSE
#' )
#'
#' data <- do.call(what = "rbind", args = data)
#'
#' ggplot2::ggplot(
#'   data = data,
#'   mapping = ggplot2::aes(x = value, fill = condition)
#' ) +
#'   ggplot2::facet_grid(
#'     rows = ggplot2::vars(.data$dispersion2),
#'     cols = ggplot2::vars(.data$dispersion1),
#'     labeller = ggplot2::labeller(
#'       .rows = ggplot2::label_both,
#'       .cols = ggplot2::label_both
#'     )
#'   ) +
#'   ggplot2::geom_density(alpha = 0.3) +
#'   ggplot2::coord_cartesian(xlim = c(0, 50)) +
#'   ggplot2::labs(
#'     x = "Value",
#'     y = "Density",
#'     fill = "Condition",
#'     caption = "Mean1=10, Mean2=15, ratio=1.5"
#'   )
#'
#' @importFrom stats rnbinom
#'
#' @export
sim_nb <- function(
    n1,
    n2 = n1,
    mean1,
    mean2,
    ratio,
    dispersion1,
    dispersion2 = dispersion1,
    nsims = 1L,
    return_type = "list",
    max_zeros = 0.99,
    ncores = 1L
) {
  #-----------------------------------------------------------------------------
  # Check arguments
  #-----------------------------------------------------------------------------
  if(!is.numeric(n1) || any(n1 < 2L)) {
    stop("Argument 'n1' must be an integer vector from [2, Inf).")
  }
  if(!is.numeric(n2) || any(n2 < 2L)) {
    stop("Argument 'n2' must be an integer vector from [2, Inf).")
  }
  if(!is.numeric(mean1) || any(mean1 <= 0)) {
    stop("Argument 'mean1' must be a positive numeric vector.")
  }

  missing_mean2 <- missing(mean2)
  missing_ratio <- missing(ratio)
  if(missing_mean2 && missing_ratio) {
    stop("You must specify one of the arguments: 'mean2' or 'ratio'.")
  }
  if(!missing_mean2 && !missing_ratio) {
    stop("Arguments 'mean2' and 'ratio' were both specified. You may only specify one.")
  }
  if(missing_ratio) {
    if(!is.numeric(mean2) || any(mean2 <= 0)) {
      stop("Argument 'mean2' must be a positive numeric vector.")
    }
    ratio_keep <- NULL # For filtering step below
    ratio <- NULL # For checking length(ratio)
  }
  if(missing_mean2) {
    if(!is.numeric(ratio) || any(ratio <= 0)) {
      stop("Argument 'ratio' must be a positive numeric vector.")
    }
    mean2 <- as.numeric(tcrossprod(ratio, mean1))
    # For filtering step below
    # Numerical accuracy should be fine to 5 decimals...
    ratio_keep <- round(ratio, 5)
  }

  if(!is.numeric(dispersion1) || any(dispersion1 <= 0)) {
    stop("Argument 'dispersion1' must be a positive numeric vector.")
  }
  if(!is.numeric(dispersion2) || any(dispersion2 <= 0)) {
    stop("Argument 'dispersion2' must be a positive numeric vector.")
  }
  if(!is.numeric(nsims) || length(nsims) != 1L || nsims < 1L) {
    stop("Argument 'nsims' must be a positive scalar integer.")
  }
    if(length(return_type) != 1L) {
    stop("Argument 'return_type' must be one of 'list' or 'data.frame'.")
  }
  data.frame <- switch(return_type,
    "list" = FALSE,
    "data.frame" = TRUE,
    stop("Argument 'return_type' must be one of 'list' or 'data.frame'.")
  )
  if(!is.numeric(max_zeros) || length(max_zeros) != 1L || max_zeros < 0 || max_zeros > 1) {
    stop("Argument 'max_zeros' must be a scalar numeric from [0,1].")
  }
  if(!is.numeric(ncores) || length(ncores) != 1L || ncores < 1L) {
    stop("Argument 'ncores' must be a positive scalar integer.")
  }
  if(ncores > 1L) {
    if(isTRUE(ncores > parallel::detectCores())) {
      max <- parallel::detectCores()
      warning("Argument 'ncores' should not be greater than ", max, ".")
    }
    cluster <- multidplyr::new_cluster(ncores)
    multidplyr::cluster_library(cluster, 'depower')
  }

  needs_grid <- any(c(nsims, lengths(list(n1, n2,
                mean1, mean2, ratio, dispersion1, dispersion2))) > 1L)

  #-----------------------------------------------------------------------------
  # Simulate data
  #-----------------------------------------------------------------------------
  res <- if(needs_grid) {
    grid_nb(
      n1 = n1,
      n2 = n2,
      mean1 = mean1,
      mean2 = mean2,
      ratio = ratio,
      dispersion1 = dispersion1,
      dispersion2 = dispersion2,
      nsims = nsims,
      data.frame = data.frame,
      max_zeros = max_zeros,
      ncores = ncores,
      ratio_keep = ratio_keep
    )
  } else {
    sim_nb_two_sample(
      n1 = n1,
      n2 = n2,
      mean1 = mean1,
      mean2 = mean2,
      dispersion1 = dispersion1,
      dispersion2 = dispersion2,
      nsims = nsims,
      data.frame = data.frame
    )
  }

  #-----------------------------------------------------------------------------
  # Return
  #-----------------------------------------------------------------------------
  res
}

grid_nb <- function(
    n1,
    n2,
    mean1,
    mean2,
    ratio,
    dispersion1,
    dispersion2,
    nsims,
    data.frame,
    max_zeros,
    ncores,
    ratio_keep
) {
  #-----------------------------------------------------------------------------
  # Unique combinations for simulating data.
  # Compute ratio in second step to filter unwanted rows for the case when
  # mean1 and ratio are provided in function.
  #-----------------------------------------------------------------------------
  grid_sim <- expand.grid(
    n1 = n1,
    n2 = n2,
    mean1 = mean1,
    mean2 = mean2,
    dispersion1 = dispersion1,
    dispersion2 = dispersion2,
    nsims = nsims,
    distribution = "Independent two-sample NB",
    stringsAsFactors = FALSE
  ) |>
    dplyr::mutate(ratio = round(mean2 / mean1, 5), .after = "mean2") |>
    {\(.) if(!is.null(ratio_keep)) {dplyr::filter(.data = ., ratio %in% ratio_keep)} else {.}}()

  #-----------------------------------------------------------------------------
  # Simulate data
  #-----------------------------------------------------------------------------
  if(ncores > 1L) {
    cluster <- multidplyr::new_cluster(ncores)
    multidplyr::cluster_library(cluster, 'depower')
  }

  # Simulate data
  res <- grid_sim |>
    dplyr::rowwise() |>
    {\(.) if(ncores > 1L) {multidplyr::partition(data = ., cluster = cluster)} else {.}}() |>
    dplyr::mutate(
      data = list(
        sim_nb_two_sample(
          n1 = n1,
          n2 = n2,
          mean1 = mean1,
          mean2 = mean2,
          dispersion1 = dispersion1,
          dispersion2 = dispersion2,
          nsims = nsims,
          data.frame = data.frame
        )
      )
    ) |>
    {\(.) if(ncores > 1L) {dplyr::collect(x = .)} else {.}}() |>
    dplyr::ungroup()

  # Check if simulated data is all zeros or if the vast majority of data is all
  # zeros. This primarily occurs if you've selected a small dispersion (<0.1)
  # and have a small sample size.
  if(any_zeros(res[["data"]], max_zeros)) {
    res <- res |>
      dplyr::rowwise() |>
      dplyr::mutate(
        data = list(data[not_zeros(data, max_zeros)]),
        nsims = length(data)
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
    "Dispersion1" = "dispersion1",
    "Dispersion2" = "dispersion2",
    "N Simulations" = "nsims",
    "Distribution" = "distribution",
    "Data" = "data"
  )
  idx <- match(names(res), vars)
  if(anyNA(idx)) {stop("Unknown variable found while labeling data frame.")}
  for(i in seq_len(ncol(res))) {
    attr(res[[i]], "label") <- names(vars)[idx][i]
  }

  # Class
  class(res) <- c("depower", "nb", class(res))

  #-----------------------------------------------------------------------------
  # Return
  #-----------------------------------------------------------------------------
  res
}

sim_nb_two_sample <- function(
    n1,
    n2,
    mean1,
    mean2,
    dispersion1,
    dispersion2,
    nsims,
    data.frame
) {
  #-----------------------------------------------------------------------------
  # Simulate
  #-----------------------------------------------------------------------------
  if(nsims > 1L) {
    res <- lapply(seq_len(nsims), function(x) {
      list(
        value1 = rnbinom(n = n1, mu = mean1, size = dispersion1),
        value2 = rnbinom(n = n2, mu = mean2, size = dispersion2)
      )
    })
    if(data.frame) {
      res <- lapply(res, list_to_df)
    }
  } else { # nsims == 1
    res <- list(
      value1 = rnbinom(n = n1, mu = mean1, size = dispersion1),
      value2 = rnbinom(n = n2, mu = mean2, size = dispersion2)
    )
    if(data.frame) {
      res <- list_to_df(res)
    }
  }

  #-----------------------------------------------------------------------------
  # Return
  #-----------------------------------------------------------------------------
  res
}

utils::globalVariables(c("data"))
