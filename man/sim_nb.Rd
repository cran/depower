% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sim_nb.r
\name{sim_nb}
\alias{sim_nb}
\title{Simulate data from a NB distribution}
\usage{
sim_nb(
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
)
}
\arguments{
\item{n1}{(integer: \verb{[2, Inf)})\cr
The sample size(s) of group 1.}

\item{n2}{(integer: \code{n1}; \verb{[2, Inf)})\cr
The sample size(s) of group 2.}

\item{mean1}{(numeric: \verb{(0, Inf)})\cr
The mean(s) of group 1 \eqn{(\mu_1)}.}

\item{mean2, ratio}{(numeric: \verb{(0, Inf)})\cr
Only specify one of these arguments.
\itemize{
\item \code{mean2}: The mean(s) of group 2 \eqn{(\mu_2)}.
\item \code{ratio}: The ratio(s) of means for group 2 with respect to group 1
\eqn{\left( r = \frac{\mu_2}{\mu_1} \right)}.
}

\code{mean2 = ratio * mean1}}

\item{dispersion1}{(numeric: \verb{(0, Inf)})\cr
The dispersion parameter(s) of group 1 \eqn{(\theta_1)}. See 'Details'
and 'Examples'.}

\item{dispersion2}{(numeric: \code{dispersion1}; \verb{(0, Inf)})\cr
The dispersion parameter(s) of group 2 \eqn{(\theta_2)}. See 'Details'
and 'Examples'.}

\item{nsims}{(Scalar integer: \code{1L}; \verb{[1,Inf)})\cr
The expected number of simulated data sets. If \code{nsims > 1}, the data
is returned in a list-column of a depower simulation data frame.
\code{nsims} may be reduced depending on \code{max_zeros}.}

\item{return_type}{(string: \code{"list"}; \code{c("list", "data.frame")})\cr
The data structure of the simulated data. If \code{"list"} (default), a
list object is returned. If \code{"data.frame"} a data frame in tall format
is returned. The list object provides computational efficiency and the
data frame object is convenient for formulas. See 'Value'.}

\item{max_zeros}{(Scalar numeric: \code{0.99}; \verb{[0, 1]})\cr
The maximum proportion of zeros each group in a simulated dataset is
allowed to have. If the proportion of zeros is greater than this
value, the corresponding data is excluded from the set of simulations.
This is most likely to occur when the sample size is small and the
dispersion parameter is small.}

\item{ncores}{(Scalar integer: \code{1L}; \verb{[1,Inf)})\cr
The number of cores (number of worker processes) to use. Do not set
greater than the value returned by \code{\link[parallel:detectCores]{parallel::detectCores()}}. May be
helpful when the number of parameter combinations is large and \code{nsims}
is large.}
}
\value{
If \code{nsims = 1} and the number of unique parameter combinations is
one, the following objects are returned:
\itemize{
\item If \code{return_type = "list"}, a list:
\tabular{lll}{
Slot \tab Name \tab Description \cr
1 \tab \tab Simulated counts from group 1. \cr
2 \tab \tab Simulated counts from group 2.
}
\item If \code{return_type = "data.frame"}, a data frame:
\tabular{lll}{
Column \tab Name \tab Description \cr
1 \tab \code{item} \tab Subject/item indicator. \cr
2 \tab \code{condition} \tab Group/condition indicator. \cr
3 \tab \code{value} \tab Simulated counts.
}
}

If \code{nsims > 1} or the number of unique parameter combinations is greater than
one, each object described above is returned in a list-column named \code{data} in
a depower simulation data frame:
\tabular{lll}{
Column \tab Name \tab Description \cr
1 \tab \code{n1} \tab Sample size of group 1. \cr
2 \tab \code{n2} \tab Sample size of group 2. \cr
3 \tab \code{mean1} \tab Mean for group 1. \cr
4 \tab \code{mean2} \tab Mean for group 2. \cr
5 \tab \code{ratio} \tab Ratio of means (group 2 / group 1). \cr
6 \tab \code{dispersion1} \tab Dispersion parameter for group 1. \cr
7 \tab \code{dispersion2} \tab Dispersion parameter for group 2. \cr
8 \tab \code{nsims} \tab Number of valid simulation iterations. \cr
9 \tab \code{distribution} \tab Distribution sampled from. \cr
10 \tab \code{data} \tab List-column of simulated data.
}
}
\description{
Simulate data from two independent negative binomial (NB) distributions. For
paired data, see \code{\link[=sim_bnb]{sim_bnb()}}.
}
\details{
The negative binomial distribution has many parameterizations. In the
regression modeling context, it is common to specify the distribution in
terms of its mean and dispersion. We use the following probability
mass function:

\deqn{
\begin{aligned}
P(X = x) &= \dbinom{x + \theta - 1}{x}
            \left( \frac{\theta}{\theta + \mu} \right)^{\theta}
            \left( \frac{\mu}{\mu + \theta} \right)^x \\
         &= \frac{\Gamma(x + \theta)}{x! \Gamma(\theta)}
            \left( \frac{\theta}{\theta + \mu} \right)^{\theta}
            \left( \frac{\mu}{\mu + \theta} \right)^{x} \\
         &= \frac{\Gamma(x + \theta)}{(\theta + \mu)^{\theta + x}}
            \frac{\theta^{\theta}}{\Gamma(\theta)} \frac{\mu^{x}}{x!}
\end{aligned}
}

where \eqn{x \in \mathbb{N}^{\geq 0}}, \eqn{\theta \in \mathbb{R}^{> 0}}
is the dispersion parameter, and \eqn{\mu \in \mathbb{R}^{> 0}} is the mean.
This is analogous to the typical formulation where \eqn{X} is counting
\eqn{x} failures given \eqn{\theta} successes and
\eqn{p = \frac{\theta}{\theta + \mu}} is the probability of success on each
trial. It follows that \eqn{E(X) = \mu} and
\eqn{Var(X) = \mu + \frac{\mu^2}{\theta}}. The \eqn{\theta} parameter
describes the 'dispersion' among observations. Smaller values of \eqn{\theta}
lead to overdispersion and larger values of \eqn{\theta} decrease the
overdispersion, eventually converging to the Poisson distribution.

Described above is the 'indirect quadratic parameterization' of the negative
binomial distribution, which is commonly found in the R ecosystem. However, it
is somewhat counterintuitive because the smaller \eqn{\theta} gets, the larger
the overdispersion. The 'direct quadratic parameterization' of the negative
binomial distribution may be found in some R packages and other languages
such as SAS and Stata. The direct parameterization is defined by substituting
\eqn{\alpha = \frac{1}{\theta}} (\eqn{\alpha > 0}) which results in
\eqn{Var(X) = \mu + \alpha\mu^2}. In this case, the larger \eqn{\alpha} gets
the larger the overdispersion, and the Poisson distribution is a special case
of the negative binomial distribution where \eqn{\alpha = 0}.

A general class of negative binomial models may be defined with mean
\eqn{\mu} and variance \eqn{\mu + \alpha\mu^{p}}. The 'linear
parameterization' is then found by setting \eqn{p=1}, resulting in
\eqn{Var(X) = \mu + \alpha\mu}. It's common to label the linear
parameterization as 'NB1' and the direct quadratic parameterization as 'NB2'.

See 'Details' in \code{\link[=sim_bnb]{sim_bnb()}} for additional information on the
gamma-Poisson mixture formulation of the negative binomial distribution.
}
\examples{
#----------------------------------------------------------------------------
# sim_nb() examples
#----------------------------------------------------------------------------
library(depower)

# Independent two-sample NB data returned in a data frame
sim_nb(
  n1 = 10,
  mean1 = 5,
  ratio = 1.6,
  dispersion1 = 0.5,
  dispersion2 = 0.5,
  nsims = 1,
  return_type = "data.frame"
)

# Independent two-sample NB data returned in a list
sim_nb(
  n1 = 10,
  mean1 = 5,
  ratio = 1.6,
  dispersion1 = 0.5,
  dispersion2 = 0.5,
  nsims = 1,
  return_type = "list"
)

# Two simulations of independent two-sample data
# returned as a list of data frames
sim_nb(
  n1 = 10,
  mean1 = 5,
  ratio = 1.6,
  dispersion1 = 0.5,
  dispersion2 = 0.5,
  nsims = 2,
  return_type = "data.frame"
)

# Two simulations of independent two-sample data
# returned as a list of lists
sim_nb(
  n1 = 10,
  mean1 = 5,
  ratio = 1.6,
  dispersion1 = 0.5,
  dispersion2 = 0.5,
  nsims = 2,
  return_type = "list"
)

#----------------------------------------------------------------------------
# Visualization of the NB distribution as dispersion varies between groups.
#----------------------------------------------------------------------------
disp <- expand.grid(c(1, 10, 100), c(1, 10, 100))
set.seed(1234)
data <- mapply(
  FUN = function(disp1, disp2) {
    d <- sim_nb(
      n1 = 1000,
      mean1 = 10,
      ratio = 1.5,
      dispersion1 = disp1,
      dispersion2 = disp2,
      nsims = 1,
      return_type = "data.frame"
    )
    cbind(dispersion1 = disp1, dispersion2 = disp2, d)
  },
  disp1 = disp[[1]],
  disp2 = disp[[2]],
  SIMPLIFY = FALSE
)

data <- do.call(what = "rbind", args = data)

ggplot2::ggplot(
  data = data,
  mapping = ggplot2::aes(x = value, fill = condition)
) +
  ggplot2::facet_grid(
    rows = ggplot2::vars(.data$dispersion2),
    cols = ggplot2::vars(.data$dispersion1),
    labeller = ggplot2::labeller(
      .rows = ggplot2::label_both,
      .cols = ggplot2::label_both
    )
  ) +
  ggplot2::geom_density(alpha = 0.3) +
  ggplot2::coord_cartesian(xlim = c(0, 50)) +
  ggplot2::labs(
    x = "Value",
    y = "Density",
    fill = "Condition",
    caption = "Mean1=10, Mean2=15, ratio=1.5"
  )

}
\references{
\insertRef{yu_2017}{depower}

\insertRef{rettiganti_2012}{depower}

\insertRef{aban_2009}{depower}

\insertRef{hilbe_2011}{depower}

\insertRef{hilbe_2014}{depower}

\insertRef{cameron_2013}{depower}
}
\seealso{
\code{\link[=sim_bnb]{sim_bnb()}}, \code{\link[stats:NegBinomial]{stats::rnbinom()}}
}
