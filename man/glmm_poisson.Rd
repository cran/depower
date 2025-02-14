% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/glmm_poisson.r
\name{glmm_poisson}
\alias{glmm_poisson}
\title{GLMM for Poisson ratio of means}
\usage{
glmm_poisson(data, test = "wald", ci_level = NULL, ...)
}
\arguments{
\item{data}{(list)\cr
A list whose first element is the vector of Poisson values from sample
1 and the second element is the vector of Poisson values from sample 2.
Each vector must be sorted by the subject/item index and must be the
same sample size. \link[base]{NA}s are silently excluded. The default
output from \code{\link[=sim_bnb]{sim_bnb()}}.}

\item{test}{(String: \code{"wald"}; \verb{"c("wald", "lrt")})\cr
The statistical method used for the test results. \code{test = "wald"}
performs a Wald test and optionally the Wald confidence intervals.
\code{test = "lrt"} performs a likelihood ratio test and optionally
the profile likelihood confidence intervals (means and ratio). The
Wald interval is always used for the standard deviation of the item
(subject) random intercept.}

\item{ci_level}{(Scalar numeric: \code{NULL}; \verb{(0, 1)})\cr
If \code{NULL}, confidence intervals are set as \code{NA}. If in \verb{(0, 1)},
confidence intervals are calculated at the specified level.
Profile likelihood intervals are computationally intensive, so
intervals from \code{test = "lrt"} may be slow.}

\item{...}{Optional arguments passed to \code{\link[glmmTMB:glmmTMB]{glmmTMB::glmmTMB()}}.}
}
\value{
A list with the following elements:
\tabular{llll}{
Slot \tab Subslot \tab Name \tab Description \cr

1 \tab \tab \code{chisq} \tab \eqn{\chi^2} test statistic for the ratio of means. \cr
2 \tab \tab \code{df}    \tab Degrees of freedom. \cr
3 \tab \tab \code{p}     \tab p-value. \cr

4 \tab   \tab \code{ratio}    \tab Estimated ratio of means (sample 2 / sample 1). \cr
4 \tab 1 \tab \code{estimate} \tab Point estimate. \cr
4 \tab 2 \tab \code{lower}    \tab Confidence interval lower bound. \cr
4 \tab 3 \tab \code{upper}    \tab Confidence interval upper bound. \cr

5 \tab   \tab \code{mean1}    \tab Estimated mean of sample 1. \cr
5 \tab 1 \tab \code{estimate} \tab Point estimate. \cr
5 \tab 2 \tab \code{lower}    \tab Confidence interval lower bound. \cr
5 \tab 3 \tab \code{upper}    \tab Confidence interval upper bound. \cr

6 \tab   \tab \code{mean2}    \tab Estimated mean of sample 2. \cr
6 \tab 1 \tab \code{estimate} \tab Point estimate. \cr
6 \tab 2 \tab \code{lower}    \tab Confidence interval lower bound. \cr
6 \tab 3 \tab \code{upper}    \tab Confidence interval upper bound. \cr

7 \tab   \tab \code{item_sd}  \tab Estimated standard deviation of the item
(subject) random intercept. \cr
7 \tab 1 \tab \code{estimate} \tab Point estimate. \cr
7 \tab 2 \tab \code{lower}    \tab Confidence interval lower bound. \cr
7 \tab 3 \tab \code{upper}    \tab Confidence interval upper bound. \cr

8 \tab \tab \code{n1}           \tab Sample size of sample 1. \cr
9 \tab \tab \code{n2}           \tab Sample size of sample 2. \cr
10  \tab \tab \code{method}     \tab Method used for the results. \cr
11 \tab \tab \code{test}        \tab Type of hypothesis test. \cr
12 \tab \tab \code{alternative} \tab The alternative hypothesis. \cr
13 \tab \tab \code{ci_level}    \tab Confidence level of the interval. \cr
14 \tab \tab \code{hessian}     \tab Information about the Hessian matrix. \cr
15 \tab \tab \code{convergence} \tab Information about convergence.
}
}
\description{
Generalized linear mixed model for two dependent Poisson outcomes.
}
\details{
Uses \code{\link[glmmTMB:glmmTMB]{glmmTMB::glmmTMB()}} in the form

\if{html}{\out{<div class="sourceCode r">}}\preformatted{glmmTMB(
  formula = value ~ condition + (1 | item),
  data = data,
  family = stats::poisson
)
}\if{html}{\out{</div>}}

to model dependent Poisson outcomes \eqn{X_1 \sim \text{Poisson}(\mu)} and
\eqn{X_2 \sim \text{Poisson}(r \mu)} where \eqn{\mu} is the mean of sample 1
and \eqn{r} is the ratio of the means of sample 2 with respect to sample 1.

The hypotheses for the LRT and Wald test of \eqn{r} are

\deqn{
\begin{aligned}
H_{null} &: log(r) = 0 \\
H_{alt} &: log(r) \neq 0
\end{aligned}
}

where \eqn{r = \frac{\bar{X}_2}{\bar{X}_1}} is the population ratio of
arithmetic means for sample 2 with respect to sample 1 and
\eqn{log(r_{null}) = 0} assumes the population means are identical.

When simulating data from \code{\link[=sim_bnb]{sim_bnb()}}, the mean is a function of the
item (subject) random effect which in turn is a function of the dispersion
parameter. Thus, \code{glmm_poisson()} has biased mean estimates. The bias
increases as the dispersion parameter gets smaller and decreases as the
dispersion parameter gets larger. However, estimates of the ratio and
standard deviation of the random intercept tend to be accurate. In summary,
the Poisson mixed-effects model fit by \code{glmm_poisson()} is not recommended
for the BNB data simulated by \code{sim_bnb()}. Instead, \code{wald_test_bnb()} or
\code{lrt_bnb()} should typically be used instead.
}
\examples{
#----------------------------------------------------------------------------
# glmm_poisson() examples
#----------------------------------------------------------------------------
library(depower)

set.seed(1234)
d <- sim_bnb(
  n = 40,
  mean1 = 10,
  ratio = 1.2,
  dispersion = 2
)

lrt <- glmm_poisson(d, test = "lrt")
lrt

wald <- glmm_poisson(d, test = "wald", ci_level = 0.95)
wald

#----------------------------------------------------------------------------
# Compare results to manual calculation of chi-square statistic
#----------------------------------------------------------------------------
# Use the same data, but as a data frame instead of list
set.seed(1234)
d <- sim_bnb(
  n = 40,
  mean1 = 10,
  ratio = 1.2,
  dispersion = 2,
  return_type = "data.frame"
)

mod_alt <- glmmTMB::glmmTMB(
  formula = value ~ condition + (1 | item),
  data = d,
  family = stats::poisson,
)
mod_null <- glmmTMB::glmmTMB(
  formula = value ~ 1 + (1 | item),
  data = d,
  family = stats::poisson,
)

lrt_chisq <- as.numeric(-2 * (logLik(mod_null) - logLik(mod_alt)))
lrt_chisq
wald_chisq <- summary(mod_alt)$coefficients$cond["condition2", "z value"]^2
wald_chisq

anova(mod_null, mod_alt)

}
\references{
\insertRef{hilbe_2011}{depower}

\insertRef{hilbe_2014}{depower}
}
\seealso{
\code{\link[glmmTMB:glmmTMB]{glmmTMB::glmmTMB()}}
}
