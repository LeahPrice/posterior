% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/convergence.R
\name{ess_bulk}
\alias{ess_bulk}
\alias{ess_bulk.default}
\alias{ess_bulk.rvar}
\title{Bulk effective sample size (bulk-ESS)}
\usage{
ess_bulk(x, ...)

\method{ess_bulk}{default}(x, ...)

\method{ess_bulk}{rvar}(x, ...)
}
\arguments{
\item{x}{(multiple options) One of:
\itemize{
\item A matrix of draws for a single variable (iterations x chains). See
\code{\link[=extract_variable_matrix]{extract_variable_matrix()}}.
\item An \code{\link{rvar}}.
}}

\item{...}{Arguments passed to individual methods (if applicable).}
}
\value{
If the input is an array, returns a single numeric value. If any of the draws
is non-finite, that is, \code{NA}, \code{NaN}, \code{Inf}, or \code{-Inf}, the returned output
will be (numeric) \code{NA}. Also, if all draws within any of the chains of a
variable are the same (constant), the returned output will be (numeric) \code{NA}
as well. The reason for the latter is that, for constant draws, we cannot
distinguish between variables that are supposed to be constant (e.g., a
diagonal element of a correlation matrix is always 1) or variables that just
happened to be constant because of a failure of convergence or other problems
in the sampling process.

If the input is an \code{\link{rvar}}, returns an array of the same dimensions as the
\code{\link{rvar}}, where each element is equal to the value that would be returned by
passing the draws array for that element of the \code{\link{rvar}} to this function.
}
\description{
Compute a bulk effective sample size estimate (bulk-ESS) for a single
variable. Bulk-ESS is useful as a diagnostic for the sampling efficiency in
the bulk of the posterior. It is defined as the effective sample size for
rank normalized values using split chains. For the tail effective sample size
see \code{\link[=ess_tail]{ess_tail()}}. See Vehtari (2021) for an in-depth
comparison of different effective sample size estimators.
}
\examples{
mu <- extract_variable_matrix(example_draws(), "mu")
ess_bulk(mu)

d <- as_draws_rvars(example_draws("multi_normal"))
ess_bulk(d$Sigma)

}
\references{
Aki Vehtari, Andrew Gelman, Daniel Simpson, Bob Carpenter, and
Paul-Christian Bürkner (2021). Rank-normalization, folding, and
localization: An improved R-hat for assessing convergence of
MCMC (with discussion). \emph{Bayesian Data Analysis}. 16(2), 667-–718.
doi:10.1214/20-BA1221

Aki Vehtari (2021). Comparison of MCMC effective sample size estimators.
Retrieved from https://avehtari.github.io/rhat_ess/ess_comparison.html
}
\seealso{
Other diagnostics: 
\code{\link{ess_basic}()},
\code{\link{ess_quantile}()},
\code{\link{ess_sd}()},
\code{\link{ess_tail}()},
\code{\link{mcse_mean}()},
\code{\link{mcse_quantile}()},
\code{\link{mcse_sd}()},
\code{\link{rhat_basic}()},
\code{\link{rhat}()},
\code{\link{rstar}()}
}
\concept{diagnostics}
