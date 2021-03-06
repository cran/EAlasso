% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PBsampler.R
\name{PB.CI}
\alias{PB.CI}
\title{Provide \code{(1-alpha)\%} confidence interval of each coefficients}
\usage{
PB.CI(object, alpha = 0.05, method = "debias", parallel = FALSE,
  ncores = 2L)
}
\arguments{
\item{object}{bootstrap samples of class \code{PB} from \code{\link{PBsampler}}}

\item{alpha}{significance level.}

\item{method}{bias-correction method. Either to be "none" or "debias".}

\item{parallel}{logical. If \code{TRUE}, use parallelization. Default is \code{FALSE}.}

\item{ncores}{integer. The number of cores to use for parallelization.}
}
\value{
\code{(1-alpha)\%} confidence interval of each coefficients
}
\description{
Using samples drawn by \code{\link{PBsampler}}, computes
\code{(1-alpha)\%} confidence interval of each coefficient.
}
\details{
If \code{method==none}, \code{\link{PB.CI}} simply compute
the two-soded \code{(1-alpha)} quantile of the sampled coefficients.
If \code{method==debias}, we use
debiased estimator to compute confidence interval.
}
\examples{
set.seed(1234)
n <- 40
p <- 50
Niter <-  10
Group <- rep(1:(p/10), each = 10)
Weights <- rep(1, p/10)
X <- matrix(rnorm(n*p), n)
object <- PBsampler(X = X, PE_1 = c(1,1,rep(0,p-2)), sig2_1 = 1, lbd_1 = .5,
niter = 100, type = "lasso")
parallel <- (.Platform$OS.type != "windows")
PB.CI(object = object, alpha = .05, method = "none")

}
\references{
Zhang, C., Zhang, S. (2014), "Confidence intervals for low dimensional
parameters in high dimensional linear models," Journal of the Royal
Statistical Society: Series B, 76, 217–242.

Dezeure, R., Buehlmann, P., Meier, L. and Meinshausen, N. (2015),
"High-Dimensional Inference: Confidence Intervals, p-values and R-Software hdi,"
Statistical Science, 30(4), 533-558
}
