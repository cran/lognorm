% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/coefLognorm.R
\name{estimateParmsLognormFromSample}
\alias{estimateParmsLognormFromSample}
\alias{estimateStdErrParms}
\title{Estimate lognormal distribution parameters from a sample}
\usage{
estimateParmsLognormFromSample(x, na.rm = FALSE)

estimateStdErrParms(x, na.rm = FALSE)
}
\arguments{
\item{x}{numeric vector of sampled values}

\item{na.rm}{a logical value indicating whether 
NA values should be stripped before the computation proceeds.}
}
\value{
numeric vector with components \code{mu} and \code{sigma},
i.e., the center parameter (mean at log scale, log(median)) and 
the scale parameter (standard deviation at log scale)
}
\description{
Estimate lognormal distribution parameters from a sample
}
\details{
The expected value of a can be determined with
  higher accuracy the larger the sample. Here, the uncorrelated
  assumption is applied at the log scale and distribution parameters
  are returned with the same expected value as the sample, but with
  uncertainty (sigma) decreased by sqrt(nfin - 1). 
  
  Since with low relative error, the lognormal becomes very close
  to the normal distribution, the distribution of the mean can be
  well approximated by a normal with sd(mean(x)) ~ sd(x)/sqrt(n-1).
}
\section{Functions}{
\itemize{
\item \code{estimateParmsLognormFromSample}: Estimate lognormal distribution parameters from a sample

\item \code{estimateStdErrParms}: Estimate parameters of the lognormal distribution of the mean from an uncorrelated sample
}}

\examples{
.mu <- log(1)
.sigma <- log(2)
n = 200
x <- exp(rnorm(n, mean = .mu, sd = .sigma))
exp(pL <- estimateParmsLognormFromSample(x)) # median and multiplicative stddev
c(mean(x), meanx <- getLognormMoments(pL["mu"],pL["sigma"])[,"mean"])
c(sd(x), sdx <- sqrt(getLognormMoments(pL["mu"],pL["sigma"])[,"var"]))

# stddev decreases (each sample about 0.9) to about 0.07
# for the mean with n replicated samples
se <- estimateStdErrParms(x)
sqrt(getLognormMoments(se["mu"],se["sigma"])[,"var"])
sd(x)/sqrt(n-1) # well approximated by normal
# expected value stays the same
c(meanx, getLognormMoments(se["mu"],se["sigma"])[,"mean"])
}
