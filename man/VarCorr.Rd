\name{VarCorr.aqmm}
\alias{VarCorr}
\alias{VarCorr.aqmm}

\title{
Extract Variance-Covariance Matrix
}
\description{
This function extracts the variance-covariance matrix of the random effects from a fitted \code{aqmm} object.
}
\usage{
\method{VarCorr}{aqmm}(x, sigma = NULL, rescaled = TRUE, ...)
}

\arguments{
  \item{x}{
an object of \code{\link{class}} "aqmm".
}

 \item{sigma}{
not used.
}

 \item{rescaled}{
logical flag. If \code{TRUE} the variance-covariance matrix of random effects is multiplied by \code{sigma} if provided or by the estimated scale parameter of the asymmetric Laplace. See details.
}

 \item{...}{
not used.
}

}

\details{
This function returns the variance-covariance matrix of the individual random effects as well the variances of the spline coefficients in the smooth terms.

The fitting algorithm rescales the variance-covariance matrices by \code{sigma} for estimation purposes. If estimates are required on the original scale, then arguments must be left at their default values.

The generic function \code{VarCorr} is imported from the \code{nlme} package (Pinheiro et al, 2014).
}


\value{
\code{VarCorr.aqmm} returns a list containing the following components:

\item{group}{the estimated variance--covariance matrix of the random effects.}
\item{smooth}{a vector of estimated variances of the spline coefficients in the smooth terms (variance components).}

}

\references{
Pinheiro J, Bates D, DebRoy S, Sarkar D and R Core Team (2014). nlme: Linear and Nonlinear Mixed Effects Models. R package version 3.1-117, \url{https://CRAN.R-project.org/package=nlme}.
}

\author{
Marco Geraci
}

\seealso{
\code{\link{aqmm}}
}

\keyword{covariance}
\keyword{coefficients}
