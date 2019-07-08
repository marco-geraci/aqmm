\name{summary.aqmm}
\alias{summary.aqmm}

\title{Summary for an \code{aqmm} Object}

\description{
Summary method for class \code{aqmm}.}

\usage{
\method{summary}{aqmm}(x, R = 200, seed = NULL)
}

\arguments{
  \item{object}{
an object of \code{\link{class}} \code{aqmm}.
	}

  \item{R}{
number of bootstrap replications.
	}
  \item{seed}{
seed for the random number generator.
	}
}

\details{
\code{print.summary.aqmm} formats the coefficients, standard errors, etc. and additionally gives `significance stars'.}

\value{
an object of class \code{summary.aqmm}. The function \code{summary.aqmm} computes and returns a list of summary statistics of the fitted additive linear quantile mixed model given in \code{object}, that is

\item{tTable}{a matrix with estimates, standard errors, etc.}
\item{B}{the matrix of all bootstrapped parameters.}
}

\author{
Marco Geraci}

\seealso{
\code{\link{print.summary.aqmm}}
\code{\link{aqmm}}
}

\keyword{bootstrap}
\keyword{standard errors}
