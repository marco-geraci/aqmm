\name{blb}
\alias{blb}

\title{
Bag of little bootstraps for AQMM
}
\description{
This function is used to obtain a bootstrap sample of a fitted AQMM
}
\usage{
blb(object, R = 50, seed = round(runif(1, 1, 10000)), partition = NULL, b = NULL, ...)
}

\arguments{
  \item{object}{
an object of \code{\link{class}} \code{"aqmm"}.
 }
  \item{R}{
number of bootstrap replications.
  }
  \item{seed}{
optional random number generator seed.
  }
  \item{partition}{
pre-defined partition for 'bag of little bootstraps' (see details).
	}
  \item{b}{
size of each bag.
	}
  \item{...}{
not used.
	}
	
}

\value{
An object of class \code{boot.aqmm}. The function \code{blb} specifically returns an array of dimension \code{c(npars, R, s)} with the bootstrap estimates of \code{theta} and \code{sigma}, where \code{npars} is the number of parameters, \code{R} is the number of replications, and \code{s} is the number of bootstrap 'bags' (i.e., the number of subsets).

The following attributes are available:

\item{tau}{index of the quantile.}
\item{estimated}{the estimated parameter as given by \code{object}.}
\item{std.error}{the estimated boostrap standard errors.}
\item{R}{number of bootstrap replications.}
\item{seed}{the random number generator seed used to produce the bootstrap sample.}
\item{npars}{total numer of parameters.}
\item{indices}{the bootstrap sample of independent data units (clusters or groups).}
\item{rdf}{the number of residual degrees of freedom.}
}

\author{
Marco Geraci
}

\references{
Kleiner A, Talwalkar A, Sarkar P, and Jordan MI. (2014). A scalable bootstrap for massive data. Journal of the Royal Statistical Society B, 76(4), 795-816.
}

\keyword{bootstrap}
\keyword{standard errors}
