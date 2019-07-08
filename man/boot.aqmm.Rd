\name{boot}
\alias{boot.aqmm}
\alias{boot}

\title{
Bootstrap functions for AQMM
}
\description{
This function is used to obtain a bootstrap sample of a fitted AQMM.
}
\usage{
boot.aqmm(x, R = 200, seed = round(runif(1, 1, 10000)))
}

\arguments{
  \item{x}{
an object of \code{\link{class}} "aqmm".
 }
  \item{R}{
number of bootstrap replications.
  }
  \item{seed}{
optional random number generator seed.
  }
}

\value{
An object of class \code{boot.aqmm} is a data frame with \code{R} rows and \code{npars} columns containing the bootstrap estimates of \code{theta_x}, \code{theta_z}, and \code{scale}.

The following attributes are available:

\item{tau}{index of the quantile(s).}
\item{estimated}{the estimated parameter as given by \code{object}.}
\item{std.err}{bootstrap standard errors.}
\item{R}{number of bootstrap replications.}
\item{seed}{the random number generator seed used to produce the bootstrap sample.}
\item{npars}{total numer of parameters.}
\item{rdf}{the number of residual degrees of freedom.}
}

\author{
Marco Geraci
}


\keyword{bootstrap}
\keyword{standard errors}
