\name{predict.aqmm}
\alias{predict.aqmm}

\title{
Predictions from an \code{aqmm} Object
}

\description{
The predictions at level 0 correspond to predictions based on the terms specified in the \code{fixed} formula argument in \code{\link{aqmm}}. The predictions at level 1 are obtained by adding the estimated random effects to the predictions at level 0.}

\usage{
\method{predict}{aqmm}(object, level = 0, newdata = NULL, ...)
}

\arguments{
  \item{object}{
an \code{aqmm} object.
}
  \item{level}{
the level of grouping to be used in obtaining the predictions.
}
  \item{newdata}{
	not used.
}

  \item{\dots}{
not used.
}
}


\value{
a vector of predictions/fitted values.
}

\references{
Geraci M (2018). Additive quantile regression for clustered data with an application to children's physical activity. ArXiv e-prints, 1803.05403. Available at \url{https://arxiv.org/abs/1803.05403}.
}

\author{
Marco Geraci
}


\seealso{
\code{\link{aqmm}}
}

\keyword{prediction}
