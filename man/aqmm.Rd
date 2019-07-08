\name{aqmm}
\alias{aqmm}

\title{Fitting Additive Quantile Mixed Models}

\description{
\code{aqmm} is used to fit additive quantile mixed models based on the asymmetric Laplace distribution.}

\usage{
aqmm(fixed, random, group, knots = NULL, covariance = "pdDiag", tau = 0.5,
	data = sys.frame(sys.parent()), subset, weights, na.action = na.fail,
	control = list(), gamm = TRUE, start = NULL, gradHess = FALSE, fit = TRUE) 

}
\arguments{
  \item{fixed}{
an object of class \code{\link{formula}} for fixed effects: a symbolic description of the model to be fitted.
}
  \item{random}{
a one-sided formula of the form \code{~x1 + x2 + ... + xn} for random effects: a symbolic description of the model to be fitted.
}
  \item{group}{
grouping factor.
}
  \item{knots}{
this is an optional list containing user specified knot values to be used for basis construction. Different terms can use different numbers of knots, unless they share a covariate. This argument is passed to \code{\link{gam}} or \code{\link{gamm}} to obtain the spline design matrix.
}
  \item{covariance}{
variance-covariance matrix of the random effects. Default is \code{pdDiag} (see details).
}
  \item{tau}{
the quantile(s) to be estimated.
}
  \item{data}{
an optional data frame containing the variables named in
   \code{fixed}, \code{random} and \code{group}. By default the variables are taken from the environment from which \code{aqmm} is called.
}
  \item{subset}{
an optional vector specifying a subset of observations to be used in the fitting process.
}
  \item{weights}{
an optional vector specifying likelihood weights.
}
 \item{na.action}{a function that indicates what should happen when the
   data contain \code{NA}s.  The default action (\code{na.fail}) causes
   \code{lqmm} to print an error message and terminate if there are any
   incomplete observations.}

\item{control}{
  list of control parameters of the fitting process. See \code{\link{aqmmControl}}.
}
\item{gamm}{
logical flag. If \code{TRUE}, starting values for the fitting algorithm are
obtained from an additive mixed model (see \code{\link{gam}}).
}
  \item{start}{
  not yet implemented.
}
\item{gradHess}{
	logical flag. If \code{TRUE} gradient and Hessian of the log-likelihood are
	calculated numerically after convergence with code{numDeriv}. This can be time
	consuming.
}
 \item{fit}{
logical flag. If \code{FALSE} the function returns a list of arguments for fitting.
}
}


\details{
The function computes an estimate on the tau-th quantile function of the response, conditional on the covariates, as specified by the \code{fixed} argument, and on random effects, as specified by the \code{random} argument. The quantile predictor is assumed to be additive. Smooth terms can be included using the \code{\link{smooth.terms}} specification in \code{mgcv}.

Different types of positive--definite matrices for the random effects can be specified: \code{pdIdent} multiple of an identity; \code{pdCompSymm} compound symmetry structure (constant diagonal and constant off--diagonal elements); \code{pdDiag} diagonal; \code{pdSymm} general positive--definite matrix, with no additional structure.

Weights are given to clusters, therefore it is expected that these are constant within cluster. When the weights are specified in the main call, then the first value by \code{group} in the vector \code{weights} will be replicated for the same length of each group.
}

\value{
\code{aqmm} returns an object of \code{\link{class}} \code{aqmm}.
}

\references{
Geraci, M. (2018). Additive quantile regression for clustered data with an application to children's physical activity. Journal of the Royal Statistical Society - Series C. doi:10.1111/rssc.12333 (pre-print available at \url{https://arxiv.org/abs/1803.05403}).
}

\author{
Marco Geraci
}

\examples{
\dontrun{

library(aqmm)
set.seed(109)
M <- 100
n <- 10
z <- runif(M*n)
id <- rep(1:M, each = n)
y <- plogis(z, location = 0.5, scale = .1) + rep(rnorm(M, 0, .1), each = n) + rnorm(M*n, 0, .1)
dd <- data.frame(y = y, z = z, id = factor(id))
FIT <- aqmm(y ~ s(z), random = ~ 1, group = id, data = dd, tau = 0.9, gamm = TRUE)
FIT

plot(z, y)
points(z, predict(FIT, level = 0), col = 2, pch = 16)

}


}

\note{
Updates/FAQ/news are published here \url{http://marcogeraci.wordpress.com/}. New versions are usually published here \url{https://r-forge.r-project.org/R/?group_id=1396} before going on CRAN.
}

\seealso{
\code{\link{lqmm}}
}

\keyword{quantile regression}
