\name{aqmmControl}
\alias{aqmmControl}

\title{
Control parameters for aqmm estimation
}
\description{
A list of parameters for controlling the fitting process.
}

\usage{
aqmmControl(method = "Nelder-Mead", LL_tol = 1e-5, theta_tol = 1e-5,
	check_theta = FALSE, beta = 0.5, max_iter = 500, analytic = FALSE,
	REoptimizer = "optim", REcontrol = list(), verbose = FALSE)
}

\arguments{
	\item{method}{
		character vector that specifies the optimization algorithm in \code{\link{optim}} to maximize the log-likelihood. Only \code{"Nelder-Mead"} and \code{"BFGS"} allowed.
	}

	\item{LL_tol}{
		tolerance expressed as absolute change of the log-likelihood.
	}

	\item{theta_tol}{
		tolerance expressed as absolute change of the parameter \code{theta}.
	}

	\item{check_theta}{
		logical flag. If TRUE the algorithm performs an additional check on the change in the estimates.
	}

	\item{beta}{
		decreasing step factor for smoothing parameter \code{omega}. See details.
	}

	\item{max_iter}{
		maximum number of iterations.
	}
	\item{analytic}{
		logical flag. If \code{TRUE}, the gradient and Hessian for computing the modal random effects with \code{\link{nlm}} are calculated analytically. Otherwise, they are computed with numerical differentiation (default).
	}
	\item{REoptimizer}{
		optimizer for the modal random effects. The options are \code{"optim"} for BFGS in \code{\link{optim}} or \code{"nlm"} for Newton-type minimization in \code{\link{nlm}}.
	}
	\item{REcontrol}{
		a list of arguments to be passed to the optimizer for the modal random effects. See arguments in \code{\link{optim}} and \code{\link{nlm}}.
	}
	\item{verbose}{
		logical flag.
	}

}

\details{
The parameter \code{omega} controls the quadratic approximation of the absolute deviation function at the kink 0 (Chen, 2007). In \code{\link{aqmm}} this is provided via the argument \code{start} or using a naive estimate. The smaller \code{omega}, the smaller the approximation error. At each iteration, \code{omega} is decreased by a factor \code{beta}.
}

\value{
a list of control parameters.
}

\references{
Chen C. (2007). A finite smoothing algorithm for quantile regression. Journal of Computational and Graphical Statistics, 16(1), 136-164.
}

\author{
Marco Geraci}

\seealso{
\code{\link{aqmm}}
}
\keyword{fitting}
\keyword{control}
