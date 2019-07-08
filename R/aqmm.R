### Fit an additive quantile mixed model (aqmm)
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License or
#  any later version.
#
#  This program is distributed without any warranty,
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/
#

##########################################################################################
# aqmm functions

#' @export
aqmm <- function(fixed, random, group, knots = NULL, covariance = "pdDiag", tau = 0.5, data = sys.frame(sys.parent()), subset, weights, na.action = na.fail, control = list(), gamm = TRUE, start = NULL, gradHess = FALSE, fit = TRUE) {

Call <- match.call()
  
if(any(tau <= 0) | any(tau >= 1)) stop("Quantile index out of range")
if(length(tau) > 1) stop("One 'tau' at the time")

## check arguments
if(!is.data.frame(data)) stop("`data' must be a data frame")
if(!inherits(fixed, "formula") || length(fixed) != 3) {
	stop("\nFixed-effects model must be a formula of the form \"y ~ x\"")
}
if(!inherits(random, "formula") || length(random) != 2) {
	stop("\nRandom-effects model must be a formula of the form \" ~ x\"")
}

## use mgcv to interpret fixed formula
gp <- mgcv::interpret.gam(fixed)
sm <- gp$smooth.spec

## extract data frame with all necessary information
groupFormula <- asOneSidedFormula(Call[["group"]])
group <- groupFormula[[2]]
cov_name <- covariance

mfArgs <- list(formula = nlme::asOneFormula(random, gp$fake.formula, group), data = data, na.action = na.action)
if(!missing(subset)) {
	mfArgs[["subset"]] <- asOneSidedFormula(Call[["subset"]])[[2]]
}
#if(!missing(weights)) {
#	mfArgs[["weights"]] <- weights
#}
mfArgs$drop.unused.levels <- TRUE
dataMix <- do.call("model.frame", mfArgs)
origOrder <- row.names(dataMix)	# preserve the original order
for(i in names(contrasts)) contrasts(dataMix[[i]]) = contrasts[[i]]

## sort the model.frame by groups
grp <- model.frame(groupFormula, dataMix)

## ordering data by groups
ord <- order(unlist(grp, use.names = FALSE))
grp <- grp[ord,,drop = TRUE]
dataMix <- dataMix[ord, ,drop = FALSE]
revOrder <- match(origOrder, row.names(dataMix)) # putting in orig. order
ngroups <- length(unique(grp))

## use mgcv to get design matrices

if(!is.null(start)) gamm <- FALSE # coerce

if(gamm){
	reStruct <- list(group = nlme::pdMat(random, pdClass = cov_name))
	names(reStruct) <- as.character(group)
	gammFit <- mgcv::gamm(formula = fixed, random = reStruct, data = dataMix, knots = knots)
	s <- length(gammFit$gam$smooth)
	
	# response
	y <- gammFit$gam$model$y

	# fixed effects
	mmf <- gammFit$gam$model$X
	feTrms <- attr(mmf, "terms") <- colnames(mmf)

	# smooth terms
	Z.spline <- gammFit$gam$model[substr(names(gammFit$gam$model), 1, 2) == "Xr"]
	class(Z.spline) <- "list"
	names(Z.spline) <- sapply(gammFit$gam$smooth, function(x) x$label)
} else {
	gammFit <- 1
	class(gammFit) <- "try-error"
	gamFit <- mgcv::gam(formula = fixed, data = dataMix, knots = knots, fit = FALSE)
	allTrms <- gamFit$term.names
	X <- gamFit$X
	colnames(X) <- allTrms
	s <- length(gamFit$smooth)

	# response
	y <-  model.response(model.frame(gamFit$pterms, dataMix))

	# smooth terms
	Z.spline <- X[,grep("s\\(", allTrms, value = TRUE)]
	zcol <- sapply(gamFit$smooth, function(x) x$df)
	keep <- Z.spline[,cumsum(zcol),drop=FALSE]
	Z.spline <- Z.spline[,-cumsum(zcol)]
	zcol <- cbind(c(1, zcol[-s]), cumsum(zcol - 1))
	zcol <- split(zcol, row(zcol))
	Z.spline <- lapply(zcol, function(y, x) x[,seq(from = y[1], to = y[2])], x = Z.spline)
	names(Z.spline) <- sapply(gamFit$smooth, function(x) x$label)

	# fixed effects
	mmf <- X[,grep("s\\(", allTrms, value = TRUE, invert = TRUE),drop=FALSE]
	mmf <- cbind(mmf, keep)
	feTrms <- attr(mmf, "terms") <- colnames(mmf)
}

# random effects
mmr <- model.frame(random, dataMix)
mmr <- model.matrix(random, data = mmr)

# matrix Z as list
reTrms <- c(names(Z.spline), "group")
attr(reTrms, "group") <- colnames(mmr)
Z <- c(Z.spline, list(mmr)) # group is last
names(Z) <- reTrms

# likelihood weights
#if(!missing(weights)) weights <- model.weights(dataMix)[!duplicated(grp)]
#if(!missing(weights) && is.null(weights)) weights <- rep(1,ngroups)
if(missing(weights))  weights <- rep(1,ngroups)

## define dimensions
P <- length(feTrms)
Q <- ncol(mmr)
H <- sum(sapply(Z.spline, ncol))

dim_theta <- integer(3) # group is last
dim_theta[1] <- ncol(mmf)
dim_theta[2] <- s
dim_theta[3] <- theta.z.dim(type = cov_name, n = Q)
names(dim_theta) <- c("fixed","smooth","group")
dim_Z <- sapply(Z, ncol)
names(dim_Z) <- reTrms

## model info
model <- gp
model$dim_theta <- dim_theta
model$dim_Z <- dim_Z
model$feTrms <- feTrms
model$reTrms <- reTrms
model$fixed <- fixed

## Control
if(is.null(names(control))) control <- aqmmControl()
    else {
    control_default <- aqmmControl();
    control_names <- intersect(names(control), names(control_default));
    control_default[control_names] <- control[control_names];
    control <- control_default
    }

## Initialize

# if gamm is TRUE, initialize from gamm object
if(gamm & !inherits(gammFit, "try-error")){
	tmp <- getPars_gamm(gammFit, tau)
	theta_0 <- tmp$theta
	sigma_0 <- tmp$sigma
	RE_0 <- tmp$ranef
	omega_0 <- max(kronecker(abs(residuals(gammFit$lme, level = s + 1)), c(tau, 1 - tau), "/"))/2
	FIT_ARGS <- InitialPar <- list(theta = theta_0, sigma = sigma_0, ranef = RE_0, omega = omega_0)
}

# if gamm is FALSE and starting values not provided, initialize with naive estimates
if(is.null(start) & inherits(gammFit, "try-error")){
	QR <- qr(mmf)
	theta_0 <- c(qr.coef(QR, y), rep(0, sum(dim_theta[2:3])))
	sigma_0 <- invvarAL(mean(qr.resid(QR, y)*qr.resid(QR, y)), tau)
	RE_0 <- InitializeRE_aqmm(H = H, Q = Q, ngrp = ngroups)
	omega_0 <- sd(y)/2
	FIT_ARGS <- InitialPar <- list(theta = theta_0, sigma = sigma_0, ranef = RE_0, omega = omega_0)
}

# initialize with starting values if provided (gamm is coerced to FALSE whenever start is not null)
if(!is.null(start)){
	if(!setequal(names(start), c("theta", "sigma", "ranef", "omega"))) stop("Provide starting values as a named list with 'theta', 'sigma', 'ranef', 'omega'")
	if(length(start$theta) != sum(dim_theta)) stop(paste("Dimension of 'theta' must be", sum(dim_theta)))
	if(length(start$sigma) != 1) stop(paste("Dimension of 'sigma' must be", 1))
	if(length(start$ranef) != (Q*ngroups + H)) stop(paste("Dimension of 'ranef' must be", Q*ngroups + H))
	if(length(start$omega) != 1) stop(paste("Dimension of 'omega' must be", 1))
	if(is.null(attr(start$ranef, "hfirst"))){
		attr(start$ranef, "hfirst") <- rep(0, Q*ngroups + H)
	}
	if(is.null(attr(start$ranef, "hsec"))){
		attr(start$ranef, "hsec") <- 0
	}
	FIT_ARGS <- InitialPar <- start
}

## Complete FIT_ARGS list with all necessary arguments

FIT_ARGS <- c(FIT_ARGS, list(model = model, data = dataMix, y = y, X = mmf, Z = Z, weights = weights, cov_name = cov_name, group = grp, tau = tau, analytic = control$analytic, REoptimizer = control$REoptimizer, REcontrol = control$REcontrol))

if(!fit) return(FIT_ARGS)

## Estimation

iter <- 0
FIT_ARGS$ranef <- do.call(modalRe, args = FIT_ARGS[match(names(formals(modalRe)), names(FIT_ARGS))])
ll_0 <- do.call(loglik_aqmm, args = FIT_ARGS[match(names(formals(loglik_aqmm)), names(FIT_ARGS))])
LL_delta <- NA
while(iter < control$max_iter){

if(control$verbose){
	cat("iter:", iter, "\n")
	cat("loglikelihood:", ll_0, "\n")
	cat("delta:", LL_delta, "\n")
	cat("theta:", FIT_ARGS$theta, "\n")
	cat("sigma:", FIT_ARGS$sigma, "\n")
	cat("omega:", FIT_ARGS$omega, "\n")
}
# control = list(reltol = control$LP_tol_ll)
tmp <- do.call(optim, args = c(list(fn = loglik_aqmm, par = FIT_ARGS$theta, method = control$method), FIT_ARGS[-c(match(c("theta","analytic","REoptimizer","REcontrol"), names(FIT_ARGS)))]))

	if(control$check_theta){
		theta_delta <- abs((tmp$par - FIT_ARGS$theta)/FIT_ARGS$theta)
		FLAG <- max(theta_delta) < control$theta_tol
	} else {FLAG <- TRUE}

ll_1 <- tmp$value
LL_delta <- abs((ll_1 - ll_0)/(ll_0))
FIT_ARGS$theta <- tmp$par

	if(LL_delta < control$LL_tol & FLAG){
		cat("Convergence reached after ", iter + 1, "iteration(s)", "\n")
		break
	} else {
		ll_0 <- ll_1
		FIT_ARGS$omega <- FIT_ARGS$omega*control$beta
		iter <- iter + 1
	}
}

# update sigma, random effects and loglikelihood
FIT_ARGS$ranef <- do.call(modalRe, args = FIT_ARGS[match(names(formals(modalRe)), names(FIT_ARGS))])
if(control$verbose) cat(attr(FIT_ARGS$ranef, "msg"), "\n")
logLik <- do.call(loglik_aqmm, args = FIT_ARGS[match(names(formals(loglik_aqmm)), names(FIT_ARGS))])
FIT_ARGS$sigma <- attributes(logLik)$sigma

# Gradient and Hessian of loglike is time consuming
if(gradHess){
	attr(logLik, "grad") <- -do.call(numDeriv::grad, args = c(list(func = loglik_aqmm, x = FIT_ARGS$theta), FIT_ARGS[-c(match(c("theta","analytic","REoptimizer","REcontrol"), names(FIT_ARGS)))]))
	attr(logLik, "hessian") <- -do.call(numDeriv::hessian, args = c(list(func = loglik_aqmm, x = FIT_ARGS$theta), FIT_ARGS[-c(match(c("theta","analytic","REoptimizer","REcontrol"), names(FIT_ARGS)))]))
}

fit <- FIT_ARGS
fit$theta_x <- FIT_ARGS$theta[1:dim_theta[1]]
names(fit$theta_x) <- feTrms
fit$theta_z <- FIT_ARGS$theta[-c(1:dim_theta[1])]
names(fit$theta_z) <- c(paste0("sm", 1:dim_theta[2]), paste0("gr",1:dim_theta[3]))
fit$logLik <- as.numeric(-logLik)
fit$residuals <- attr(logLik, "resid")
fit$fitted <- attr(logLik, "fitted")
fit$call <- Call
fit$nobs <- length(y)
fit$dim_theta <- dim_theta
fit$edf <- sum(dim_theta)
fit$rdf <- fit$nobs - fit$edf
fit$df <- sum(dim_theta) + 1
fit$tau <- tau
fit$revOrder <- revOrder
fit$weights <- weights
fit$group <- grp
fit$ngroups <- ngroups
fit$InitialPar <- InitialPar
fit$control <- control
fit$opt <- list(iter = iter, delta = as.numeric(LL_delta))
fit$cov_name <- cov_name
fit$mfArgs <- mfArgs

class(fit) <- "aqmm"
fit
}

InitializeRE_aqmm <- function(H, Q, ngrp){

# Starting values of random effects. Random sample is not a good idea 

#val <- c(smooth = mapply("rnorm", n = H, SIMPLIFY = FALSE), group = list(matrix(rnorm(Q*ngrp, 0, 1), ngrp, Q)))
#val <- as.numeric(unlist(val))
val <- rep(0, H + Q*ngrp)

attr(val, "hfirst") <- rep(0, length(val))
attr(val, "hsec") <- 0

return(val)

}

#' @export
aqmmControl <- function(method = "Nelder-Mead", LL_tol = 1e-5, theta_tol = 1e-5, check_theta = FALSE, beta = 0.5, max_iter = 500, analytic = FALSE, REoptimizer = "optim", REcontrol = list(), verbose = FALSE){

if(!method %in% c("Nelder-Mead", "BFGS")) {method <- "Nelder-Mead"; cat("Switching to Nelder-Mead optimization \n")}

if(beta > 1 || beta < 0) stop("Beta must be a decreasing factor in (0,1)")
if(max_iter < 0) stop("Number of iterations cannot be negative")

list(method = method, LL_tol = LL_tol, theta_tol = theta_tol, check_theta = check_theta, beta = beta, max_iter = as.integer(max_iter), analytic = analytic, REoptimizer = REoptimizer, REcontrol = REcontrol, verbose = verbose)

}

getPars_gamm <- function(object, tau){

ObjClass <- names(object)[1]

if(!ObjClass %in% c("lme","mer")) stop("object must be a 'gamm' object")

if(ObjClass == "lme"){
	theta_x <- as.numeric(nlme::fixef(object$lme))
	sigma <- invvarAL(object$lme$sigma^2, tau)
	ranef <- as.numeric(unlist(nlme::ranef(object$lme)))
	attr(ranef, "hfirst") <- rep(0, length(ranef))
	attr(ranef, "hsec") <- 0
	theta_z <- coef(object$lme[[1]])
# theta_z from lme is scaled by sigma2
# reorder theta_z with group last (note that the order of the smooth terms needs to be inverted)
# divide smooth terms parameters by 2 to go from variance to standard deviation (i.e., from pdIdnot to pdIdent) so that all parameters for the random effects are on the same scale
	cov_name <- class(object$lme[[1]]$reStruct[[1]])
	sel <- cov_name %in% c("pdIdent","pdDiag","pdCompSymm","pdSymm")
	if(any(sel)) cov_name <- cov_name[sel] else stop("pdMat not recognized")
	Q <- length(attr(object$lme[[1]]$reStruct[[1]], "Dimnames")[[1]])
	s <- length(object$gam$smooth)
	m <- theta.z.dim(type = cov_name, n = Q)
	th1 <- theta_z[-c(1:m)]
	th1 <- th1[s:1]/2
	th2 <- theta_z[1:m]
	theta_z <- c(th1, th2)
}

if(ObjClass == "mer"){
	stop("Function not yet implemented for lme4 objects")
}

ans <- list(theta = c(theta_x, theta_z), sigma = sigma, ranef = ranef)
return(ans)
}

hfun_aqmm <- function(ranef, theta_x, Phi_inv, Sigma_inv, y, X, Z, weights, group, dim_Z, tau, omega, analytic){

eps <- .Machine$double.eps^(2/3)

id <- unique(group)
ni <- as.vector(table(group))
M <- length(id)
N <- length(y)
re <- length(dim_Z)
s <- re - 1
H <- sum(dim_Z[-re])
Q <- dim_Z['group']
P <- ncol(X)

B <- do.call(cbind, Z[-re])
Z <- Z[[re]]
v <- as.vector(ranef[c(1:H)])
u <- matrix(ranef[-c(1:H)], M, Q)

if(analytic){
	val <- C_hfunD_aqmm(y, X, B, Z, weights, theta_x, v, u, Phi_inv, Sigma_inv, M, N, ni, P, Q, H, s, tau, omega)

	A <- val[['weights']]
	Z <- lapply(split(Z, group), matrix, ncol = Q)
	sel <- rep(weights, table(group)) != 0
	G <- cbind(B[sel, , drop = FALSE], Matrix::bdiag(Z[weights != 0]))
	M <- M - sum(weights == 0)
	weights <- weights[weights != 0]
	Psi_inv <- Matrix::bdiag(diag(Phi_inv), kronecker(diag(weights, M, M), Sigma_inv))
	hessian <- 2*(Matrix::t(G) %*% Matrix::Diagonal(x = A[sel]) %*% G + Psi_inv)

	ans <- val[['val']]
	attr(ans, "gradient") <- val[['gradient']]
	attr(ans, "hessian") <- as.matrix(hessian)
} else {
	val <- C_hfun_aqmm(y, X, B, Z, weights, theta_x, v, u, Phi_inv, Sigma_inv, M, N, ni, P, Q, H, s, tau, omega)
	ans <- val[['val']]
}

return(ans)

}

#' @export
modalRe <- function(ranef, theta, sigma, model, y, X, Z, weights, cov_name, group, tau, omega, analytic, REoptimizer, REcontrol){

id <- unique(group)
M <- length(id)

if(length(weights) != M) stop("Length of \"weights\" does not match number of groups")

feTrms <- model$feTrms
reTrms <- model$reTrms

dim_theta <- model$dim_theta
dim_Z <- model$dim_Z
re <- length(dim_Z)
P <- dim_theta['fixed']
s <- dim_theta['smooth']
m <- dim_theta['group']
Q <- dim_Z['group']
N <- length(y)

theta_x <- theta[1:P]
theta_phi <- theta[(P + 1):(P + s)]
theta_csi <- theta[-c(1:(P + s))]

# smooth terms
Phi <- rep(NA, s)
for(i in 1:s){
	#Phi[i] <- as.matrix(pdIdnot(theta_phi[i], nam = "phi"))/sigma # !!!
	Phi[i] <- as.matrix(nlme::pdIdent(theta_phi[i], nam = "phi"))/sigma
}
Phi_inv <- 1/Phi
if(any(Phi_inv == Inf)) return(Inf)
Phi_inv <- rep(Phi_inv, dim_Z[-re])

# random effects
Sigma <- switch(cov_name,
		pdIdent = nlme::pdIdent(theta_csi, nam = attr(reTrms, "group")),
		pdDiag = nlme::pdDiag(theta_csi, nam = attr(reTrms, "group") ),
		pdSymm = nlme::pdSymm(theta_csi, nam = attr(reTrms, "group") ),
		pdCompSymm = nlme::pdCompSymm(theta_csi, nam = attr(reTrms, "group"))
)
Sigma <- as.matrix(Sigma)/sigma # !!!
Sigma_inv <- try(solve(Sigma), silent = TRUE)
if(inherits(Sigma_inv, "try-error")) Sigma_inv <- diag(100, Q, Q)

## nlm

if(REoptimizer == "nlm"){
	control_default <- list(hessian = FALSE, typsize = rep(1, length(ranef)), fscale = 1, print.level = 0, ndigit = 12, gradtol = 1e-3, stepmax = max(1000 * sqrt(sum((ranef/rep(1, length(ranef)))^2)), 1000), steptol = 1e-4, iterlim = 100, check.analyticals = FALSE)

	if(is.null(names(REcontrol))){
		REcontrol <- control_default} else {
		control_names <- intersect(names(REcontrol), names(control_default));
		control_default[control_names] <- REcontrol[control_names];
		REcontrol <- control_default
	}

	nlmArgs <- c(list(f = hfun_aqmm, p = ranef, theta_x = theta_x, Phi_inv = Phi_inv, Sigma_inv = Sigma_inv, y = y, X = X, Z = Z, weights = weights, group = group, dim_Z = dim_Z, tau = tau, omega = omega, analytic = analytic), REcontrol)

	fit <- do.call(nlm, nlmArgs)

	if(!REcontrol$hessian){
		fit$hessian <- attr(do.call(hfun_aqmm, list(ranef = fit$estimate, theta_x = theta_x, Phi_inv = Phi_inv, Sigma_inv = Sigma_inv, y = y, X = X, Z = Z, weights = weights, group = group, dim_Z = dim_Z, tau = tau, omega = omega, analytic = TRUE)), "hessian")
	}

	msg <- switch(as.character(fit$code),
		"1" = "relative gradient is close to zero, current iterate is probably solution",
		"2" = "successive iterates within tolerance, current iterate is probably solution",
		"3" = "last global step failed to locate a point lower than estimate. Either estimate is an approximate local minimum of the function or steptol is too small",
		"4" = "iteration limit exceeded",
		"5" = "maximum step size stepmax exceeded five consecutive times. Either the function is unbounded below, becomes asymptotic to a finite value from above in some direction or stepmax is too small"
	)
	
}

## optim

if(REoptimizer == "optim"){
	control_default <- list(method = "BFGS", hessian = FALSE)

	if(is.null(names(REcontrol))){
		REcontrol <- control_default} else {
		control_names <- intersect(names(REcontrol), names(control_default));
		control_default[control_names] <- REcontrol[control_names];
		REcontrol <- control_default
	}
	
	hfunD <- function(ranef, theta_x, Phi_inv, Sigma_inv, y, X, Z, weights, group, dim_Z, tau, omega, analytic) {attr(hfun_aqmm(ranef = ranef, theta_x = theta_x, Phi_inv = Phi_inv, Sigma_inv = Sigma_inv, y = y, X = X, Z = Z, weights = weights, group = group, dim_Z = dim_Z, tau = tau, omega = omega, analytic = TRUE), "gradient")}
	
	optimArgs <- c(list(par = ranef, fn = hfun_aqmm, theta_x = theta_x, Phi_inv = Phi_inv, Sigma_inv = Sigma_inv, y = y, X = X, Z = Z, weights = weights, group = group, dim_Z = dim_Z, tau = tau, omega = omega, analytic = FALSE, gr = hfunD), REcontrol)
	
	fit <- do.call(optim, optimArgs)
	names(fit)[names(fit) == "par"] <- "estimate"
	
	tmp <- do.call(hfun_aqmm, list(ranef = fit$estimate, theta_x = theta_x, Phi_inv = Phi_inv, Sigma_inv = Sigma_inv, y = y, X = X, Z = Z, weights = weights, group = group, dim_Z = dim_Z, tau = tau, omega = omega, analytic = TRUE))
	fit$gradient <- attr(tmp, "gradient")
	fit$hessian <- attr(tmp, "hessian")

	msg <- switch(as.character(fit$convergence),
		"0" = "successful completion",
		"1" = "iteration limit maxit had been reached",
		"10" = "degeneracy of the Nelder-Mead simplex",
		"51" = "warning from the 'L-BFGS-B' method",
		"52" = "error from the 'L-BFGS-B' method"
	)

}


val <- fit$estimate
hfirst <- fit$gradient
hsec <- as.numeric(determinant(fit$hessian/2, logarithm = TRUE)$modulus)

attr(val, "hfirst") <- hfirst
attr(val, "hsec") <- hsec
attr(val, "msg") <- paste("Algorithm for random effects:", msg)
return(val)

}

loglik_aqmm <- function(theta, sigma, ranef, model, data, y, X, Z, weights, cov_name, group, tau, omega) {
# Argument ranef is a (H + q) x 1 vector of random effects for each smooth
# and the grouping structure. Gradient and logDet(Hessian/2) of h function w.r.t. 'ranef' must be provided
# as attributes

eps <- .Machine$double.eps^(2/3)
id <- unique(group)
ni <- as.vector(table(group))
M <- length(id)

if(length(weights) != M) stop("Length of \"weights\" does not match number of groups")

feTrms <- model$feTrms
reTrms <- model$reTrms

dim_theta <- model$dim_theta

dim_Z <- model$dim_Z
re <- length(dim_Z)

P <- dim_theta['fixed']
s <- dim_theta['smooth']
m <- dim_theta['group']
Q <- dim_Z['group']
H <- sum(dim_Z[-re])
N <- length(y)

theta_x <- theta[1:P]
theta_phi <- theta[(P + 1):(P + s)]
theta_csi <- theta[-c(1:(P + s))]

# smooth terms
Phi <- rep(NA, s)
for(i in 1:s){
    #Phi[i] <- as.matrix(mgcv::pdIdnot(theta_phi[i], nam = "phi"))/sigma # !!!
	Phi[i] <- as.matrix(nlme::pdIdent(theta_phi[i], nam = "phi"))/sigma

}
Phi_inv <- 1/Phi
if(any(Phi_inv == Inf)) return(Inf)
Phi_inv <- rep(Phi_inv, dim_Z[-re])

# random effects
Sigma <- switch(cov_name,
                pdIdent = nlme::pdIdent(theta_csi, nam = attr(reTrms, "group")),
                pdDiag = nlme::pdDiag(theta_csi, nam = attr(reTrms, "group") ),
                pdSymm = nlme::pdSymm(theta_csi, nam = attr(reTrms, "group") ),
                pdCompSymm = nlme::pdCompSymm(theta_csi, nam = attr(reTrms, "group"))
)
Sigma <- as.matrix(Sigma)/sigma # !!!
Sigma_inv <- try(solve(Sigma), silent = TRUE)
if(inherits(Sigma_inv, "try-error")) return(Inf)

# log determinants
hsec <- attr(ranef, "hsec")
logDetH <- sum(hsec)
detH_detPsi <- sum(dim_Z[-re]*log(Phi)) + M*determinant(Sigma, logarithm = TRUE)$modulus + logDetH

# (negative) smoothed loss function
B <- do.call(cbind, Z[-re])
Z <- Z[[re]]
v <- as.vector(ranef[c(1:H)])
u <- matrix(ranef[-c(1:H)], M, Q)

val <- C_hfun_aqmm(y, X, B, Z, weights, theta_x, v, u, Phi_inv, Sigma_inv, M, N, ni, P, Q, H, s, tau, omega)
hsum <- val[['val']]

ans <- N * log(tau * (1- tau)/sigma) - 0.5*(detH_detPsi + hsum/sigma)
sigma <- hsum/(2*N)

attr(ans, "resid") <- as.numeric(val[['resid']])
attr(ans, "fitted") <- as.numeric(y - val[['resid']])
attr(ans, "sigma") <- sigma
return(-ans)
}

#' @export
VarCorr.aqmm <- function(x, sigma = NULL, rescaled = TRUE, ...){

reTrms <- x$model$reTrms
s <- x$dim_theta['smooth']
m <- x$dim_theta['group']

theta_phi <- x$theta_z[1:s]
theta_csi <- x$theta_z[-c(1:s)]
if(rescaled) {sigma <- if(is.null(sigma)) x$sigma else sigma} else {sigma <- 1}

Phi <- rep(NA, s)
for(i in 1:s){
    #Phi[i] <- as.matrix(mgcv::pdIdnot(theta_phi[i], nam = "phi")) * sigma
	Phi[i] <- as.matrix(nlme::pdIdent(theta_phi[i], nam = "phi")) * sigma
}
names(Phi) <- reTrms[1:s]

Sigma <- switch(x$cov_name,
		pdIdent = nlme::pdIdent(theta_csi, nam = attr(x$model$reTrms, "group")),
		pdDiag = nlme::pdDiag(theta_csi, nam = attr(x$model$reTrms, "group")),
		pdSymm = nlme::pdSymm(theta_csi, nam = attr(x$model$reTrms, "group")),
		pdCompSymm = nlme::pdCompSymm(theta_csi, nam = attr(x$model$reTrms, "group"))
)

return(list(group = as.matrix(Sigma) * sigma, smooth = Phi))

}

#' @export
print.aqmm <- function(x, digits = max(3, getOption("digits") - 3), ...){
tau <- x$tau
feTrms <- x$model$feTrms
reTrms <- x$model$reTrms

theta_x <- x$theta_x
sigma <- x$sigma
names(theta_x) <- feTrms

P <- x$dim_theta['fixed']
s <- x$dim_theta['smooth']
m <- x$dim_theta['group']

theta_phi <- x$theta[(P + 1):(P + s)]
theta_csi <- x$theta[-c(1:(P + s))]

Phi <- rep(NA, s)
for(i in 1:s){
    #Phi[i] <- as.matrix(mgcv::pdIdnot(theta_phi[i], nam = "phi")) * sigma
	Phi[i] <- as.matrix(nlme::pdIdent(theta_phi[i], nam = "phi")) * sigma
}
names(Phi) <- reTrms[1:s]

Sigma <- switch(x$cov_name,
			pdIdent = nlme::pdIdent(theta_csi, nam = attr(reTrms, "group")),
			pdDiag = nlme::pdDiag(theta_csi, nam = attr(reTrms, "group") ),
			pdSymm = nlme::pdSymm(theta_csi, nam = attr(reTrms, "group") ),
			pdCompSymm = nlme::pdCompSymm(theta_csi, nam = attr(reTrms, "group"))
)
Sigma <- as.matrix(Sigma) * sigma


cat("Call: ")
dput(x$call)
cat("\n")
cat(paste("Quantile", tau, "\n"))
cat("\n")
cat("Fixed effects:\n")
print.default(format(theta_x, digits = digits), print.gap = 2, 
	quote = FALSE)
cat("\n")
cat("Covariance matrix of the random effects:\n")
print.default(format(as.matrix(Sigma), digits = digits), quote = FALSE)
cat("\n")
cat("Variances of the (random) smooth terms:\n")
print.default(format(Phi, digits = digits), quote = FALSE)
cat("\n")
cat(paste("Residual scale parameter: ", format(sigma, 
	digits = digits)), "\n")
cat(paste("Log-likelihood:", format(x$logLik, digits = digits), 
	"\n"))
cat(paste("Tuning parameter:", format(x$omega, digits = digits), 
	"\n"))
	
cat(paste("\nNumber of observations:", x$nobs, "\n"))
cat(paste("Number of groups:", x$ngroups, "\n"))

invisible(x)
}

#' @export
predict.aqmm <- function(object, level = 0, newdata = NULL, ...){

# level = 0 includes fixed effects and smooths, but not random effects

#if(!is.null(newdata)){
#	if(level > 0) stop("Only level-0 predictions with newdata")
#	mt <- model.frame(delete.response(terms(object$model$fake.formula)), newdata)
	#G <- mgcv:::gamm.setup(object$model$fixed, attr(mt, "terms"), data = mt)
#	G <- getFromNamespace("gamm.setup", "mgcv")(object$model$fixed, attr(mt, "terms"), data = mt)
#	X <- G$X[,c(match("(Intercept)", colnames(G$X)), grep("s\\(", colnames(G$X)))]
#	Z <- G$data[substr(names(G$data), 1, 2) == "Xr"]
#	class(Z) <- "list"
#} else {
	X <- object$X
	Z <- object$Z
#}

id <- unique(object$group)
M <- length(id)
re <- length(Z)
s <- re - 1
theta_x <- object$theta_x
ranef <- object$ranef
Zu <- rep(0, nrow(X))
dim_Z <- object$model$dim_Z
H <- sum(dim_Z[-re])
Q <- dim_Z['group']
stop <- cumsum(dim_Z[-re])
start <- c(1, (stop + 1)[-s])

for(i in 1:s){
	Zu <- Zu + Z[[i]] %*% matrix(ranef[start[i]:stop[i]])
}
if(level > 0){
	u <- split(matrix(ranef[-c(1:H)], M, Q), id)
	Z <- split(Z[[re]], object$group)
	Z <- lapply(Z, function(x) as.matrix(x))
	Z <- lapply(Z, function(x) matrix(x, ncol = Q))
	Zu <- Zu + do.call(rbind, Map('%*%', Z, u))
}

val <- X %*% theta_x + Zu
if(is.null(newdata)){
	val <- val[object$revOrder]
}
return(val)

}

blb_part <- function(object, b = NULL, seed = NULL){

# b is the size (number of clusters) of each subset (except for last if remainder > 0). Subsets are disjoint (i.e., partition)
if(!is.null(seed)) set.seed(seed)
M <- object$ngroups
if(is.null(b)) {b <- floor(M/min(c(M, 5)))}
	else {if(b > M) stop("b is > M"); if(b == 1) stop("b must be > 1")}

if(b/M > 0.5) print(paste0("Size is more than half ", M, ". Length of partition is one"))

if(!is.factor(object$group)) stop("Grouping factor must be a factor. Re-run aqmm")
id <- levels(object$group)
flag <- M %% b != 0
nb <- if(flag) floor(M/b) else M/b
I <- rep(b, nb)
if(flag) I[nb] <- M - b*(nb - 1)

Sets <- list()
Sets[[1]] <- sample(id, I[1], replace = F)
if(nb > 1){
for(i in 2:nb){
remain <- !id %in% unlist(Sets[1:(i - 1)])
Sets[[i]] <- if(I[i] > 1) sample(id[remain], I[i], replace = F) else id[remain]
}}

return(list(partition = Sets, size = I))

}

#' @export
blb <- function(object, R = 50, seed = round(runif(1, 1, 10000)), partition = NULL, b = NULL, ...){

if(!is.null(seed)) set.seed(seed)

tau <- object$tau
nt <- length(tau)

group <- object$group
id <- levels(group)

est <- c(object$theta, object$sigma)
npars <- length(est)
names(est) <- c(names(object$theta_x), names(object$theta_z), "sigma")
data <- object$mfArgs$data

if(is.null(partition)){
	blbout <- blb_part(object, b = b, seed = seed)
	} else {blbout <- partition}
	
partition <- blbout$partition
size <- blbout$size
s <- length(partition)

if(length(id) != sum(size)) stop("Something went wrong with the partition size")
bootmat <- array(NA, dim = c(npars, R, s), dimnames = list(par = names(est), replicate = 1:R, partition = 1:s))
rows <- split(1:length(group), group)

for(j in 1:s){
	print(paste("bag", j, "of", s))
	obsS <- rmultinom(R, size = sum(size), prob = rep(1/size[j],size[j]))
	data.sub <- subset(data, group %in% partition[[j]])
	for(i in 1:R){
		if((i %% 10) == 0) print(paste(i, "in bag", j))
		w <- obsS[,i]
		fit <- try(update(object, data = data.sub, weights = w), silent = TRUE)
			if(!inherits(fit, "try-error")){
				bootmat[,i,j] <- c(fit$theta, fit$sigma)
			}
		}
	}

rowVars <- function(x) apply(x, 1, var, na.rm = TRUE)


class(bootmat) <- "boot.aqmm"
attr(bootmat, "tau") <- tau
attr(bootmat, "estimated") <- est
attr(bootmat, "std.err") <- sqrt(rowMeans(apply(bootmat, 3, rowVars)))
attr(bootmat, "R") <- R
attr(bootmat, "seed") <- seed
attr(bootmat, "npars") <- npars
attr(bootmat, "indices") <- blbout
attr(bootmat, "rdf") <- object$rdf

return(bootmat)

}

#' @export
boot.aqmm <- function(x, R = 200, seed = round(runif(1, 1, 10000))){

if (!is.null(seed)) set.seed(seed)
group <- factor(x$group)
id <- levels(group)
est <- c(x$theta, x$sigma)
npars <- length(est)

M <- length(id)
bootmat <- array(NA, dim = c(R, npars), dimnames = list(NULL, NULL))

for(r in 1:R){
	sel <- sample(id, M, replace = TRUE)
	newd <- NULL
	for(i in 1:M){
		tmp <- subset(x$data, group == sel[i])
		tmp$ID <- i
		newd <- rbind(newd, tmp)
	}
	newfit <- try(update(x, data = newd, group = ID), silent = TRUE)
	if(!inherits(newfit, "try-error")){
		bootmat[r,] <- c(newfit$theta, newfit$sigma)
	}
}

class(bootmat) <- "boot.aqmm"
attr(bootmat, "tau") <- x$tau
attr(bootmat, "estimated") <- est
attr(bootmat, "std.err") <- apply(bootmat, 2, sd, na.rm = TRUE)
attr(bootmat, "R") <- R
attr(bootmat, "seed") <- seed
attr(bootmat, "npars") <- npars
attr(bootmat, "rdf") <- x$rdf
return(bootmat)

}

#' @export
summary.aqmm <- function(object, R = 200, seed = NULL){

B <- boot.aqmm(object, R = R, seed = seed)
est <- attr(B, "estimated")
stds <- attr(B, "std.err")

tP <- 2 * pt(-abs(est/stds), R - 1)
tTable <- data.frame(est, stds, est/stds, tP)
rownames(tTable) <- c(names(object$theta_x), names(object$theta_z), "sigma")
colnames(tTable) <- c('Estimate', 'Std. Error', 't value', 'Pr(>|t|)')
object$tTable <- tTable
object$B <- B
class(object) <- "summary.aqmm"

object

}

#' @export
print.summary.aqmm <- function(x, digits = max(3, getOption("digits") - 3), ...){
tau <- x$tau
feTrms <- x$model$feTrms
reTrms <- x$model$reTrms

theta_x <- x$theta_x
sigma <- x$sigma
names(theta_x) <- feTrms

P <- x$dim_theta['fixed']
s <- x$dim_theta['smooth']
m <- x$dim_theta['group']

theta_phi <- x$theta[(P + 1):(P + s)]
theta_csi <- x$theta[-c(1:(P + s))]

Phi <- rep(NA, s)
for(i in 1:s){
    #Phi[i] <- as.matrix(mgcv::pdIdnot(theta_phi[i], nam = "phi")) * sigma
	Phi[i] <- as.matrix(nlme::pdIdent(theta_phi[i], nam = "phi")) * sigma
}
names(Phi) <- reTrms[1:s]

Sigma <- switch(x$cov_name,
			pdIdent = nlme::pdIdent(theta_csi, nam = attr(reTrms, "group")),
			pdDiag = nlme::pdDiag(theta_csi, nam = attr(reTrms, "group") ),
			pdSymm = nlme::pdSymm(theta_csi, nam = attr(reTrms, "group") ),
			pdCompSymm = nlme::pdCompSymm(theta_csi, nam = attr(reTrms, "group"))
)
Sigma <- as.matrix(Sigma) * sigma


cat("Call: ")
dput(x$call)
cat("\n")
cat(paste("Quantile", tau, "\n"))
cat("\n")
cat("Fixed effects:\n")
printCoefmat(x$tTable[1:x$dim_theta[1], ,drop = FALSE], signif.stars = TRUE, P.values = TRUE)
cat("\n")
cat("Covariance matrix of the random effects:\n")
print.default(format(as.matrix(Sigma), digits = digits), quote = FALSE)
cat("\n")
cat("Variances of the (random) smooth terms:\n")
print.default(format(Phi, digits = digits), quote = FALSE)
cat("\n")
cat(paste("Residual scale parameter: ", format(sigma, 
	digits = digits)), "\n")
cat(paste("Log-likelihood:", format(x$logLik, digits = digits), 
	"\n"))
cat(paste("Tuning parameter:", format(x$omega, digits = digits), 
	"\n"))
	
cat(paste("\nNumber of observations:", x$nobs, "\n"))
cat(paste("Number of groups:", x$ngroups, "\n"))

invisible(x)
}

##################################################
### Generics and methods (to be included in lqmm)
##################################################

