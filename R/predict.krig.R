"predict.krig"<-
function(out, x, lambda = NA, model = NA)
{
	if(missing(x)) {
		x <- out$x
	}
	x <- as.matrix(x)	#
# scale the x values 
# using information from the output object
# scaling is (0,1) by default
#
	if(is.null(out$transform)) {
		xc <- rep(0, ncol(x))
		xs <- rep(1, ncol(x))
	}
	else {
		xc <- out$transform$x.center
		xs <- out$transform$x.scale
	}
	x <- scale(x, xc, xs)
	knots <- scale(out$knots, xc, xs)	#
# find out if fast multiplication routine for covarinace is available
#
# if cov.function allows an argument for C then use this typ of call 
# if the C argument is not present in teh function do it the long way 
# using matrix multiplcation explicitly
#
#
	if(!is.na(model)) {
		lambda <- model[1]
	}
	if(!is.na(lambda)) {
# use a different lambda so we need to get the new out$d and out$c
#coefficietns
		beta <- out$matrices$G %*% ((1/(1 + lambda * out$matrices$D)) * 
			out$matrices$u)
		nt <- out$nt
		np <- out$np
		out$d <- beta[1:nt]
		temp <- c(rep(0, nt), beta[(nt + 1):np])	#
#
# tranform the beta into the parameter associated with the covariance
# function  basis set into the c parameter vector.
#
		out$c <- c(qr.qy(out$matrices$qr.T, temp))
	}
##
###
# decide whether to use the fast multiplication routines for the
#covariance function
#	if(is.null(out$cov.function[["C"]])) {
	if(is.null(formals(out$cov.function)$C)) {
		c(make.tmatrix(x, out$m) %*% out$d + out$cov.function(x, out$
			knots) %*% out$c)
	}
	else {
		c(make.tmatrix(x, out$m) %*% out$d + out$cov.function(x, out$
			knots, C = out$c))
	}
}
