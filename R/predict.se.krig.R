"predict.se.krig"<-
function(out, x, cov.function, rho, sigma2, stationary = T)
{
	if(missing(x))
		x <- out$x
	x <- as.matrix(x)	
	########## if covariance function and parameters are missing
########## extract them from the krig object
	if(missing(cov.function)) {
		fun <- out$cov.function
	}
	else {
		fun <- cov.function
	}
	if(missing(sigma2)) {
		sigma2 <- out$sigma2
	}
#
# OK now fix the right value for sigma
#
	if(missing(rho)) {
		rho <- out$rho
	}
#
	if(!is.null(fun$marginal)) stationary <- F	#
	lambda <- sigma2/rho
	if(out$lambda != lambda) {
		warning("lambda value used is different from the one in the krig object"
			)
	}
	nx <- nrow(out$x)	
	# wght.vec are the linear combinations of the data that give the 
# correpsonding estimates of the function at the points x
	wght.vec <- t(make.Amatrix(out, x, lambda))	
	# Cy is the observed covariance matrix of the data vector
	Cy <- rho * fun(out$x, out$x) + sigma2 * diag(1/out$weights)
	temp2 <- c(t(wght.vec * (Cy %*% wght.vec)) %*% rep(1, nx))
	temp1 <- rho * c(t(wght.vec * fun(out$x, x)) %*% rep(1, nx))	#
#
#
	if(stationary) {
		x0 <- matrix(0, ncol = ncol(x), nrow = 1)
		return(sqrt(rho * fun(x0, x0) - 2 * temp1 + temp2))
	}
	else {
		if(is.null(fun$marginal)) {
#
#
# if covariance is not stationary then loop through each point to get
# the variance of field at that point. 
#
#
			temp <- rep(0, nrow(x))
			for(k in 1:nrow(x)) {
				x0 <- matrix(x[k,  ], nrow = 1)
				temp[k] <- rho * fun(x0, x0) - 2 * temp1[k] + 
				  temp2[k]
			}
		}
		else {
#
# marginal variances available by single call
#
			temp <- rho * fun(x, marginal = T) - 2 * temp1 + temp2
		}
		return(sqrt(temp))
	}
}
