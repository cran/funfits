"predict.netfit"<-
function(fit, x, derivative = 0, type = "full")
{
	nx <- nrow(x)
	d <- fit$d
	if(is.null(fit$xm)) {
		fit$xm <- rep(0, fit$d)
		fit$xsd <- 1
	}
	if(is.null(fit$ym)) {
		fit$ym <- 0
		fit$ysd <- 1
	}
	if(d != ncol(x)) stop(" columns of X not equal to d!")	
	# standardize the X's 
	u <- (x - matrix(fit$xm, ncol = d, nrow = nx, byrow = T))
	u <- u/matrix(fit$xsd, ncol = d, nrow = nx, byrow = T)
	k <- fit$k
	theta <- fit$theta
	beta <- theta[1:(k + 1)]
	mu <- theta[(1:k) + k + 1]
	gamma <- matrix(theta[(1:(d * k)) + 2 * k + 1], ncol = d, nrow = k, 
		byrow = T)
	pu <- cbind(rep(1, nx), u) %*% t(cbind(mu, gamma))
	if(derivative == 0) {
# find predicted values and transform to original scale
		if(type == "full") {
			yhat <- fit$ysd * (squasher.nnreg(pu) %*% beta[2:(k + 1
				)] + beta[1]) + fit$ym
		}
		if(type == "terms") {
			temp <- dim(pu)
			yhat <- list(u = pu, yhat = fit$ysd * (squasher.nnreg(
				pu)) * matrix(beta[2:(k + 1)], nrow = temp[1], 
				ncol = temp[2], byrow = T), constant = fit$ysd * 
				beta[1] + fit$ym)
		}
		return(yhat)
	}
	else {
		if(type == "terms")
			stop("derviative not available for individual\nhidden units"
				)
		jac <- fit$ysd * (d.squasher.nnreg(pu) %*% (gamma/matrix(fit$
			xsd, ncol = d, nrow = k, byrow = T) * matrix(beta[2:(k + 
			1)], ncol = d, nrow = k)))
		return(jac)
	}
}
