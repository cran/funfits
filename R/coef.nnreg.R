"coef.nnreg"<-
function(out, model = out$best.model)
{
	fit <- out$model[[model]]
	d <- fit$d
	if(is.null(fit$xm)) {
		fit$xm <- rep(0, fit$d)
		fit$xsd <- 1
	}
	if(is.null(fit$ym)) {
		fit$ym <- 0
		fit$ysd <- 1
	}
	k <- fit$k
	theta <- fit$theta
	beta <- theta[1:(k + 1)]
	mu <- theta[(1:k) + k + 1]
	gamma <- matrix(theta[(1:(d * k)) + 2 * k + 1], ncol = d, nrow = k, 
		byrow = T)
	dimnames(gamma) <- list(paste("gamma", 1:k, sep = ""), NULL)
	names(beta) <- paste("beta", 0:k, sep = "")
	names(mu) <- paste("mu", 1:k, sep = "")
	return(model, x.center = fit$xm, x.scale = fit$xsd, y.center = fit$ym, 
		y.scale = fit$ysd, d, k, beta, mu, gamma)
}
