"make.Amatrix.tps"<-
function(out, x0 = out$x, lambda, diagonal = F)
{
	if(missing(lambda)) {
		lambda <- out$lambda
	}
	xc <- out$transform$x.center
	xs <- out$transform$x.scale
	x <- scale(out$x, xc, xs)
	knots <- scale(out$knots, xc, xs)
	X <- cbind(make.tmatrix(x, out$m), qr.yq2(out$matrices$qr.T, make.rb(x, 
		knots, p = out$power, with.constant = out$with.constant)))
	if(!diagonal) {
		x0 <- scale(x0, xc, xs)
		temp <- (out$matrices$G) %*% diag(1/(1 + lambda * out$matrices$
			D))
		temp <- temp %*% t(out$matrices$G) %*% t(X)
		temp <- temp %*% diag(out$weights)	#
#
		temp <- cbind(make.tmatrix(x0, out$m), qr.yq2(out$matrices$qr.T,
			make.rb(x0, knots, p = out$power, with.constant = out$
			with.constant))) %*% temp	#
#
	}
	else {
		temp <- X %*% out$matrices$G %*% sqrt(diag(1/(1 + lambda * out$
			matrices$D)))
		temp <- c((temp^2) %*% rep(1, ncol(X))) * out$weights
	}
	return(temp)
}
