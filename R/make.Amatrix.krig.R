"make.Amatrix.krig"<-
function(out, x0 = out$x, lambda)
{
	if(missing(lambda)) {
		lambda <- out$lambda
	}
	xc <- out$transform$x.center
	xs <- out$transform$x.scale
	x <- scale(out$x, xc, xs)
	knots <- scale(out$knots, xc, xs)
	x0 <- scale(x0, xc, xs)
	X <- cbind(make.tmatrix(x, out$m), qr.yq2(out$matrices$qr.T, out$
		cov.function(x, knots)))
	temp <- (out$matrices$G) %*% diag(1/(1 + lambda * out$matrices$D))
	temp <- temp %*% t(out$matrices$G) %*% t(X)
	temp <- temp %*% diag(out$weights)	#
#
	temp <- cbind(make.tmatrix(x0, out$m), qr.yq2(out$matrices$qr.T, out$
		cov.function(x0, knots))) %*% temp	#
#
	return(temp)
}
