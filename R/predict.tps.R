"predict.tps"<-
function(out, x, y, lambda, df, omega, derivative = 0, model = NA)
{
# model is a generic argument that may be used to pass a different lambda
	if(!is.na(model)) lambda <- model
	if(out$tag != 1)
		stop("This is an old tps object please rerun\n\ntps to get right coefficients"
			)
	if(missing(x))
		x <- out$x
	x <- as.matrix(x)
	xc <- out$transform$x.center
	xs <- out$transform$x.scale
	n <- nrow(x)
	p <- ncol(x)
	x <- scale(x, xc, xs)
	knots <- scale(out$knots, xc, xs)
	nt <- out$nt
	np <- out$np	#
# these are the estimated coefficients to use from the tps object
#
	dtemp <- out$d
	ctemp <- out$c	#
# recompute the omega vector if somehting is different
#
	if(!missing(lambda) | !missing(df) | !missing(y)) {
		if(missing(lambda))
			lambda <- out$lambda
		if(!missing(df))
			lambda <- tps.df.to.lambda(df, out$matrices$D)
		if(!missing(y)) {
			u <- t(out$matrices$X %*% out$matrices$G) %*% (y * out$
				weights)
		}
		else {
			u <- out$matrices$u
		}
		omega <- out$matrices$G %*% ((1/(1 + lambda * out$matrices$D)) * 
			u)
		dtemp <- omega[1:nt]
		temp <- c(rep(0, nt), omega[(nt + 1):np])
		ctemp <- c(qr.qy(out$matrices$qr.T, temp))
	}
	if(!missing(omega)) {
		dtemp <- omega[1:nt]
		temp <- c(rep(0, nt), omega[(nt + 1):np])
		ctemp <- c(qr.qy(out$matrices$qr.T, temp))
	}
#
# at this point dtemp and ctemp are the right coefficients for the splines
#
	if(derivative == 0) {
		return(make.tmatrix(x, out$m) %*% dtemp + make.Kc(x, knots, 
			ctemp, p = out$power, with.constant = out$with.constant
			))
	}
	if(derivative == 1) {
		temp <- matrix(1/xs, ncol = p, nrow = n, byrow = T)
		return((make.DTd(x, dtemp, m = out$m) + make.DKc(x, knots, 
			ctemp, p = out$power, with.constant = out$with.constant
			)) * temp)
	}
}
