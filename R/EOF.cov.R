"EOF.cov"<-
function(x1, x2, marginal = F, obj)
{
	symetric <- F | marginal
	if(missing(x2)) {
		symetric <- T
		x2 <- x1
	}
	M <- obj$M
	PHI1 <- matrix(0, ncol = M, nrow(x1))
	PHI2 <- matrix(0, ncol = M, nrow(x2))
	for(k in 1:M) {
		PHI1[, k] <- interp.surface(obj$fits[[k]], x1)
		if(!symetric) PHI2[, k] <- interp.surface(obj$fits[[k]], x2)	
	#		cat(k, " ")
	}
	if(marginal) {
#	print(PHI1)
		temp1 <- (t(sqrt(obj$delta[1:M]) * t(PHI1))^2) %*% rep(1, M)
		return(temp1)
	}
	if(!symetric) {
		temp2 <- PHI1 %*% (obj$delta[1:M] * t(PHI2))
		return(temp2)
	}
	else {
		temp2 <- PHI1 %*% ((obj$delta[1:M]) * t(PHI1))
		return(temp2)
	}
}
