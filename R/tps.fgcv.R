"tps.fgcv"<-
function(lam, obj)
{
	lD <- obj$matrices$D * lam
	RSS <- obj$pure.ss + sum(((obj$matrices$u * lD)/(1 + lD))^2)
	MSE <- RSS/obj$N	#
	trA <- sum(1/(1 + lD))
	den <- (1 - (obj$cost * (trA - obj$nt) + obj$nt)/obj$N)	#
# If the denominator is negative then flag this as a bogus case
# by making the GCV fucntion "infinity"
#
	ifelse(den > 0, MSE/den^2, 1e+20)
}
