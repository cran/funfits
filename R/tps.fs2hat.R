"tps.fs2hat"<-
function(lam, obj)
{
	lD <- obj$matrices$D * lam
	RSS <- obj$pure.ss + sum(((obj$matrices$u * lD)/(1 + lD))^2)
	trA <- sum(1/(1 + lD))
	RSS/(obj$N - trA)	#
}
