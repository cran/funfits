"wt.huber"<-
function(u, c = 1.345)
{
	U <- abs(u)
	Ugtc <- (U > c)
	w <- u
	w[!Ugtc] <- 1
	w[Ugtc] <- c/U[Ugtc]
	w
}
