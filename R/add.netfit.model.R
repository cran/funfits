"add.netfit.model"<-
function(fita, fitb)
{
	if(fita$d != fitb$d)
		stop(" mismatch in dimensions of two models")
	d <- fita$d
	ka <- fita$k
	kb <- fitb$k
	temp.a <- summary.netfit(fita, F)
	temp.b <- summary.netfit(fitb, F)
	temp <- list()
	class(temp) <- "netfit"
	temp$np <- fita$np + fitb$np - 1
	temp$d <- d
	temp$k <- ka + kb
	theta <- c(temp.a$beta[1] + temp.b$beta[1], temp.a$beta[2:(ka + 1)], 
		temp.b$beta[2:(kb + 1)], temp.a$mu, temp.b$mu, c(t(rbind(temp.a$
		gamma, temp.b$gamma))))
	temp$theta <- theta
	temp$xm <- rep(0, d)
	temp$xsd <- rep(1, d)
	temp$ym <- 0
	temp$ysd <- 1
	temp
}
