"qr.q2y"<-
function(qr, y)
{
	c(qr.qy(qr, c(rep(0, qr$rank), y)))
}
