"summary.lle"<-
function(obj, digits = 5)
{
	cat("estimated global exponent", obj$glb, fill = T)
	cat("summary of QR estimate", fill = T)
	temp <- t(stats(obj$local.qr))
	temp <- signif(temp, digits)
	temp <- temp[, c(1, 2, 3, 5, 6, 7)]
	temp
}
