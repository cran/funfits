"print.tpsreg"<-
function(fit)
{
	cat(" Call: ", "\n", as.character(fit$call), fill = T)
	cat("Dimension of surface:", fit$parameter[2], fill = T)
	print.text(fit$summary)
	invisible()
}
