"plot.lle"<-
function(out)
{
	bplot(out$local.qr, xpos = log10(out$nprod), label = format(out$nprod), 
		xlab = "Number of Steps")
	yline(out$glb)
}
