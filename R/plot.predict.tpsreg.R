"plot.predict.tpsreg"<-
function(x, model = 1, digits = 2)
{
	lims <- range(c(x$y, x$fit))
	n <- length(x$y)
	plot(x$y, x$fit, ylab = "Predicted", xlab = "Observed", xlim = lims, 
		ylim = lims, type = "n")
	abline(0, 1)
	points(x$y, x$fit, mark = 4)
	rmse <- sqrt(sum(x$res^2)/(n - x$eff.df))
	title(paste("TPSREG Results", "; Obs. =", length(x$y), "; R.sq =", 
		round(100 * cor(x$y, x$fit)^2, 1), "%", "\n", "DFr =", round(x$
		eff.df, 1), "; DFe =", round(length(x$y) - x$eff.df, 1), 
		"; RMSE =", signif(rmse, digits)))
}
