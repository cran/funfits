"plot.krig"<-
function(out, main = NA, digits = 4, graphics.reset = T, ...)
{
	old.par <- par("mfrow", "oma")
	if(graphics.reset) {
		on.exit(par(old.par))
		par(xpd = T)
	}
	set.panel(2, 2, T)
	temp <- summary(out)
	plot(out$fitted.values, out$y, ylab = "Y", xlab = " predicted values", 
		bty = "n", ...)
	abline(0, 1)
	hold <- par("usr")
	text(hold[1], hold[4], paste(" R**2 = ", format(round(100 * temp$
		covariance, 2)), "%", sep = ""), cex = 0.80000000000000004, adj
		 = 0)
	plot(out$fitted.values, out$residuals, ylab = "residuals", xlab = 
		" predicted values", bty = "n", ...)
	yline(0)
	hold <- par("usr")
	text(hold[1], hold[4], paste(" RMSE =", format(signif(sqrt(sum(out$
		residuals^2)/(temp$num.observation - temp$enp)), digits))), cex
		 = 0.80000000000000004, adj = 0)
	if(nrow(out$gcv.grid) > 1) {
# trim off + infinity due to pole in the denominator of GCV function
#with cost
		ind <- out$gcv.grid[, 3] < 1e+19
		out$gcv.grid <- out$gcv.grid[ind,  ]
		plot(out$gcv.grid[, 2], (out$gcv.grid[, 3]), xlab = 
			"Effective number of parameters", ylab = 
			" estimated (EASE) + sigma**2 ", bty = "n")
		xline(out$eff.df)
		hold <- par("usr")
		text(out$eff.df, hold[4], paste(" Eff. df. =", format(round(out$
			eff.df, 1)), "\n Res. df. =", format(round(temp$
			num.observation - temp$enp, 1))), cex = 
			0.80000000000000004, adj = 0)
		title("GCV", cex = 0.59999999999999998)
		plot(out$gcv.grid[, 1], out$gcv.grid[, 3], xlab = "Lambda", 
			ylab = "GCV", log = "x", bty = "n")
		xline(out$lambda)
		hold <- par("usr")
		text(out$lambda, hold[4], paste(" Lambda =", format(round(out$
			lambda, 2))), cex = 0.80000000000000004, adj = 0)
		title("GCV", cex = 0.59999999999999998)
	}
	if(is.na(main))
		mtext(deparse(out$call), cex = 1.3, outer = T, line = -2)
	else mtext(main, cex = 1.3, outer = T, line = -2)
}
