"plot.tps"<-
function(out, main = NA, digits = 4, graphics.reset = T, ...)
{
## pdh 8/16/96 - added Q^2, pure error, changed some labels and
##   locations of text, add a gcvmin line if fit.pure.error=T
## DWN 9/26/96 changed arguments to fit in with new tps object
## pdh 10/29/96 - added drop=F to GCV matrix
## pdh 11/14/96 - x and y axes turned around on fit vs obs plot
## pdh 11/25/96 - uses r.squared instead of covariance if available
	old.par <- par("mfrow", "oma")
	if(graphics.reset) {
		on.exit(par(old.par))
		par(xpd = T)
	}
	set.panel(2, 2, T)
	temp <- summary(out)
	par1 <- par(pty = "s")
	lims <- range(out$fitted.values, out$y)
	plot(out$fitted.values, out$y, xlim = lims, ylim = lims, ylab = 
		"Observed Values", xlab = "Predicted Values", bty = "n", ...)
	abline(0, 1)
	hold <- par("usr")
	if(!is.null(temp$r.square))
		r.square <- temp$r.square
	else r.square <- temp$covariance
	text(hold[1], hold[4], paste(" R^2 = ", format(round(100 * r.square, 2)
		), "%", "\n", " Q^2 = ", format(round(100 * temp$q2, 2)), "%", 
		sep = ""), cex = 0.80000000000000004, adj = 0)
	par(par1)
	maxres <- max(abs(out$residuals))
	plot(out$fitted.values, out$residuals, ylim = c( - maxres, maxres), 
		ylab = "Residuals", xlab = "Predicted values", bty = "n", ...)
	yline(0)
	hold <- par("usr")
	if(!is.na(out$shat.pure.error))
		text(hold[1], hold[4], paste(" RMSE =", format(signif(out$shat, 
			digits)), "\n", "Pure Error =", format(signif(out$
			shat.pure.error, digits))), cex = 0.80000000000000004, 
			adj = 0)
	else text(hold[1], hold[4], paste(" RMSE =", format(signif(out$shat, 
			digits))), cex = 0.80000000000000004, adj = 0)
	if(nrow(out$gcv.grid) > 1) {
## trim off + infinity due to pole in the denominator of GCV function
##with cost
		ind <- out$gcv.grid[, 3] < 1e+19
		out$gcv.grid <- out$gcv.grid[ind,  ]
		plot(out$gcv.grid[, 2], (out$gcv.grid[, 3]), xlab = 
			"Effective number of parameters", ylab = 
			"Estimated (EASE) + sigma**2 ", bty = "n")
		xline(out$eff.df)
		hold <- par("usr")	##    text(out$eff.df, hold[4], 
		text(hold[1], hold[4], paste(" Eff. df. =", format(round(out$
			eff.df, 1)), "\n Res. df. =", format(round(temp$
			num.observation - temp$enp, 1))), cex = 
			0.80000000000000004, adj = 0)
		title("GCV", cex = 0.59999999999999998)
		plot(out$gcv.grid[, 1], out$gcv.grid[, 3], xlab = "Lambda", 
			ylab = "GCV", log = "x", bty = "n")
		temp <- out$lambda.est[!is.na(out$lambda.est[, "lambda"]),  , 
			drop = F]	## temp <- out$lambda.est
		hold <- par("usr")
		lam <- temp[, 1]
		names(lam) <- row.names(temp)	##		print(lam)
		xline(lam)
		points(lam, temp[, "GCV"], mark = 1, cex = 1.1000000000000001)
		title(paste("GCV", "\n", " Lambda =", format(round(out$lambda, 
			5))), cex = 0.59999999999999998)
	}
	if(is.na(main))
		mtext(deparse(out$call), cex = 1.3, outer = T, line = -2)
	else mtext(main, cex = 1.3, outer = T, line = -2)
	invisible()
}
