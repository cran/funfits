"plot.nnreg"<-
function(out, model = out$best.model, main = NA, digits = 4, graphics.reset = T,
	...)
{
	old.par <- par("mfrow", "oma")
	if(graphics.reset) {
		on.exit(par(old.par))
		par(xpd = T)
	}
	if(model == out$best.model & (ncol(out$residuals) > 1)) {
		cat("Note: there is more than one model in the nnreg output object",
			fill = T)
	}
	set.panel(2, 2, T)
	temp <- summary(out, noprint = T)
	plot(out$fitted.values[, model], out$y, ylab = "Y", xlab = 
		"predicted values", bty = "n")
	abline(0, 1)
	hold <- par("usr")
	text(hold[1], hold[4], paste(" R**2 =", format(round(100 * cor(out$
		fitted.values[, model], out$y)^2, 2)), "%", sep = ""), cex = 
		0.80000000000000004, adj = 0)
	plot(out$fitted.values[, model], out$residuals[, model], ylab = 
		"residuals", xlab = "predicted values", bty = "n")
	yline(0)
	hold <- par("usr")
	text(hold[1], hold[4], paste(" RMSE =", format(signif(temp[model, 4], 
		digits))), cex = 0.80000000000000004, adj = 0)
	if(ncol(out$residuals) > 1) {
		matplot(temp[, 2], temp[, 5:6], ylab = "GCV (1) and GCV2 (2)", 
			xlab = "Number of Parameters", bty = "n", col = 1)
		xline(temp[model, 2])
		hold <- par("usr")
		text(temp[model, 2], hold[4], paste(" # par =", format(temp[
			model, 2]), "\n # units =", format(temp[model, 1])), 
			cex = 0.80000000000000004, adj = 0)
		title("GCV and GCV2", cex = 0.59999999999999998)
	}
	else {
		matplot(rep(temp[, 2], 2), temp[, 5:6], ylab = 
			"GCV (1) and GCV2 (2)", xlab = "Number of Parameters", 
			pty = 7, bty = "n", col = 1, type = "n")
		text(rep(temp[, 2], 2), temp[, 5:6], labels = c("1", "2"))
		xline(temp[model, 2])
		hold <- par("usr")
		text(temp[model, 2], hold[4], paste(" # par =", format(temp[
			model, 2]), "\n # units =", format(temp[model, 1])), 
			cex = 0.80000000000000004, adj = 0)
		title("GCV and GCV2", cex = 0.59999999999999998)
	}
	matplot(temp[, 2], temp[, 4], ylab = "RMSE", xlab = 
		"Number of Parameters", pch = "*", bty = "n")
	title("Root Mean Squared Error", cex = 0.59999999999999998)
	if(is.na(main))
		mtext(deparse(out$call), cex = 1.3, outer = T, line = -2)
	else mtext(main, cex = 1.3, outer = T, line = -2)
	invisible()
}
