"plot.tpsreg"<-
function(out)
{
	old.par <- par()
        # error in R 0.61.3
        old.par$fin <- NULL
	if(length(out$gcv.grid[, 1]) > 1) {
		set.panel(2, 1)
		plot(out$gcv.grid[, 3], out$gcv.grid[, 2], type = "l", xlab = 
			"Effective number of parameters", ylab = "GCV function"
			)
		title(" Generalized Cross-validation function")
	}
	plot(out$fitted.value, out$residual, xlab = " Predicted values", ylab
		 = "Residuals")
	par(old.par)
	invisible()
}
