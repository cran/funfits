"summary.tps"<-
function(x, digits = 4, ...)
{
## pdh - 7/16/96: added q2 and press
## pdh - 8/16/96: added rmse.press, pure.error, pure.df, fit.pure.error
## pdh - 11/25/97: added r.square as a pass through if available
##      and changed name of covariance to r.square
	summary <- list(call = x$call, num.observation = length(x$residuals), 
		enp = x$trace, nt = x$nt, res.quantile = quantile(x$residuals, 
		seq(0, 1, 0.25)), shat = x$shat, m = x$m, lambda = x$lambda, 
		form = x$form, power = x$power, cost = x$cost, gcvmin = min(x$
		gcv.grid[, 3]), press = x$press, r.square = x$r.square, q2 = x$
		q2, pure.error = x$shat.pure.error, pure.df = x$pure.df, method
		 = x$method, gcv.pure.error = x$gcv.pure.error, gcv.rmse = x$
		GCV, rmse.press = x$rmse.press)
	if(is.null(summary$method))
		summary$method <- "gcvmin"
	class(summary) <- "summary.tps"
	if(is.null(summary$r.square))
		summary$r.square <- cor(x$fitted.values * sqrt(x$weights), (x$y
			) * sqrt(x$weights))^2
	hold <- (sum((x$y - mean(x$y))^2) - sum(x$residuals^2))/(sum((x$y - 
		mean(x$y))^2))
	summary$adjr2 <- 1 - ((length(x$residuals) - 1)/(length(x$residuals) - 
		x$eff.df)) * (1 - hold)
	summary$digits <- digits
	summary
}
