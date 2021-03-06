"summary.krig"<-
function(x, digits = 4, ...)
{
	summary <- list(call = x$call, num.observation = length(x$residuals), 
		enp = x$trace, nt = x$nt, res.quantile = quantile(x$residuals, 
		seq(0, 1, 0.25)), shat.MLE = x$shat.MLE, shat.GCV = x$shat.GCV, 
		rhohat = x$rhohat, m = x$m, lambda = x$lambda, cost = x$cost, 
		gcvmin = min(x$gcv.grid[, 3]), rho = x$rho, sigma2 = x$sigma2)
	class(summary) <- "summary.krig"
	summary$covariance <- cor(x$fitted.values * sqrt(x$weights), (x$y) * 
		sqrt(x$weights))^2
	hold <- (sum((x$y - mean(x$y))^2) - sum(x$residuals^2))/(sum((x$y - 
		mean(x$y))^2))
	summary$adjr2 <- 1 - ((length(x$residuals) - 1)/(length(x$residuals) - 
		x$eff.df)) * (1 - hold)
	summary$digits <- digits
	summary
}
