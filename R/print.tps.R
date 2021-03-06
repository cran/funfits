"print.tps"<-
function(x, ...)
{
	digits <- 4
	c1 <- "Number of Observations:"
	c2 <- length(x$residuals)
	c1 <- c(c1, "Degree of polynomial null space ( base model):")
	c2 <- c(c2, x$m - 1)
	c1 <- c(c1, "Number of parameters in the null space")
	c2 <- c(c2, x$nt)
	c1 <- c(c1, "Effective degrees of freedom:")
	c2 <- c(c2, format(round(x$eff.df, 1)))
	c1 <- c(c1, "Residual degrees of freedom:")
	c2 <- c(c2, format(round(length(x$residuals) - x$eff.df, 1)))
	c1 <- c(c1, "Root Mean Square Error:")
	c2 <- c(c2, format(round(x$shat, digits)))
	c1 <- c(c1, "Log10(lambda)")
	c2 <- c(c2, format(round(log10(x$lambda), 2)))
	sum <- cbind(c1, c2)
	dimnames(sum) <- list(rep("", dim(sum)[1]), rep("", dim(sum)[2]))
	cat("Call:\n")
	dput(x$call)
	print(sum, quote = F)
	invisible(x)
}
