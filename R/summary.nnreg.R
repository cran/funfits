"summary.nnreg"<-
function(out, noprint = F, digits = 4)
{
	n <- out$n
	nfits <- out$nfits
	temp <- matrix(0, ncol = 6, nrow = nfits)
	for(j in 1:nfits) {
		fit <- out$model[[j]]
		temp[j, 1] <- fit$k
		ss <- sum((out$y - predict(fit, out$x))^2)
		temp[j, 2] <- fit$np
		temp[j, 3] <- n - fit$np
		temp[j, 4] <- signif(sqrt(ss/temp[j, 3]), digits)
		temp[j, 5] <- signif((ss/n)/(1 - temp[j, 2]/n)^2, digits)
		temp[j, 6] <- signif((ss/n)/(1 - (2 * temp[j, 2])/n)^2, digits)
	}
	dimnames(temp) <- list(format(1:nfits), c("# hidden units", "DF model", 
		"DF residuals", "Root MSE", "GCV", "GCV cost=2"))
	if(!noprint) {
		cat("Summary of outout from neural net fit", fill = T)
		cat("see the component summary in the output list for", fill = 
			T)
		cat("more details of the fitting process", fill = T)
		cat("call to nnreg :", fill = T)
		print(out$call)
		temp
	}
	else {
		return(temp)
	}
}
