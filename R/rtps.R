"rtps"<-
function(..., max.iter = 20, acc = 10 * .Machine$single.eps^0.5, finish = F, 
	conv.method = "lambda")
{
## take the standard tps function but use iterative reweighted
## least squares
## first fit the regular tps to get initial weights
## then do huber until convergence
## if finish=T, then finish off with two rounds of bisquare
	fit1 <- tps(...)
	resid0 <- fit1$resid
	lambda0 <- fit1$lambda
	converged <- F
	status <- "converged"
	conv <- NULL
	if(max.iter > 0) {
		for(i in 1:max.iter) {
			if(i == 1) {
				scale <- median(abs(resid0))/
				  0.67449999999999999
				w <- w0 <- wt.huber(resid0/scale)
			}
			else {
				w0 <- w
				lambda0 <- lambda
			}
			fit1 <- tps(..., weights = w)
			resid <- fit1$resid
			scale <- median(abs(resid))/0.67449999999999999
			w <- wt.huber(resid/scale)
			lambda <- fit1$lambda
			if(conv.method == "lambda")
				convi <- abs(lambda - lambda0)
			if(conv.method == "weights")
				convi <- sum((w - w0)^2)
			conv <- c(conv, convi)
			converged <- convi <= acc
			if(converged)
				break
		}
		if(!converged)
			warning(status <- paste("failed to converge in", 
				max.iter, "steps"))
	}
	if(finish) {
## now finish off with two iterations of the bisquare
		for(i in 1:2) {
			scale <- median(abs(resid))/0.67449999999999999
			wts <- wt.bisquare(resid/scale)
			fit1 <- tps(..., weights = wts)
			resid <- fit1$resid
		}
	}
	class(fit1) <- c("rtps", "tps", "funfits")
	fit1$call <- match.call()
	invisible(fit1)
}
