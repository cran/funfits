"gcv.krig"<-
function(out, lambda = NA, cost = 1, nstep.cv = 80, verbose = F)
{
	nt <- out$nt
	np <- out$np
	N <- out$N
	D <- out$matrices$D
	u <- out$matrices$u
	pure.ss <- out$pure.ss
	if(is.na(lambda)) {
		l1 <- 1/D[np - nt - 1]
		tr <- np	########## find upper value of lambda
##########
		for(k in 1:8) {
			tr <- sum(1/(1 + l1 * D))
			if(tr < (nt + 0.050000000000000003))
				break
			l1 <- l1 * 2
		}
########## find lower lambda
##########
		l2 <- 1/D[1]
		for(k in 1:8) {
			tr <- sum(1/(1 + l2 * D))
			if(tr > (np * 0.94999999999999996))
				break
			l2 <- l2/2
		}
		lambda <- exp(seq(log(l2), log(l1),  , nstep.cv))
	}
	nl <- length(lambda)
	nd <- length(D)	#
	if(verbose) print(lambda)	#
## In S the fastest way to take a weighted sum of the columns of a matrix
##  is by  matrix multiplication
#
## A big matrix that is the product of the lambdas and D's
#
	big.lD <- matrix(D, nrow = nl, ncol = nd, byrow = T) * matrix(lambda, 
		ncol = nd, nrow = nl)	#
#
#
	if(verbose)
		print(pure.ss)
	RSS <- ((big.lD/(1 + big.lD))^2) %*% u^2
	if(verbose)
		print(RSS)
	RSS <- RSS + pure.ss
	MSE <- RSS/N	#
	trA <- (1/(1 + big.lD)) %*% rep(1, np)
	if(verbose)
		print(trA)	#	V <- MSE/(1 - (cost * (trA - nt) + nt)/N)^2
	denom <- (1 - (cost * (trA - nt) + nt)/N)
	V <- ifelse(denom > 0, MSE/denom^2, 1e+20)	#
## find global minimum of the GCV function
#
	gcv.grid <- data.frame(lambda, trA, V, sqrt(RSS/(N - trA)))
	names(gcv.grid) <- c("lambda", "trA", "GCV", "shat")
	il <- order(gcv.grid$GCV)[1]
	lambda.best <- gcv.grid$lambda[il]
	gcv.grid$GCV[gcv.grid$GCV == 1e+20] <- NA
	list(gcv.grid = gcv.grid, lambda.best = lambda.best)
}
