"gcv.tps"<-
function(out, lambda.grid = NA, cost = 1, nstep.cv = 80, rmse = NA, verbose = F,
	tol = 1.0000000000000001e-05)
{
	nt <- out$nt
	np <- out$np
	N <- length(out$y)
	D <- out$matrices$D
	u <- out$matrices$u
	shat.pure.error <- out$shat.pure.error
	pure.ss <- out$pure.ss	#
# create a reasonable grid for the GCV search if not supplied
#
	if(is.na(lambda.grid[1])) {
		l1 <- 1/D[np - nt - 1]
		tr <- np	#
########## find upper value of lambda
#
		for(k in 1:8) {
			tr <- sum(1/(1 + l1 * D))
			if(tr < (nt + 0.050000000000000003))
				break
			l1 <- l1 * 2
		}
#
########## find lower lambda
#
		l2 <- 1/D[1]
		for(k in 1:8) {
			tr <- sum(1/(1 + l2 * D))
			if((tr > (np * 0.94999999999999996)) | ((1 - (cost * (
				tr - nt) + nt)/N) <= 0))
				break
			l2 <- l2/2
		}
		lambda.grid <- exp(seq(log(l2), log(l1),  , nstep.cv))
	}
#
# done with finding a good default range for lambda
#
	nl <- length(lambda.grid)
	nd <- length(D)	#
#
## In S the fastest way to take a weighted sum of the columns of a matrix
##  is by  matrix multiplication
#
## Now make a big matrix that is the product of the lambdas and D's
#
	big.lD <- matrix(D, nrow = nl, ncol = nd, byrow = T) * matrix(
		lambda.grid, ncol = nd, nrow = nl)	#
#
#
	RSS <- pure.ss + ((big.lD/(1 + big.lD))^2) %*% u^2
	MSE <- RSS/N	#
	trA <- (1/(1 + big.lD)) %*% rep(1, np)
	den <- (1 - (cost * (trA - nt) + nt)/N)	#
# If the denominator is negative then flag this as a bogus case
# by making the GCV function "infinity": 10^20
#
	V <- ifelse(den > 0, MSE/den^2, 1e+20)	#
#
## find global minimum of the GCV function on the grid
#
#
	gcv.grid <- data.frame(lambda.grid, trA, V, sqrt(RSS/(N - trA)))
	names(gcv.grid) <- c("lambda", "trA", "GCV", "shat")	#
# il is the index of the smallest value in the grid 
#
	il <- order(gcv.grid$GCV)[1]
	lambda.gcv <- gcv.grid$lambda[il]	#
#
#
	gcv.raw <- min(gcv.grid$GCV)
	if(verbose) {
		cat("GCV coarse search:", gcv.raw)
	}
# Now switch the 1e20 in V to NA's ( these are cases where the 
# the denominator is negative or zero due to the cost being greater than 1. 
#
	gcv.grid$GCV[den < 0] <- NA	#
#
#
#  create a mini tps object list with the information needed
# for further refinements of this estimate and the others
#
	info <- list(matrices = list(D = D, u = u), N = N, nt = nt, cost = cost,
		pure.ss = pure.ss)	#       
	if(verbose) print(info)	#
#
# do a golden section refined search for minimizing lamdda
# if the minimum is in interior of the grid search. 
#
	lambda.est <- matrix(ncol = 4, nrow = 3, dimnames = list(c("GCV", 
		"RMSE", "pure error"), names(gcv.grid)))
	lambda.est[1, 1] <- gcv.grid[il, 1]	#
#
	if((il > 1) & (il < nstep.cv)) {
#
# now do the Golden section refinement
# tolerance for convergence scaled with respect to GCV from the coarse search
#
		out <- golden.section.search(lambda.grid[il - 1], lambda.grid[
			il], lambda.grid[il + 1], tps.fgcv, f.extra = info, tol
			 = tol * gcv.raw)
		lambda.gcv <- out$x
		lambda.est[1, 1] <- lambda.gcv
	}
	else {
		warning("GCV search gives a minumum at the endpoints of the grid search"
			)
	}
	lambda.rmse <- NA
	lambda.pure.error <- NA
	if(!is.na(rmse)) {
		guess <- max(gcv.grid$lambda[gcv.grid$shat < rmse])
		if(verbose) {
			print(rmse)
			print(guess)
		}
		if(!is.na(guess)) {
			lambda.rmse <- find.upcross(tps.fs2hat, info, 
				upcross.level = rmse^2, guess = guess, tol = 
				tol * rmse^2)
			lambda.est[2, 1] <- lambda.rmse
		}
		else {
			warning("Value of rmse is outside possible range")
		}
	}
#
##
#
	if(!is.na(shat.pure.error)) {
		guess <- max(gcv.grid$lambda[gcv.grid$shat < shat.pure.error])
		if(!is.na(guess)) {
			lambda.pure.error <- find.upcross(tps.fs2hat, info, 
				upcross.level = shat.pure.error^2, guess = 
				guess, tol = tol * shat.pure.error^2)
			lambda.est[3, 1] <- lambda.pure.error
		}
		else {
			warning("Value of pure error estimate  is outside possible range"
				)
		}
	}
#
#
# fill in other stuff for each estimate of lambda
	for(k in 1:3) {
		lam <- lambda.est[k, 1]
		if(!is.na(lam)) {
			lambda.est[k, 2] <- tps.ftrace(lam, D)
			lambda.est[k, 3] <- tps.fgcv(lam, info)
			lambda.est[k, 4] <- sqrt(tps.fs2hat(lam, info))
		}
	}
	list(gcv.grid = gcv.grid, lambda.est = lambda.est, lambda.best = 
		lambda.gcv)
}
