"tps"<-
function(x, y, lambda = NA, df = NA, cost = 1, knots, weights = rep(1, length(y
	)), m, power, scale.type = "unit.sd", x.center, x.scale, 
	return.matrices = T, nstep.cv = 80, method = "GCV", rmse = NA, 
	link.matrix = NA, verbose = F, subset = NULL, tol = 0.0001, 
	print.warning = T, yname = NULL)
{
	out <- list()
	out$tag <- 1
	class(out) <- c("tps", "funfits")
	out$call <- match.call()	##
## S function to find minizier of 
##  || y- Xb||^2 + lambda b^T H b where H is a nonnegative definite 
## matrix
## Solution for b is  b= (X^T*X + lambda*H)^(-1) X^T*Y
##  the  H matrix is consructed to be the thin plate spline roughness
## matrix. (If the power =2m-d) 
##
## First set up some constants
## and some values in the output list
##
	x <- as.matrix(x)
	y <- c(y)	# make sure y is a vector!
	if(!is.null(subset)) {
		x <- x[subset,  ]
		y <- y[subset]
		out$subset <- paste(substitute(subset))
	}
	out$x <- x
	out$y <- y
	N <- length(y)
	out$N <- N
	lambda.est <- NA
	d <- ncol(x)	##
## make sure that 2m-d > 0
##
	out$form <- T
	with.constant <- T	## refers to weird constant for radial basis
	if(missing(m)) {
		m <- max(2, ceiling(d/2 + 0.10000000000000001))
	}
	if(missing(power)) {
		power <- 2 * m - d
		if(power < 1) {
			power <- 1
			out$form <- F
			if(print.warning)
				cat("Warning: Model is not a true thin plate spline",
				  fill = T)
		}
	}
## If not a true tps then do not find the  weird constant for the basis 
## functions
	if(2 * m - d <= 0) {
		with.constant <- F
	}
	if(2 * m - d != power) {
		with.constant <- F
	}
	out$cost <- cost
	out$m <- m
	out$with.constant <- with.constant
	out$trace <- NA
	if(is.null(yname))
		out$yname <- as.character(paste(substitute(y), collapse = ""))
	else out$yname <- yname
	out$weights <- weights	##
## Now find the estimate of sigma based on replicated points if this 
## makes sense
	rep.info <- cat.matrix(x)	## integer tags to indicate replications
	if(verbose)
		print(rep.info)
	if(max(rep.info) == N | !is.na(link.matrix[1])) {
		shat.rep <- NA
		shat.pure.error <- NA
	}
	else {
##
## do a simple 1-way ANOVA to get the rep error
##
		shat.pure.error <- sqrt(fast.1way(rep.info, y, weights)$MSE)
		shat.rep <- shat.pure.error
		out$shat.pure.error <- shat.pure.error
	}
	out$shat.rep <- shat.rep
	out$shat.pure.error <- shat.pure.error
	if(missing(knots))
		knots <- x[!dup(rep.info),  ]
	knots <- as.matrix(knots)
	out$knots <- knots	##
##
## scale the X's 
	x <- transformx(x, scale.type, x.center, x.scale)
	transform <- attributes(x)
	out$transform <- transform	## scale the knots int eh same way
	knots <- scale(knots, center = transform$x.center, scale = transform$
		x.scale)	##
#######################   NOTE        #############################
############ both the x and the knots must be scaled ################
################################################
##
	just.solve <- (lambda[1] == 0)
	if(is.na(just.solve))
		just.solve <- F
	out$power <- power	## make up the T and K matrices
## find the QR decopmposition of T matrix  that spans null space with
## respect to the knots 
	qr.T <- qr(make.tmatrix(knots, m))
	tmat <- make.tmatrix(x, m)
	out$ptab <- attributes(tmat)$ptab
	X <- cbind(tmat, qr.yq2(qr.T, make.rb(x, knots, power, with.constant = 
		with.constant)))
	if(verbose) print(dim(X))	
	## transform the X evalution matrix by a linear transformation if 
## the link matrix has been passed
##
	if(!is.na(link.matrix[1])) X <- link.matrix %*% X	##
	np <- ncol(X)	## the number of parameters
	nr <- nrow(X)
	N <- nr
	nt <- qr.T$rank	## number of para. in NULL space
	nk <- np - nt
	out$np <- np
	out$nt <- nt	
	##   construct the roughness penalty matrix  using radial basis
##functions and Qr decomposition of T
##
	H <- matrix(0, ncol = np, nrow = np)
	temp <- qr.yq2(qr.T, make.rb(knots, knots, power, with.constant = 
		with.constant))
	temp <- qr.q2ty(qr.T, temp)
	H[(nt + 1):np, (nt + 1):np] <- temp	##
## if lambda = 0 then just solve the system 
	if(just.solve) {
#
##  just find the least squares fit using radial basis functions or
## the interpolation if knots are missing. 
##
		out$method <- "interpolation"
		omega <- qr.coef(qr(X), y)
	}
	else {
##
##   do all the heavy decompositions if lambda is not = 0 
##   or if it is omitted
##
##
## inverse symetric square root of X^T W  X
##
		temp <- svd(sqrt(weights) * X)[c("v", "d")]	##
		if(max(temp$d)/min(temp$d) > 10000000000) {
			if(verbose)
				print(temp$d)
			print("Must use a reduced set of\nknots because the radial basis functions are close to being singular"
				)
			out <- NULL
			return(out)
		}
##
##
		B <- temp$v %*% diag(1/(temp$d)) %*% t(temp$v)	##
##   eigenvalue eigenvector decomposition of BHB
##
		temp <- svd(B %*% H %*% B)
		U <- temp$u
		D <- temp$d
		if(verbose) print(D)	
	##   We know that H has at least nt zero singular values ( see how H is
##   filled)
##   So make these identically zero.
##   the singular values are returned from largest to smallest.
##
		D[(1:nt) + (np - nt)] <- 0
		G <- B %*% U	##
##   with these these decompositions it now follows that 
##     b= B*U( I + lambda*D)^(-1) U^T * B * X^T*Y
##      = G*( I + lambda*D)^(-1) G^T* X^T*Y
##	
##
		u <- t(X %*% G) %*% (y * weights)	##
## find the (weighted) pure error sums of squares by calculating 
## predcited values when lambda=0 
## 
		temp1 <- (X %*% G) %*% u
		out$pure.ss <- sum(weights * (y - X %*% G %*% u)^2)
		out$matrices <- list(D = D, G = G, u = u, qr.T = qr.T)	##
## find some estimates of lambda
##
		gcv.out <- gcv.tps(out, cost = cost, nstep.cv = nstep.cv, rmse
			 = rmse, verbose = verbose, tol = tol)	##
		out$gcv.grid <- gcv.out$gcv.grid	##
##
		lambda.est <- gcv.out$lambda.est	##
		if(verbose) print(lambda.est)	##
## find the one specified by the method but first fill in a 
## possible user supplied value
##
##
		if(!missing(lambda) | !missing(df)) {
			method <- "user"	
	## is the df is supplied then find the corresponding lambda
			if(!is.na(df)) {
				lambda <- tps.df.to.lambda(df, D)
			}
			temp <- c(lambda, NA, NA, NA)
			lab <- c(dimnames(lambda.est)[[1]], "user")
			lambda.est <- rbind(lambda.est, temp)
			row.names(lambda.est) <- lab
		}
## find the best one. 
##
		lambda.best <- lambda.est[method, "lambda"]
		if(verbose) print(lambda.best)	##
## To solve for the coefficients,  recall: omega= G*( I + lambda*D)^(-1)*u
## predicted values are X%*% omega
		omega <- G %*% ((1/(1 + lambda.best * D)) * u)
	}
	if(!just.solve) {
		out$eff.df <- sum(1/(1 + lambda.best * D))
		out$trA2 <- sum(1/(1 + lambda.best * D)^2)
		temp <- X %*% out$matrices$G %*% sqrt(diag(1/(1 + lambda.best * 
			out$matrices$D)))
		diagA <- c((temp^2) %*% rep(1, ncol(X))) * out$weights
		out$diagA <- diagA
	}
	if(just.solve)
		out$eff.df <- out$np
	out$fitted.values <- c(X %*% omega)
	out$residuals <- y - c(X %*% omega)
	out$trace <- out$eff.df	##
	if(verbose)
		print(out$eff.df)
	if(just.solve) {
		out$lambda <- lambda
		out$gcv.grid <- matrix(c(lambda, rep(NA, 4)), nrow = 1)
	}
	else {
		out$lambda <- lambda.best
		out$method <- method
	}
	out$best.model <- out$lambda
	out$omega <- omega
	out$d <- omega[1:nt]	
	## transform the omegas associated with the radial basis functions back
##  into the c parameter vector. 
## 
##
	temp <- omega
	temp[1:nt] <- 0
	out$c <- c(qr.qy(qr.T, temp))
	out$coefficients <- c(omega[1:nt], out$c)
	out$just.solve <- just.solve
	res.df <- (N - out$trace)	##
##
## find an estimate of the residual standard deviation
## based on fitted spline
	if(res.df > 0) {
		out$GCV <- (sum(out$residuals^2 * weights)/N)/(1 - out$eff.df/N
			)^2
		out$shat <- sqrt(sum(out$residuals^2 * weights)/(res.df))
		if(method == "user") {
## fill in the info for the lambda.est data frame
## for the user supplied value of lambda
			lambda.est["user", 2] <- out$eff.df
			lambda.est["user", 3] <- out$GCV
			lambda.est["user", 4] <- out$shat
		}
	}
	else {
		out$shat <- 0
		out$GCV <- NA
	}
	if(verbose) {
		print("shat")
		print(out$shat)
	}
	if(is.na(link.matrix[1])) {
		r.square <- cor(out$fitted.values * sqrt(out$weights), (out$y) * 
			sqrt(out$weights))^2
		out$r.square <- r.square	## calculate the q2 value
		if(!just.solve) {
#
# don't do this if interplotating
#
			cv.res <- out$residuals/(1 - diagA)
			press <- sum((cv.res)^2)
			rmse.press <- (press/length(cv.res))^0.5
			ss.tot <- sum((y - mean(y))^2)
			q2 <- (ss.tot - press)/ss.tot
			out$q2 <- q2
			out$press <- press
			out$rmse.press <- rmse.press
		}
	}
	else {
		out$press <- NA
		out$rmse.press <- NA
		out$q2 <- NA
	}
	out$lambda.est <- lambda.est
	out$best.model <- out$lambda	#
#zap matrices if no return
#
	if(!return.matrices)
		out$matrices <- NA
	out
}
