"krig"<-
function(x, Y, cov.function = exp.cov, lambda = NA, cost = 1, knots, weights = 
	rep(1, length(Y)), m = 2, return.matrices = T, nstep.cv = 80, 
	scale.type = "user", x.center = rep(0, ncol(x)), x.scale = rep(1, ncol(
	x)), rho = NA, sigma2 = NA, method = "GCV", verbose = F, ...)
{
  out<-list(call=NULL)
	out$call <- match.call()	#
	out$N <- length(Y)	# reshuffle arguments if x nd y are passed as a list
	class(out) <- c("krig", "funfits")
#
###
### add passed theta arguments to the covariance function
	temp <- list(...)
	ntemp <- names(temp)
	if(length(temp) > 0) {
		for(k in 1:length(ntemp)) {
			cov.function[ntemp[k]] <- temp[ntemp[k]]
		}
		cov.function <- as.function(cov.function)
	}
	out$cov.function <- cov.function	#
# S function to find minizier of 
#  || Y- Xb||^2 + lambda b^T H b where H is a 
# covariance matrix found from cov.function
# Solution for b is  b= (X^T*X + lambda*H)^(-1) X^T*Y
#  H is the covariance matrix 
# First set up some constants
	x <- as.matrix(x)
	Y <- c(Y)	# make sure Y is a vector!
	out$y <- Y
	out$x <- x
	out$weights <- weights	
	## if knots are missing then use the set of unique x vexctors.
	if(missing(knots))
		knots <- x[!dup.matrix(x),  ]
	knots <- as.matrix(knots)	##
	out$knots <- knots	##
#
## scale x and knots 
	x <- transformx(x, scale.type, x.center, x.scale)
	transform <- attributes(x)
	knots <- scale(knots, center = transform$x.center, scale = transform$
		x.scale)
	out$transform <- transform	##
##
##  use value of lambda implied by rho and simga2 if these are passed
##
	if(!is.na(lambda))
		method <- "user"
	if(!is.na(rho) & !is.na(sigma2)) {
		lambda <- sigma2/rho
		method <- "user"
	}
	just.solve <- (lambda[1] == 0)
	if(is.na(just.solve))
		just.solve <- F
	d <- ncol(x)	# make up the T and K matrices
# find the QR decopmposition of T matrix  that spans null space with
# respect to the knots 
	qr.T <- qr(make.tmatrix(knots, m))
	X <- cbind(make.tmatrix(x, m), qr.yq2(qr.T, cov.function(x, knots)))
	np <- ncol(X)	# the number of parameters
	nr <- nrow(X)
	N <- nr
	nt <- qr.T$rank	# number of para. in NULL space
	nk <- np - nt	#
	out$np <- np
	out$nt <- nt	#   construct the covariance matrix 
#functions and Qr decomposition of T
#
	H <- matrix(0, ncol = np, nrow = np)
	temp <- qr.yq2(qr.T, cov.function(knots, knots))
	temp <- qr.q2ty(qr.T, temp)
	mean.var <- mean(diag(temp))
	H[(nt + 1):np, (nt + 1):np] <- temp	#
#
# if lambda = 0 then just solve the system 
	if(just.solve) {
		beta <- qr.coef(qr(X), Y)
	}
	else {
#
#   do all the heavy decompositions if lambda is not = 0
#   or if it is omitted
#
#
# inverse symetric square root of X^T W  X
#
		temp <- svd(diag(sqrt(weights)) %*% (X))[c("v", "d")]
		cond.matrix <- max(temp$d)/min(temp$d)
		if(cond.matrix > 10000000000)
			stop("Covarinace matrix is clsoe\nto singular")
		B <- temp$v %*% diag(1/(temp$d)) %*% t(temp$v)	#
#   eigenvalue eigenvector decomposition of BHB
#
		temp <- svd(B %*% H %*% B)
		U <- temp$u	#	cat(diag(U %*% t(U)), fill = T)
		D <- temp$d	#
		if(verbose) {
			cat("singular values:", fill = T)
			print(D)
		}
#   We know that H has atleast nt zero singular values ( see how H is
#   filled)
#   So make these identically zero.
#   the singular values are returned from largest to smallest.
#
		D[(1:nt) + (np - nt)] <- 0
		G <- B %*% U	#
#   with these these decompositions it now follows that 
#     b= B*U( I + lambda*D)^(-1) U^T * B * X^T*Y
#      = G*( I + lambda*D)^(-1) G^T* X^T*Y
#	
# Now tranform  Y based on this last equation
#
		u <- t(G) %*% t(X) %*% (weights * Y)	#
#
#   So now we have   
#
#    b= G*( I + lambda*D)^(-1)*u 
#   Note how in this form we can rapidly solve for b for any lambda
#
# save matrix decopositions in out list
#
# find the pure error sums of sqaures. 
#
		out$pure.ss <- sum(weights * (Y - X %*% G %*% u)^2)
		if(verbose) {
			cat("pure.ss", fill = T)
			print(out$pure.ss)
		}
		out$matrices <- list(B = B, U = U, u = u, D = D, G = G, qr.T = 
			qr.T)
		gcv.out <- gcv.krig(out, nstep.cv = nstep.cv, verbose = verbose
			)
		gcv.grid <- gcv.out$gcv.grid	
	# To solve for the coefficients,  recall: b= G*( I + lambda*D)^(-1)*u
		if(method == "user") {
			lambda.best <- lambda
		}
		else {
			lambda.best <- gcv.out$lambda.best
		}
		beta <- G %*% ((1/(1 + lambda.best * D)) * u)
		out$gcv.grid <- gcv.grid
		if(verbose) {
			print(out$gcv.grid)
		}
	}
	out$cost <- cost
	out$m <- m
	if(!just.solve) {
		out$eff.df <- sum(1/(1 + lambda.best * D))
	}
	else {
		out$eff.df <- out$np
	}
	out$fitted.values <- c(X %*% beta)
	out$residuals <- Y - out$fitted.values
	if(just.solve)
		out$lambda <- lambda
	else out$lambda <- lambda.best
	out$yname <- substitute(Y)	##
##
## wipe out big matrices if they are not to be returned
##
##
## add some more stuff to put object
##
	out$beta <- beta
	out$d <- beta[1:nt]	#
# funny conversions are in case nt is equal to 1 and X is just a vector
#
	out$fitted.values.null <- as.matrix(X[, 1:nt]) %*% out$d	#
##
#
# tranform the beta into the parameter associated with the covariance
# function
# basis set. 
#  into the c parameter vector. 
#
	out$trace <- out$eff.df
	if(verbose) {
		cat("trace of A", fill = T)
		print(out$trace)
	}
	temp <- c(rep(0, nt), beta[(nt + 1):np])
	if(verbose)
		print(temp)
	out$c <- c(qr.qy(qr.T, temp))	
	#	out$coefficients <- c(beta[1:nt], out$c)
	if(verbose)
		print(out$c)
	out$just.solve <- just.solve	#
	out$shat.GCV <- sqrt(sum(out$weights * out$residuals^2)/(length(Y) - 
		out$trace))	#
# fill in the linear parameters of the covariance function in 
# the output object
#
# the next formula is pretty strange. It follows from solving the
# system of equations for the basis coefficients. 
#       
	out$rhohat <- sum(out$c * out$y)/(N - nt)	#
	if(is.na(rho)) {
		out$rho <- out$rhohat
	}
	else out$rho <- rho
	if(is.na(sigma2))
		sigma2 <- out$rho * out$lambda
	out$sigma2 <- sigma2	#
	out$shat.MLE <- sqrt(out$rhohat * out$lambda)	#
	out$best.model <- c(out$lambda, out$sigma2, out$rho)
	if(!return.matrices) {
		out$x <- NULL
		out$y <- NULL
		out$matrices <- NULL
		out$weights <- NULL
	}
##
##
	out
}
