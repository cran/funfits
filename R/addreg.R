"addreg"<-
function(x, y, lam, nback = 20, tol = 1.0e-05, start, cost = 1)
{
#	if(!is.loaded(symbol.For("addreg"))) {
#		temp <- dyn.load("addreg.sub.o", 2)
#		cat("addreg  FORTRAN subroutines dynamically loaded:", temp, 
#			fill = T)
#	}
	if(!is.loaded(symbol.For("css"))) {
          stop("Compiled code has not been dynamically loaded")
	}
	if(!is.loaded(symbol.For("addreg"))) {
          stop("Compiled code has not been dynamically loaded")
        }
	out <- NULL
	x <- as.matrix(x)
	M <- ncol(x)
	N <- length(y)
	maxit <- 8
	nstep <- 10
	if(nrow(x) != length(y))
		stop(" X and Y do not match")
	if(!missing(lam)) {
		if(length(lam) == 1) {
			lam <- rep(lam, M)
		}
		lam <- ifelse(is.na(lam), -1, lam)
	}
	if(missing(lam)) {
		lam <- rep(-1, M)
	}
	if(missing(start)) {
		use.start <- 0
		start <- rep(0, N)
	}
	else {
		use.start <- 1
	}
	gcv.grid <- matrix(0, nrow = 4, ncol = M)
	a <- .Fortran("addreg",
		x = as.double(x),
		as.integer(N),
		as.integer(M),
		as.integer(N),
		y = as.double(y),
		as.double(rep(1, N)),
		as.double(lam),
		trace = as.double(rep(0, M)),
		sxy = as.double(matrix(0, ncol = M, nrow = N)),
		dsxy = as.double(matrix(0, ncol = M, nrow = N)),
		sy = as.double(start),
		din = as.double(c(tol, nback, maxit, nstep, cost)),
		dout = as.double(rep(0, 2 + M)),
		job = as.integer(use.start),
		ierr = as.integer(0))
	if(a$ierr != 0) {
		cat("Error in call to addreg")
		return(a)
	}
	if(sum(a$trace) > length(y))
		cat(" WARNING: Effective number of parameters exceeds the nmbe of observations",
			fill = T)
	out$x <- matrix(a$x, ncol = M)
	out$y <- a$y
	out$residuals <- a$y - a$sy
	out$fitted.values <- a$sy
	out$predicted.comp <- matrix(a$sxy, ncol = M)
	out$trace <- a$trace
	out$lambda <- a$dout[(1:M) + 2]
	out$converge <- c(a$dout[1:2])	
	#	list(x = matrix(a$x, ncol = M), y = a$y, residuals = a$y - a$sy, 
#	predicted = a$sy, predicted.comp = matrix(a$sxy, ncol = M), 
#	predicted.comp.d = matrix(a$dsxy, ncol = M), trace = a$trace, 
#	lambda = lam, converge = c(a$dout[1:2]))
	class(out) <- "addreg"
        out
}
