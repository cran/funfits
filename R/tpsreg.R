"tpsreg"<-
function(x, y, lambda, m, clean = T)
{
### run the tps program while in S
### data file is tps.dat
	out <- list()
	out$call <- match.call()
	class(out) <- "tpsreg"
	if(missing(lambda)) {
		gcv <- 1
	}
	else gcv <- 101
	ncov <- 0	## swh scale x matrix according to user specified options
	x <- as.matrix(x)
	d <- ncol(x)	#print(d)
### pdh 10/12/94 - set a better default for m
	if(missing(m))
		m <- ceiling((d + 2)/2)
	write(t(cbind(x, y)), "tps.dat", ncol = (d + 1))
	write(c(d, m, ncov, gcv), "tps.in")
	if(!missing(lambda))
		write(lambda, "tps.in", append = T)
	write("tps.dat", "tps.in", append = T)
	cat("Running thinplate spline program in the shell", fill = T)
	#TPSREG <- paste(FUNFITS.BIN, "tpsreg.x", sep = "/")
	#unix(paste(TPSREG, "< tps.in  > tps.sum"))
        #unix(paste(.Library,"/funfits/exec", "/nnreg.x  > ", fout, sep = ""))
	unix(paste(system.file("exec/tpsreg.x")," < tps.in > tps.sum",
                   sep = ""))
	cat(" Output from tps is in the file  tps.sum", fill = T)	
	## nobs,dim,m,ncov1,(iout(k),k=1,4)
	parms <- scan("tps.par")
	dim <- parms[2]
	spar <- parms[9]	#	print(parms)
	cat(" Reading in estimated spline at the data points", fill = T)
	tps.spline <- matrix(scan("tps.spl"), ncol = 2, byrow = T)
	gcvf <- scan("tps.gcv")
	gcvf <- matrix(gcvf, ncol = 3, byrow = T)	# add to the outout  object 
#out$x <- tps.spline[, 1:dim]
## pdh 5/25/94
	out$x <- x
	out$y <- y
	out$yname <- paste(substitute(y))
	out$summary <- scan("tps.sum", what = "a", sep = "\n")
	class(out$summary) <- "text"
	out$residual <- tps.spline[, 2]
	out$fitted.values <- tps.spline[, 1]
	out$gcv.grid <- gcvf
	out$parameters <- parms	
	##         read in the coefficients of the  spline          
	out$coefficients <- scan("tps.ev")
	out$spar <- spar
	out$eff.df <- parms[10]
	if(clean) {
		remove.file("tps.ev")
		remove.file("tps.spl")
		remove.file("tps.sum")
		remove.file("tps.gcv")
		remove.file("tps.in")
		remove.file("tps.dat")
		remove.file("tps.in")
	}
	out
}
