"nnreg"<-
function(x, y, k1, k2, start, ngrind = 250, ntries = 100, npol = 20, tol1 = 
	9.9999999999999995e-07, tol2 = 1.0000000000000001e-09, itmax1 = 250, 
	itmax2 = 10000, derivative = F, fout = "nnreg.out", run = T, just.setup
	 = F, just.read = F, fitted.values = F, all.fits = F, greedy = F, seed, 
	clean = T)
{
	lags <- NA
	call <- match.call()
	y <- c(y)	# make sure y is just a vector!
	if(missing(seed))
		seed <- as.integer(runif(1) * 125000)
	if(is.list(x)) {
		k1 <- y
		k2 <- k1
		lags <- x$lags
		y <- x$y
		x <- x$x
	}
	x <- as.matrix(x)
	d <- ncol(x)
	jac.list <- NA
	if(!just.read) {
		write(t(cbind(y, x)), "nnreg.dat", ncol = 1)	#
# A negative grid number is the switch to indicate a single fit based
# on the start values. In this case only one specification of the
#hidden units makes sense.
#
		if(!missing(start)) {
			ngrind <- -1
		}
		if(all.fits)
			iprint <- 1
		else iprint <- 0
		if(greedy)
			igreed <- 1
		else igreed <- 0
		write(c("nnreg.dat", "nnreg.sum"), "nnreg.par", ncol = 1)
		temp <- c(length(y), ncol(x), ngrind, ntries, npol, iprint, 
			igreed, seed, tol1, tol2, itmax1, itmax2, k1, k2)
		write(temp, "nnreg.par", ncol = 1, append = T)
		if(ngrind < 0) {
#
# extract the parameters from the object passed as start. 
#
			if(class(start)[1] == "nnreg") start <- start$model[[
				  start$best.model]]$theta
			if(class(start)[1] == "netfit")
				start <- start$theta
			write(unlist(start), "nnreg.str", ncol = 1)
		}
		if(just.setup) {
			cat("Input file and data files have been constructed for nnreg",
				fill = T)
			cat("Run nnreg in the UNIX shell or in DOS by:     nnreg.x > outfile",
				fill = T)
			cat("the S data set FUNFITS has the path to nnreg.x", 
				fill = T)
			return()
		}
		if(run == T) {
			cat("Running nnreg in the shell", fill = T)
			remove.file(fout)
			#unix(paste(.Library,"/funfits/exec", "/nnreg.x  > ", fout, sep = "")
			unix(paste(system.file("exec/nnreg.x")," > ", fout, sep = "")
				)
		}
	}
# end of run and set block
	if(!just.setup) {
		temp <- list()
		cat("Reading in results from output file", fill = T)
		temp$model <- read.nnreg(fout)	# read in summary from LENNS
		class(temp) <- c("nnreg", "funfits")
		if(greedy)
			temp$model <- netfit.reformat.greedy(temp$model)
		temp$summary <- scan("nnreg.sum", what = "a", sep = "\n")
		class(temp$summary) <- "text"	#end
		nfits <- length(temp$model)
		cat(nfits, " models read in from ", fout, fill = T)
		temp$fitted.values <- matrix(NA, ncol = nfits, nrow = length(y)
			)
#                print(temp)
		if(clean) {
			remove.file("nnreg.dat")
			remove.file("nnreg.par")
			remove.file("nnreg.sum")
			remove.file("nnreg.str")
			remove.file(fout)
		}
		if(!all.fits | fitted.values) {
			for(k in 1:nfits) {
				temp$fitted.values[, k] <- predict(temp$model[[k
				  ]], x)
			}
			temp$residuals <- y - temp$fitted.values
		}
		temp$call <- call
		temp$x <- x
		temp$y <- y
		temp$n <- length(y)
		temp$nfits <- nfits
		temp$lags <- lags
		temp$seed <- seed
		hold <- summary(temp, noprint = T)
		temp$best.model <- order(hold[, 6])[1]
		return(temp)
	}
	else {
		invisible()
	}
}
