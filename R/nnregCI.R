"nnregCI"<-
function(fit, model = fit$best.model, ngrind = 250, ntries = 100, npol = 20, 
	clevel = 0.94999999999999996, cut1 = NA, cut2 = NA, nfits = 500, tol1
	 = 9.9999999999999995e-07, tol2 = 1.0000000000000001e-09, itmax1 = 250, 
	itmax2 = 10000, fdata, fout = "nnci.out", seed)
{
	call <- match.call()
	if(missing(seed))
		seed <- as.integer(runif(1) * 125000)
	y <- fit$y
	x <- as.matrix(fit$x)
	k <- fit$model[[model]]$k
	d <- ncol(x)
	rms <- fit$model[[model]]$rms
	np <- fit$model[[model]]$np
	if(missing(fdata)) {
		fdata <- "nnci.dat"
		write(t(cbind(y, x)), "nnci.dat", ncol = 1)
	}
	if(is.na(cut1))
		cut1 <- rms * sqrt((1 + qf(clevel, np, length(y) - np) * (np/(
			length(y) - np))))
	if(is.na(cut2))
		cut2 <- cut1 - 0.20000000000000001 * (cut1 - rms)
	remove.file("nnci.sum")
	write(c(fdata, "nnci.sum"), "nnci.par", ncol = 1)
	temp <- c(length(y), ncol(x), ngrind, ntries, npol, rms, cut1, cut2, 
		nfits, seed, tol1, tol2, itmax1, itmax2, k)
	write(temp, "nnci.par", ncol = 1, append = T)
	cat("Running nnregci in the shell", fill = T)
	remove.file(fout)
	#unix(paste(FUNFITS.BIN, "/nnregci.x  > ", fout, sep = ""))	
        unix(paste(system.file("exec/nnregci.x")," > ", fout, sep = ""))
	# end of run and set block
	temp <- list(model=NULL)
	class(temp) <- c("nnreg", "funfits")
	cat("Reading in results from output file", fill = T)
	temp$model <- read.nnreg(fout)	# read in summary from LENNS
	temp$summary <- scan("nnci.sum", what = "a", sep = "\n")
	class(temp$summary) <- "text"
	nfits <- length(temp$model)
	cat(nfits, " models read in from ", fout, fill = T)
	temp$call <- call
	temp$x <- x
	temp$y <- y
	temp$n <- length(y)
	temp$nfits <- nfits
	temp$seed <- seed
	return(temp)
}
