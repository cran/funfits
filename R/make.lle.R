"make.lle"<-
function(jac, nprod, statevector = F, verbose = T, clean = T)
{
	nc <- ncol(jac)
	remove.file("lle.par")
	remove.file("temp.lle")
	remove.file("lle.out")
	if(verbose)
		cat(" removed UNIX temp files", fill = T)
	if(statevector) {
		state <- 1
		emd <- sqrt(nc)
		test <- emd - round(emd)
		if(test != 0)
			cat("Stop, no. of col of jac is not a perfect \nsquare",
				fill = T)
	}
	else {
		state <- 0
		emd <- nc
	}
	write(c(nc, emd, nprod, state), "lle.par")
	write(t(jac), "temp.lle", ncol = nc)
	#unix(paste(FUNFITS.BIN, "/lle.x < temp.lle > lle.out", sep = ""))
	unix(paste(system.file("exec/lle.x")," < temp.lle > lle.out", sep = ""))
	if(verbose)
		cat(" reading in matrix of LLE's", fill = T)
	temp <- matrix(scan("lle.out"), ncol = 3, byrow = T)
	if(clean) {
		remove.file("lle.par")
		remove.file("temp.lle")
		remove.file("lle.out")
		remove.file("lle.warnings")
	}
	temp
}
