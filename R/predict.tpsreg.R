"predict.tpsreg"<-
function(out, grid, clean = T)
{
	if(missing(grid)) grid <- out$x	
	#	cat("Writing grid out to an input file for ev.x", fill = T)
	write(t(grid), "tpsgrid.in", ncol = 1)
	write(out$coefficients, "tps.ev", ncol = 1)	
	#	cat("Running ev.x to evalute spline on  grid", fill = T)
	#unix(paste(FUNFITS.BIN, "/tpsregev.x ", " < tpsgrid.in > tpsev.out", 
	unix(paste(system.file("exec/tpsregev.x"), " < tpsgrid.in > tpsev.out", 
		sep = ""))	# read output of ev.x back into S
	fit <- scan("tpsev.out")
	if(clean) {
		remove.file("tpsgrid.in")
		remove.file("tpsev.out")
		remove.file("tps.ev")
	}
	fit
}
