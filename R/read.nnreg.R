"read.nnreg"<-
function(fname = "nnreg.out")
{
	temp <- scan(fname)
	ntemp <- length(temp)
	out <- list()	#	cmd <- paste("grep '#' ", fname, " > nnreg.temp")
#	unix(cmd)
#	out$summary <- c(scan("nnreg.temp", what = "a", sep = "\n"), 
#scan("nnreg.summary", what = "a", sep = "\n"))
	i <- 1
	loc <- 0
	out <- list()
	while(loc < ntemp) {
#print(" start model loc=")
#print(loc)
		if(loc + 3 > ntemp) break
		d <- temp[1 + loc]
		k <- temp[2 + loc]
		rms <- temp[3 + loc]
		loc <- loc + 3
		if(loc + 2 * (d + 1) > ntemp)
			break
		xm <- temp[(1:d) + loc]
		loc <- loc + d
		xsd <- temp[(1:d) + loc]
		loc <- loc + d
		ym <- temp[loc + 1]
		ysd <- temp[loc + 2]
		loc <- loc + 2
		np <- 1 + k * (d + 2)
		if(loc + np > ntemp)
			break
		theta <- temp[loc + (1:np)]
		loc <- loc + np
		out[[i]] <- list(d = d, k = k, xm = xm, ym = ym, xsd = xsd, ysd
			 = ysd, np = np, theta = theta, rms = rms)
		class(out[[i]]) <- "netfit"
		i <- i + 1	#	print(loc)
	}
	if(loc != length(temp))
		cat("incomplete information in ", fname, " output file", fill
			 = T)
	out
}
