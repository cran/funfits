"FEXP"<-
function(x1, x2, C, theta = rep(1, ncol(x1)), load = T)
{
	if(!is.loaded(symbol.For("m2deb")) | load) {
          stop("Compiled code has not been dynamically loaded")
          #temp2 <- dyn.load(paste(FUNFITS.BIN.NEW, "FEXP.o", sep = "/"), 
          #                  2)
	}
	if(length(theta) == 1)
		theta <- rep(theta, ncol(x1))
	d <- ncol(x1)
	n1 <- nrow(x1)
	n2 <- nrow(x2)	# scale both X's to reflect theta parameter
	x1 <- t(t(x1/theta))
	x2 <- t(t(x2/theta))
	par <- rep(1, d)
	if(!is.matrix(C)) {
		C <- matrix(C, ncol = 1)
	}
	nc <- ncol(C)
	if(nrow(C) != n2) stop("number of rows in C does not match\nnumber of locations in x2"
			)	#8
# return the covaraince matrix multiplied by the vector C
#
	matrix(.Fortran("m2deb",
		nd = as.integer(d),
		x1 = as.double(x1),
		n1 = as.integer(n1),
		x2 = as.double(x2),
		n2 = as.integer(n2),
		par = as.double(par),
		c = as.double(C),
		nc = as.integer(nc),
		h = as.double(rep(0, n1 * nc)),
		work = as.double(rep(0, nc)))$h, ncol = nc)
}
