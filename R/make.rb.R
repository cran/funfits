"make.rb"<-
function(x1, x2, p = 1, with.log = T, with.constant = T)
{
	if(!is.loaded(symbol.For("radbas"))) {
#		temp <- dyn.load(paste(FUNFITS.BIN, "TPS.o", sep = "/"), 2)
          stop("Compiled code has not been dynamically loaded")
	}
	if(!is.matrix(x1))
		x1 <- as.matrix(x1)
	if(!is.matrix(x2))
		x2 <- as.matrix(x2)
	d <- ncol(x1)
	n1 <- nrow(x1)
	n2 <- nrow(x2)	# 2m-d =p so m= ( d+p)/2
	m <- (d + p)/2
	par <- c(p/2, ifelse((d %% 2 == 0) & (with.log), 1, 0))
	temp <- .Fortran("radbas",
		nd = as.integer(d),
		x1 = as.double(x1),
		n1 = as.integer(n1),
		x2 = as.double(x2),
		n2 = as.integer(n2),
		par = as.double(par),
		k = as.double(rep(0, n1 * n2)))
	if(with.constant) {
		Amd <- radbas.constant(m, d)
	}
	else {
		Amd <- 1
	}
	Amd * matrix(temp$k, ncol = n2, nrow = n1)
}
