"make.tmatrix"<-
function(x, m = 2)
{
	if(!is.loaded(symbol.For("radbas"))) {
	#	temp <- dyn.load(paste(FUNFITS.BIN, "TPS.o", sep = "/"), 2)
          stop("dynamic code not loaded!")
	}
	d <- ncol(x)	#
	n <- nrow(x)	#cat(d, n, fill = T)	
#      subroutine dmaket(m,n,dim,des,lddes,npoly,t,ldt, wptr,info)
#      integer m,n,dim,lddes,npoly,ldt,wptr(dim),info
#      double precision des(lddes,dim),t(ldt,*)
#   m			order of the derivatives in the penalty
#   n			number of rows in des
	nterms <- .Fortran("mkpoly",
		as.integer(m),
		as.integer(d),
		nterms = as.integer(0))$nterms
	temp <- .Fortran("dmaket",
		m = as.integer(m),
		n = as.integer(n),
		dim = as.integer(d),
		des = as.double(x),
		lddes = as.integer(n),
		npoly = as.integer(nterms),
		tmatrix = as.double(rep(0, n * (nterms))),
		ldt = as.integer(n),
		wptr = as.integer(rep(0, d * m)),
		info = as.integer(0),
		ptab = as.integer(rep(0, nterms * d)),
		ldptab = as.integer(nterms))
	temp2 <- matrix(temp$tmatrix, nrow = n)
	attr(temp2, "ptab") <- matrix(temp$ptab, nrow = nterms, ncol = d)
	temp2
}
