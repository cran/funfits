      double precision function dtrace(lgnlam,svals,npsing)
      integer npsing
      double precision lgnlam,svals(npsing)
c
c Purpose: compute the trace of A 
c
c On Entry:
c   lgnlam		log10(nobs*lambda) where lambda is the value of
c			lambda for which V is evaluated
c   svals(npsing)	singular values
c   npsing		number of positive svals
c
c On Exit:
c   dtrace			tr(A(lambda))
c
c $Header: /usr/local/cvsroot/funfits/src/dtrace.f,v 1.1.1.1 1998/05/24 21:50:08 agebhard Exp $
c
      integer j
      double precision nlam,denom,factor
c
      common / gcvcom / addend,rss,tria,n,h
      integer n,h
      double precision rss,tria,addend
c     			see dvlop for definition of common block
c			variables
c
      nlam = 10**lgnlam
      denom = dble(n - h - npsing)
      do 10 j = 1,npsing
         factor = 1.0d0/(1.0d0 + (svals(j)**2)/nlam)
         denom = denom + factor
   10 continue
      dtrace=  dble(n- denom)
      return
      end
