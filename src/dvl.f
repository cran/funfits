      double precision function dvl(lgnlam,svals,z,npsing)
      integer npsing
      double precision lgnlam,svals(npsing),z(npsing)
c
c Purpose: evaluate the cross-validation function with a semi-norm.
c
c On Entry:
c   lgnlam		log10(nobs*lambda) where lambda is the value of
c			lambda for which V is evaluated
c   svals(npsing)	singular values
c   z(npsing)		data vector in canonical coordinates
c   npsing		number of positive svals
c
c On Exit:
c   dvl			V(lambda)
c
c $Header: /usr/local/cvsroot/funfits/src/dvl.f,v 1.1.1.1 1998/05/24 21:50:08 agebhard Exp $
c
      integer j
      double precision nlam,numrtr,denom,factor
c
      common / gcvcom / addend,rss,tria,n,h
      integer n,h
      double precision rss,tria,addend
c     			see dvlop for definition of common block
c			variables
c
      nlam = 10**lgnlam
      numrtr = addend
      denom = dble(n - h - npsing)
      do 10 j = 1,npsing
         factor = 1.0d0/(1.0d0 + (svals(j)**2)/nlam)
         numrtr = numrtr + (factor*z(j))**2
         denom = denom + factor
  10  continue
      rss=numrtr
      tria=denom
 1000 format( e15.8)
      dvl=dble(n)*numrtr/denom**2
      return
      end
