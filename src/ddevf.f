      subroutine ddevf(des,lddes,nobs,dim,m,desb,lddesb,ndesb,
     * ncov1,ncov2,coef,npar,f,work,lwa,iwork,info,nder)
      integer lddes,nobs,dim,m,lddesb,ndesb,ncov1,ncov2,npar,lwa,
     * iwork(dim),info
      double precision des(lddes,dim),desb(lddesb,dim),
     * coef(npar),f(lddes),work(lwa)
      double precision  dcon(10)
c
c   Purpose: determine values of spline or its derivative
c
c  On Entry:
c   des(lddes,dim) 	design for the variables to be splined
c   lddes		leading dimension of des as declared in the
c			calling	program
c   nobs		number of rows in des
c   desb(lddesb,dim) 	locations for the basis functions
c			returned from dtpss and dptpss in the
c			variable des
c   lddesb		leading dimension of desb as declared in the
c			calling	program
c   ndesb		number of rows in desb
c   dim			number of columns in des
c   m			order of the derivatives in the penalty
c			calling	program
c   ncov1		number of covariates which duplicate the
c			replication structure of des
c   coef(npar)		coefficient estimates  [delta':xi']'
c   nder                order of derivatives to be computed
c                        ( only availble when dim=1,for first derivative, m=2)
c
c On Exit:
c   f(nobs)		predicted values
c   info		error indicator
c			  0 : successful completion
c			  1 : error in npar,ncov1,ncov2,m or dim
c			  2 : lwa too small
c			  3 : error in dmaket
c			
c
c Working Storage:
c   work(lwa)		double precision work vector
c   lwa			length of work vector
c			must be at least nobs*(nct+ndesb)
c			where nct = (m+dim-1 choose dim)
c   iwork(dim)		integer work vector
c
c Subprograms Called Directly:
c    Gcvpack - dmaket dmakek
c    Blas    - ddot
c
c Subprograms Called Indirectly:
c    Gcvpack - mkpoly
c    Blas    - dcopy
c    Other   - fact
c
c $Header: /usr/local/cvsroot/funfits/src/ddevf.f,v 1.1.1.1 1998/05/24 21:50:07 agebhard Exp $
c
      double precision dummy
      double precision ddot
      integer i,nct,p1,p1p1,npoly
      integer mkpoly
c
      ncov2=0
      nct = mkpoly(m,dim)
      if (npar .ne. ndesb + nct +ncov1+ncov2) then
         info = 1
         return
      endif
      if (lwa .lt. nobs*(nct+ndesb)) then
	 info = 2
	 return
      endif
c			first nobs*nct positions of work contain t
      p1 = nobs*nct
c			next nobs*ndesb postions of work contain k
      p1p1 = p1 + 1
c
      call dmaket(m,nobs,dim,des,lddes,dummy,1,0,npoly,work(1),nobs,
     * iwork,info)
      if (info .ne. 0) then
         info = 3
         return
      endif
      call dmakek(m,nobs,dim,des,lddes,ndesb,desb,lddesb,work(p1p1),
     * nobs)
c
c derivative of basis functions
      if(nder.gt.0 )  then 
            if( (dim.eq.1) ) then 
c            alter T and K matrices
             do 100 k=1,nobs
                work(k)=0.0
                work(nobs+k)=1.0
  100        continue
c
             index= p1p1
             do 110 i=1,ndesb
                do 120 j=1,nobs
                       diff= (des(j,1)-desb(i,1))
                       if ( abs(diff).lt.1.0d-8) then 
                        work(index)=0.0
                       else
                        work(index)= (2*m-1)*work(index) /diff
                       endif
                       index= index + 1
  120           continue
  110        continue


                else
           write(*,*) 'derivative only available for m=2 dim=1'
           stop
           endif
      endif
c			compute spline
      do 10 i = 1,nobs
	     f(i) = ddot(nct,coef,1,work(i),nobs)
	     f(i) = f(i) + ddot(ndesb,coef(nct+ncov1+ncov2+1),1,
     *		work(p1+i),nobs)
   10 continue
      return
      end
 

