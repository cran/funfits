c**** driver to read in Jacobians and compute LLE
      parameter(NMAX=20000, MAXL=20)
      real*8 der(MAXL,  NMAX)
      real*8 llesv(MAXL,NMAX), lleqr(MAXL,NMAX), lle11(NMAX)
      integer state
      open(unit=2,file='lle.par',status='old')
      open(unit=3,file='lle.warnings')
      read( 2,*) nc, nlags, nprod, state
      if( nc.gt.MAXL) then
          write(3,*) "number of columns of der too large!"
          stop
      endif

      ldder=MAXL

      do 10 k=1,NMAX
      read(*,*,end=20) (der(il,k), il=1,nc)

 10   continue
 
      write(3,*) "Maximum number of rows have been read"
      stop
 20   continue
      n=k-1
      if( nprod.lt.0)  then
      nprod= n
      endif
      if( nprod.gt.n) then 
         write(3,*) " number of steps is larger than series length"
         stop
       endif 
      call llest( nlags,ldder,der,n,nprod,nll,llesv,lleqr,lle11,state)

      do 500 k=1,nll
      write(*, 1000) llesv(1,k),lleqr(1,k),lle11(k)
 1000 format( 5e16.7)
 500  continue
      stop  
      end
