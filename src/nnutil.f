       subroutine getxy(fname,xmat,y,nxmax,nx,nxc)
       implicit double precision(a-h,o-z)
       character*35 fname
       dimension xmat(nxmax,10), y(nxmax)
       integer nx,nxc
       open(unit=11,file=fname,status='old')
       rewind(11)
       m=0
       do 50 i=1,nxmax
           read(11,*,end=60)   y(i), (xmat(i,j), j=1,nxc)
           m=m+1
50      continue
60     continue
       close(11)
       if (m.lt.nxmax)  then 
        nx=m
       else 
        nx=nxmax
       endif
       return
       end

       subroutine canpar( theta, nu,nxc,ctheta)
       implicit  double precision (a-h, o-z)
       integer  ind(500),loc1,loc2,loc0   
       real*8 theta(1)  
       real*8 ctheta(*)   , sgn(500)
            np= nu*(nxc+2) +1
        if( np.gt.500) then
        write(*,*) " too m any parameters for canpar"
        stop
        endif  
       ctheta(1) = theta(1)      
        do  j=1,nu       
          ind(j)=j       
          ctheta(j+ 1)=  abs(theta(j+1))  
          if( theta(j+1).le.0)  then 
                  sgn(j)=-1
          else
                sgn(j)=1
          endif
        enddo   
                    loc1= 1+nu 
         call msort( ctheta( 2),ind,nu)  
             do j=1, nu     
     
             jp= ind(j)          
               ctheta(j +loc1)= theta( jp+ loc1)*sgn(jp)   
           enddo 

              loc0= 2*nu+1
            do j=1, nu       
             jp= ind(j)          
               loc1=  (jp-1)* nxc + loc0
              loc2=  (j-1)*nxc + loc0
            do jj=1,nxc               
               ctheta(loc2+jj)= theta( loc1+ jj)*sgn(jp)    
            enddo           
         enddo            
      return 
       
      end

        subroutine msort( k,ki,n)
C  HEAPSORT ALGORITHM FOR SORTING ON VECTOR OF KEYS K OF LENGTH N    
C  J F MONAHAN        TRANSCRIBED FROM KNUTH, VOL 2, PP 146-7.      
C integer array ki is permuted along with K
 
      REAL*8 K(1),KK 
      integer ki(1),kki
      INTEGER R             
      IF(N.LE.1) RETURN    
      L=N/2+1             
      R=N                
  2   IF(L.GT.1) GO TO 1
      KK=K(R)
      kki= ki(R)               
      K(R)=K(1)
      ki(R)=ki(1)  
      R=R-1       
      IF(R.EQ.1) GO TO 9   
      GO TO 3             
  1   L=L-1              
      KK=K(L)
      kki=ki(L) 
  3   J=L       
  4   I=J       
      J=2*J     
      IF(J-R) 5,6,8   
  5   IF(K(J).LT.K(J+1)) J=J+1   
  6   IF(KK.GT.K(J)) GO TO 8    
  7   K(I)=K(J)
      ki(I)=ki(J) 
      GO TO 4     
  8   K(I)=KK    
      ki(I)=kki 
      GO TO 2  
  9   K(1)=KK 
      ki(1)=kki   
      RETURN     
      END       

       subroutine tinit(theta,npar,eps,nreps)
       implicit double precision(a-h,o-z)
       parameter(kmax=8,jdmax=16,nxmax=3000,lmax=2)
       parameter(npmax=1+kmax*(jdmax+2))
       dimension theta(npmax),ti(npmax),g(npmax)
       twoeps=2.*eps
       fmin=10000.
       do 30 i=1,nreps
           do 10 j=1,npar
               ti(j)=twoeps*(rng(1)-0.5)
10         continue
           call objfun(0,npar,ti,fi,g,1)
           if(fi.lt.fmin) then
               call dcopy(npar,ti,1,theta,1)
               fmin=fi
           endif
30     continue
       return
       end

        subroutine setidmat(h,n)
	implicit real*8(a-h,o-z)
	dimension h(n,n)
	do 5 i=1,n
            do 3 j=1,n
               h(i,j)=0.
3           continue
            h(i,i)=1.
5       continue
	return
	end

      FUNCTION RNG(IXX)
      implicit double precision (a-h,o-z)
C  UNIFORM PSEUDORANDOM NUMBER GENERATOR
C  FORTRAN VERSION OF LEWIS, GOODMAN, MILLER
C  SCHRAGE,  ACM TOMS V.5 (1979) P132
C  FIRST CALL SETS SEED TO IXX, LATER IXX IGNORED
      INTEGER A,P,IX,B15,B16,XHI,XALO,LEFTLO,FHI,K
      DATA A/16807/,B15/32768/,B16/65536/,P/2147483647/
      DATA IX/0/
      IF(IX.EQ.0) IX=IXX
      XHI=IX/B16
      XALO=(IX-XHI*B16)*A
      LEFTLO=XALO/B16
      FHI=XHI*A+LEFTLO
      K=FHI/B15
      IX=(((XALO-LEFTLO*B16)-P)+(FHI-K*B15)*B16)+K
      IF(IX.LT.0) IX=IX+P
      RNG=FLOAT(IX)*4.656612875E-10
      RETURN
      END



C------ CONFUN & MAKEOPT are for use by NPSOL, otherwise not used 
        SUBROUTINE CONFUN(MODE,NCNLN,N,NROWJ,NEEDC,X,C,CJAC,NSTATE)
        IMPLICIT double precision (A-H,O-Z)
        INTEGER*4 MODE,NCNLN,N,NROWJ
        INTEGER*4 NEEDC(*)
        REAL*8 X(N),C(*),CJAC(NROWJ,*)
        SAVE
        RETURN
        END

        subroutine makeopt(bigbnd,itmax1,itmax2,ftol1,ftol2,msglvl,
     *                     mverify)
        implicit double precision(a-h,o-z)

      OPEN(7,FILE='options1.dat',status='unknown')
      rewind(7)
      WRITE(7,7000)
      write(7,7010)
      WRITE(7,7001) BIGBND
      WRITE(7,7002) ITMAX1
      WRITE(7,7004) FTOL1
      WRITE(7,7005) MSGLVL
      WRITE(7,7006) MVERIFY
      WRITE(7,7900)
      close(7)


      OPEN(7,FILE='options2.dat',status='unknown')
      rewind(7)
      WRITE(7,7000)
      write(7,7010)
      WRITE(7,7001) BIGBND
      WRITE(7,7002) ITMAX2
      WRITE(7,7004) FTOL2
      WRITE(7,7005) MSGLVL
      WRITE(7,7006) MVERIFY
      WRITE(7,7900)
      close(7)

      ITMAX3=100
      OPEN(7,FILE='options3.dat',status='unknown')
      rewind(7)
      WRITE(7,7000)
      write(7,7010)
      WRITE(7,7001) BIGBND
      WRITE(7,7002) ITMAX3
      WRITE(7,7004) FTOL2
      WRITE(7,7005) MSGLVL
      WRITE(7,7006) MVERIFY
      WRITE(7,7900)
      close(7)

7000  FORMAT('BEGIN',74X,' ')
7010  format('   NOLIST                    ',48x,' ')
7001  FORMAT('   INFINITE BOUND            ',D15.3,35X,' ')
7002  FORMAT('   MAJOR ITERATION LIMIT     ',I15  ,35X,' ')
7004  FORMAT('   OPTIMALITY TOLERANCE      ',D15.3,35X,' ')
7005  FORMAT('   MAJOR PRINT LEVEL         ',I15  ,35X,' ')
7006  FORMAT('   VERIFY LEVEL              ',I15  ,35X,' ')
7008  FORMAT('   LINESEARCH TOLERANCE      ',D15.3,35X,' ')
7900  FORMAT('END  ',74X,' ')
      return
      end


        subroutine sdev(x,n1,n,sd,ave)
c
c    Purpose: calculates standard deviation of values n1,n2,....n
c             in array x of length n
c
        implicit double precision(a-h,o-z)
        dimension x(n)
        xn=float(n-n1+1)
        ave=0.
        do 10 i=n1,n
            ave=ave+x(i)
10      continue
        ave=ave/xn
        var=0.
        do 20 i=n1,n
            var=var+(x(i)-ave)**2
20      continue
        var=var/(xn-1.)
        sd=sqrt(var)
        return
        end

        integer*4 function igcdiv(i1,i2)
        implicit integer*4(i-n)
        igcd=0
        imin=min(i1,i2)
        do 10 i=1,imin
            ir1=mod(i1,i)
            ir2=mod(i2,i)
            if (ir1.eq.0.and.ir2.eq.0) igcd=i
10      continue
        igcdiv=igcd
        return
        end

        subroutine recontp(xt,xdat,xlag,nxmax,jdmax,nx,ldelay,
     *     jd,jt1,m,jtp)
        implicit double precision(a-h,o-z)
        dimension xt(nxmax),xdat(nxmax),xlag(nxmax,jdmax)
        m=0
        do 10 i=1,nx-jt1
           m=m+1
           ipred=i+jt1
           xt(i)=xdat(ipred)
	   if (jd.gt.0) then	
             do 5 j=1,jd
               iback=jtp+(j-1)*ldelay
               xlag(i,j)=xdat(ipred-iback)
5            continue
	   endif
10      continue
        return
        end
