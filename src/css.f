      subroutine css(h,npoint,x,y,wght,sy,trace,diag,vlam,  
     +                  ngrid,xg,yg,job,ideriv,ierr)  
   
C COMPUTES A CUBIC SMOOTHING SPLINE FIT TO A SET OF DATA GIVEN  
C NPOINT(=NUMBER OF DATA VALUES) AND LAMBDA(=VALUE OF  
C THE SMOOTHING parameter). THE CODE IS ADAPTED FROM A  
C PROGRAM IN DEBOOR,C.(1978).A PRACTICAL GUIDE TO SPLINES.  
C SPRINGER-VERLAG : NEW YORK. AN O(NPOINT) ALGORITHM  
C SUGGESTED BY HUTCHINSON AND DEHOOG IS USED TO COMPUTE  
C LVRAGE VALUES AND CONSTRUCT CONFIDENCE INTERVALS.  
c Adapted from Randy Eubank Texas A&M  
c  
c  
c this subroutine solves the problem:  
c  
c   minimize  (1/n) sum(i=1,n) [ (y(i) - f(x(i)))/wght(i) ]**2 + lambda*J(f)  
c    over f  
c   The solution will always be a piecewise cubic polynomial with join   
c   points at the x values. Natural boundary conditions are assumed: at the   
c   x(1) and x(npoints) the second and third derivatives of the slotuion   
c   will be zero   
c   All matrix calculations are done in double precision   
c  
c   Arguments of evss:   
c    h : natural log of lambda ( more convenient scale whepassing a   
c                               real*4)  
c        if h is passed with a value less than or equal -1000 no smoothing will   
c        be done and the spline will interploate the data points  
c    npoint: number of observations  
c    (x,y) : pairs of data points to be smoothed  
c    sy : on return predicted values of f at x  
c    wght : weights used in sum of squares. If the y have unequal   
c      variance then an appropriate choice for wght is the standard deviation  
c      (These weights are not normalized so multiplying by a constant   
c      will imply solving the minimization problem with a different smoothing   
c      parameter)    
c   trace: in matrix from Yhat= A(lambda)Y  with A(lambda) an nxn matrix  
c          trace= tr(A(lambda)) = " effective number of paramters"  
c   diag: diagonal elements of A(lambda) ( this is the mostt   
c         computationally intetnsive opertation in this subroutine.)  
c   vlam: value of the generalized cross-validation functyion (Used to  
c         select an appropriate value for lambda based on the data.)  
c   ngrid,xg,yg : on return, the ith deriv of the spline estimate is  
c                 evaluated on a grid, xg, with ngrid points.  The  
c                 values are returned in yg  
c  
c   ideriv: 0 = evaluate the function according to the job code  
c           1 = evaluate the first derivative of the function according  
c               to the job code  
c           2 = evaluate the second derivative of the function according  
c               to the job code  
c  
c   job: in decimal job=(a,b,c)  (a=igcv, b=igrid)   
c        a=0  evaluate spline at x values, return predicted values in sy  
c        a=1  same as a=0 plus returning values of trace, diag and vlam  
c        a=2  do no smoothing, interpolate the data  
c        a=3  same as a=1 but use the passed value invlam argument as  
c             a cost an offset factors in the diag vector  
c  
c        b=0  do not evaluate the spline at any grid points  
c        b=1  evaluate the spline (ideriv=0) or the ith deriv (i=ideriv)  
c             of the spline at the (unique) sorted,data points, xg, return yg  
c        b=2  evaluate the spline (ideriv=0) or the ith deriv (i=ideriv)  
c             of the spline at ngrid equally spaced points between x(1)  
c             and x(npoints)      
c        b=3  evaluate the spline (ideriv=0) or the ith deriv (i=ideriv)  
c             of the spline at ngrid points on grid passed in xg.
c             NOTE: extrapolation of the estimate beyond the range of   
c                     the x's will just be a linear function.  
c  
c
c  
c         c=1    Assume that the X's are in sorted order
c         c=2 Do not sort the X's use the  current set of keys
c              Should only be used on a second call to smooth
c              same design
c
c  ierr : on return ierr>0 indicates all is not well.   
      parameter (NMAX=20000)  
      implicit real*8 (a-h,o-z)  
      REAL*8 h,trace,vlam  
      REAL*8 wght(npoint),X(npoint),Y(npoint),sy(npoint),diag(npoint)  
      REAL*8 xg(ngrid),yg(ngrid)  
      REAL*8 A(NMAX,4),V(NMAX,7)  
      REAL*8 P,SIXP,SIX1MP,cost  
      REAL*8 ux(NMAX),uy(NMAX), uw(NMAX),ud(NMAX),utr 
      integer imx(NMAX) 
      integer idx(NMAX)  
      integer npoint,igcv,igrid, isort,  job(3)  
      data eps / 1d-7/  
c    cost =1.0 corrsponds to GCV  
c        igcv= job/10  
c        igrid= job- 10*(job/10)  
        igcv= job(1)
        igrid=job(2)
        isort=job(3)
      cost=1.0  
      offset= 0
c  
c      write(*,*) ' evss npoint, job', npoint,job  
      if( npoint.gt.NMAX) then  
           ierr=1  
           return  
       endif 
c initialize unique vector and two pointers
      
      do 5 k=1,npoint
      ux(k)=x(k)

      idx(k)=k
    5 continue

c sort vector X along with the keys in imx
c
c       initialize the indices  imx
c
      if( isort.le.1) then 

          do 6 k=1,npoint
            imx(k)=k
    6     continue
       endif
c   
c    sort the X's permuting the indices
c 
       if(isort.eq.0) then 
          call sortm( ux, imx,npoint)
c   the rank of the value x( imx(k)) is k
       endif 
c put y and the weights in the  sorted order 
c  
C**** construct vector consisting of unique observations.  
        
       jj=1 
       ind= imx(1) 
       ux(jj)= x(ind)  
       uy(jj)= y(ind)  
       uw(jj)= wght(ind)  
       idx(1)=jj  
  
       do 10 k=2,npoint
c   we are looping through ranks but this is not how the 
c   order of the X are stored. The location of the kth smallest
c   is at idx(k) 
           kshuf= imx(k)  
          if( abs( ux(jj)-x(kshuf)).lt.eps) then  
c**** we have a repeat observation, update weight and the weighted   
c***** average at this point  
            temp1=  1.0d0/uw(jj)**2  
            temp2=  1.0d0/wght(kshuf)**2  
            temp3 = (temp1 + temp2)  
            uy(jj)=  ( uy(jj)*temp1 + y(kshuf)*temp2)/temp3  
            uw(jj)= 1.0d0/dsqrt(temp3)  
          else  
            jj= jj+1  
            ux(jj)= x(kshuf)  
            uy(jj)= y(kshuf)  
            uw(jj)= wght(kshuf)  
          endif
c  save the value that indexes unique values to repliacted ones.
c    x(k) corresponds to the unique X at idx(k)  
          idx(k)=jj  
c     write(*,  *) "devssr", ux(jj), x(kshuf),uy(jj), y(kshuf) 
   10 continue  
      nunq= jj  
  
        itp=0  
        if(igcv.eq.2) itp=1  
  
c   handle condition for interpolation        
  
      if(itp.eq.0) then   
          P=1.d0/(npoint*dexp(h)+ 1.d0)  
      else  
          P=1.d0  
      endif  

      
      call dSETUP(uX,uW,uY,nunq,V,A(1,4),NMAX,itp,ierr)  
C****  check for duplicate X's if so exit   
      if(ierr.gt.0) then  
c         write(*,*) '**ERROR in dSETUP: duplicate Xs '  
         return  
      endif  
  
      call dCHOLD(P,V,A(1,4),nunq,A(1,3),A(1,1),NMAX)  
  
c compute predicted values   
  
      SIX1MP=6.d0*(1.d0-P)  
      if(itp.eq.0) then   
         DO 61 I=1,Nunq  
            a(i,1)=uY(I) - SIX1MP*(uW(I)**2)*A(I,1)  
   61    continue  
  
c fill in predicted values taking into account repeated data  
         do 70 k=1,npoint  
            jj= idx(k)  
            sytemp= a(jj,1)
c 
c unscramble the smoothed Y's
           kshuf= imx(k)  

           sy(kshuf)=sytemp  
   70    continue  
      else  
         do 60 i=1,nunq  
            a(i,1)=uy(i)  
   60    continue  
      endif  
  
c return estimates on unique x's if igrid =1  
      if(igrid.eq.1) then  
         do 65 i=1,nunq  
            xg(i)=ux(i)  
             yg(i)=a(i,1)  
   65    continue  
         ngrid=nunq  
      endif  
c      write(*,*) " grid value ************", igrid  
      if(igrid.ge.1) then   
c  
c********* evaluate spline on grid  
C piecewise cubic COEFFICIENTS ARE STORED IN A(.,2-4).   

         SIXP=6.d0*P  
         DO 62 I=1,nunq  
   62       A(I,3)=A(I,3)*SIXP  
         NPM1=nunq - 1  
         DO 63 I=1,NPM1  
            A(I,4)=(A(I+1,3)-A(I,3))/V(I,4)  
            A(I,2)=(A(I+1,1)-A(I,1))/V(I,4)  
     *        -(A(I,3)+A(I,4)/3.*V(I,4))/2.*V(I,4)  
   63    continue  
c  
c  create equally spaced x's if asked for ( igrid=2)  
c  
         if( igrid.eq.2) then   
              step= (ux(nunq)-ux(1))/(ngrid-1)  
              xg(1)=ux(1)  
              do 190 j=2,ngrid-1  
                 xg(j)= xg(j-1)+step  
  190         continue  
              xg(ngrid)=ux(nunq)  
         endif  
  
c        if( (xg(1).lt.x(1)).or.(xg(ngrid).gt.ux(nunq))) then   
c           write(*,*) '**ERROR in devssr: grid out of data range'  
c           write(*,*) 'data: ',x(1),ux(nunq),' grid: ',xg(1),xg(ngrid)  
c           ierr=2  
c           return  
c        endif  
  
         uxmin= ux(1)  
         uxmax= ux(nunq)  
  
         a1= a(1,1)  
         an= a(nunq,1)  
  
         b1= a(1,2)  
  
         d= ux(nunq)- ux(nunq-1)  
         ind= nunq-1  
         bn= a(ind,2) + a(ind,3)*d + a(ind,4)*(d**2)/2.d0  
  
c evalute spline by finding the interval containing the evaluation   
c point and applying the cubic formula for the curve estiamte  
c finding the interval such that ux(ind) <=xg(j) < ux( ind+1)  
c is done using a bisection search   
  
         do 195 j=1,ngrid  
           lint= ifind(xg(j),ux,nunq)  
c           write(*,*) lint , " lint value"
          if( lint.eq.0) then   
          d= xg(j)-uxmin  
               
            if (ideriv .eq. 0)   
     -         yg(j)= a1 + b1*d   
            if (ideriv .eq. 1)  
     -         yg(j)=  b1  
            if (ideriv .eq. 2)  
     -         yg(j)= 0.0  
         endif  
          if( lint.eq.nunq) then   
          d= xg(j)-uxmax  
               
            if (ideriv .eq. 0)   
     -         yg(j)= an + bn*d   
            if (ideriv .eq. 1)  
     -         yg(j)=  bn  
            if (ideriv .eq. 2)  
     -         yg(j)= 0.0  
         endif  
            if( ((lint.ge.1 ) .and.( lint.lt.nunq))) then   
                ind=lint  
c            a1=a(ind,1)  
c            a2=a(ind,2)  
c            b=a(ind,3)/2.d0  
c            c=a(ind,4)/6.d0  
c  
            d= xg(j)-ux(ind)  
  
            if (ideriv .eq. 0)   
     -         yg(j)= a(ind,1) + a(ind,2)*d + a(ind,3)*(d**2)/2.d0  
     -                + a(ind,4)*(d**3)/6.d0  
            if (ideriv .eq. 1)  
     -         yg(j)= a(ind,2) + a(ind,3)*d + a(ind,4)*(d**2)/2.d0  
            if (ideriv .eq. 2)  
     -         yg(j)= a(ind,3) + a(ind,4)*d  
         endif  
c         write(*,*) " ev loop",xg(j), yg(j)
c  
  195  continue   
      endif  
c****end of evaluation block  
  
      if((igcv.eq.1).or.( igcv.eq.3)) then  
        if( igcv.eq.3)  then 
              cost= diag(1)
              offset=diag(2)
         endif
c*****                       begin computing gcv and trace  
C COMPUTE LVRAGE VALUES ,THE VARIANCE ESTIMATE   
C     SGHAT2 AND CONFIDENCE INTERVALS.  
c  
         call dLV(nunq,V,uw,SIX1MP,utr,ud,NMAX)  
  
         rss=0.d0  
         trace=0.d0  
         vlam=0.d0  
  
         do 100 i=1,nunq  
c            rss= rss + ((uy(i)-a(i,1))/uw(i))**2  
            trace= trace +ud(i)  
  100    continue  
  
         do 110 k=1,npoint 
            kshuf= imx(k)  
            jj= idx(k) 
            diag(kshuf)= ud(jj)
            rss= rss + ( (y(kshuf)- sy(kshuf))/wght(kshuf) )**2 
  110    continue  
         ctrace=  2+  cost*( trace-2)   
         if( (npoint -ctrace -offset) .gt. 0.d0) then  
            vlam= (rss/npoint)/( 1- (ctrace-offset)/npoint)**2   
         else  
            vlam=1e20  
         endif  
  
      endif  
      
       
  
      return  
  
      END  
  
      INTEGER FUNCTION IFIND(X,XK,N)                                   
C  FIND I SUCH THAT XK(I) LE X LT XK(I+1)                              
C  IFIND=0  IF X LT XK(1)                                              
C  IFIND=N  IF X GT XK(N)                                              
C  J F MONAHAN  JAN 1982  DEPT OF STAT, N C S U, RALEIGH, NC 27650     
      REAL*8 X,XK(N)                                                   
      IFIND=0                                                          
      IF(X.LT.XK(1)) RETURN                                            
      IFIND=N                                                          
      IF(X.GE.XK(N)) RETURN                                            
      IL=1                                                             
      IU=N                                                             
  1   IF(IU-IL.LE.1) GO TO 4                                           
      I=(IU+IL)/2                                                      
      IF(X-XK(I)) 2,5,3                                                
  2   IU=I                                                             
      GO TO 1                                                          
  3   IL=I                                                             
      GO TO 1                                                          
  4   IFIND=IL                                                         
      RETURN                                                           
  5   IFIND=I                                                          
      RETURN                                                           
      END                                                              
 
  
   
       SUBROUTINE SORTM(K,ki,N)  
C  HEAPSORT ALGORITHM FOR SORTING ON VECTOR OF KEYS K OF LENGTH N      
C  J F MONAHAN        TRANSCRIBED FROM KNUTH, VOL 2, PP 146-7.        
C integer array ki is permuted along with K  
   
      REAL*8 K(N),KK   
      integer ki(N),kki  
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
  
  
  
  
         



  
  
         



** finds value of h minimizing the  generalized cross-validation

      double precision function gcvfc(h,nobs,x,y,wght,c,offset,trace)

      parameter(mxM=20000)
      implicit double precision (a-h,o-z)
         
      REAL*8 h,trace
      REAL*8 x(nobs),y(nobs),wght(nobs)
      REAL*8 sy(mxM),diag(mxM),dumm1(1),dumm2(1)
      integer job(3),ideriv

      data ideriv/0/
      data rinf/1e20/
       job(1)=3
       job(2)=0
       job(3)=0
       diag(1)=c
       diag(2)=offset
      call css(h,nobs,x,y,wght,sy,trace,diag,vlam,ndum,dumm1,dumm2,
     -            job,ideriv,ierr)
      gcvfc= vlam
c
c      write(*, 100) h, c, offset, trace, vlam
 100   format( 5e12.4)
c
c
      return
      end

      subroutine gcvcss(n,x,y,wt,c,offset,nstep,maxit,hopt,vopt,
     -     tropt,mxstep,fout,
     -                 ierr)
      implicit double precision (a-h,o-z)
      real*8  fout(mxstep,3)
      real*8  x(n),y(n),wt(n)
      real*8  hopt,vopt,trace
      
      data tau,tausq/.6180339d0,.3819660d0/
**** coarse search in bandwidth h: this search feeds (hl,hr) to golden search
c compute range for coarse search   
      if (mxstep .lt. nstep) then
         ierr=10
         return
      endif

      xmin= x(1)
      xmax= x(1)
      do 500 i=1,n
          if( x(i).lt.xmin) xmin=x(i)
          if( x(i).gt.xmax) xmax=x(i)
 500  continue
      range = xmax- xmin

      xnobs=n
c**  c= cost value on GCV function there is a pole at trp=2+ nobs/c
      trp=(Xnobs-offset)*.97/c 
      if( trp.le.2.0)  then
c** offset is too big!
            ierr=11
       return
       endif
      hmin=   -0.755762d0 +  0.706693d0*log(xnobs)
     *  +  0.01884722d0*(log(xnobs)**2)  -4.918114d0*log(trp)
     *  +  0.0879931d0*(log(trp)**2)
c      write(*,*) hmin, range, xnobs 
      hmin= hmin + log(range)*3.d0-log(xnobs)

      trp=2.01d0
      hmax=   -0.755762d0 +  0.706693d0*log(xnobs)
     *  +  0.01884722d0*(log(xnobs)**2)  -4.918114d0*log(trp)
     *  +  0.0879931d0*(log(trp)**2)
 
      hmax= hmax + log(range)*3.d0-log(xnobs)

      hstep=(hmax-hmin)/(nstep-1)
      do 3 j=0,nstep-1
         h=hmin+j*hstep
         gcvh= gcvfc(h,n,x,y,wt,c,offset,trace)
         fout(j+1,1)=h
         fout(j+1,2)=trace
         fout(j+1,3)=gcvh
         if ((gcvh .lt. gcvmin) .or. (j .eq. 0)) then
            hopt=h
            best=h
            gcvmin=gcvh
            trbest=trace
         endif
c          write(*,5000) h, trace,gcvh
 5000   format(5e12.4)
c
  3   continue
c
c      write(*,*) ' Crude search of optimal h completed'
c
**** fast return if crude search minimum gcv is hmin or hmax
           hopt= best
           vopt=gcvmin
           tropt=trbest
      if( (best.eq.hmin) .or. (best.eq.h) ) then
           ierr=-1
           return
      endif


c**** start values for golden search
      hl=best-hstep
      hr=best+hstep

*** Golden section search for min gcv on interval (hl,hr)--maxit(=15) iterations
***   gcv(h) must be quasiconvex on initial interval (hl,hr).
***   On return, h=abscissa of minimum, v=gcv(h)=actual minimum achieved.
***   Interval of uncertainty is (hl,hr), hl <= hlm <= hrm <= hr.

      gcvhl=gcvfc(hl,n,x,y,wt,c,offset,trace)
      gcvhr=gcvfc(hr,n,x,y,wt,c,offset,trace)
      hlm=hl*tau+hr*tausq
      hrm=hl+hr-hlm
      gcvhlm=gcvfc(hlm,n,x,y,wt,c,offset,trchlm)
      gcvhrm=gcvfc(hrm,n,x,y,wt,c,offset,trchrm)

      do 5 it=1,maxit
         if( gcvhlm .ge. gcvhrm ) then
            if (gcvhl .lt. gcvhlm) then
               err= gcvhl/gcvhlm
c              write(*,10) it,err,hl,hlm
               ierr=2
               return
            endif
c
            hl=hlm
            gcvhl=gcvhlm
            hlm=hrm
            hrm=hrm+(hrm-hl)*tau
            gcvhlm=gcvhrm
            gcvhrm=gcvfc(hrm,n,x,y,wt,c,offset,trchrm) 
         else
            if (gcvhr .lt. gcvhrm) then
               err= gcvhrm/gcvhr
c              write(*,11) it,err,hr,hrm
               ierr=2
               return
            endif
c
            hr=hrm
            gcvhr=gcvhrm
            hrm=hlm
            hlm=hlm+(hlm-hr)*tau
            gcvhrm=gcvhlm
            gcvhlm=gcvfc(hlm,n,x,y, wt,c,offset,trchlm)
         endif
5     continue

**** finished -- take the best h so far
      if( gcvhlm .ge. gcvhrm) then
         hopt=hrm
         vopt=gcvhrm
         tropt=trchrm
      else
         hopt=hlm
         vopt=gcvhlm
         tropt=trchlm
      endif
c
10    format(' gcv(h) is NOT quasiconvex ',/,
     - ' iter,error,hl,hlm: ',I4,e15.5,2e15.5)
11    format(' gcv(h) is NOT quasiconvex ',/,
     - ' iter,error,hr,hrm: ',I4,e15.5,2e15.5)
c
c
c     write(*,*) 'hopt', hopt
      return
      end

      SUBROUTINE dSETUP(X,WGHT,Y,NPOINT,V,QTY,NMAX,itp,ierr)
C   PUT DELX=X(.+1)-X(.) INTO V(.,4)
C   PUT THE THREE BANDS OF THE MATRIX Q-TRANSP*D INTO
C     V(.,1-3)
C   PUT THE THREE BANDS OF (D*Q)-TRANSP*(D*Q) AT AND 
C     ABOVE THE DIAGONAL INTO V(.,5-7)
C   HERE Q IS THE TRIDIAGONAL MATRIX OF ORDER (NPOINT
C     -2,NPOINT) THAT SATISFIES Q-TRANSP*T=0 AND WGHT
C     IS THE DIAGONAL MATRIX WHOSE DIAGONAL ENTRIES 
C     ARE THE SQUARE ROOTS OF THE WEIGHTS USED IN THE 
C     PENALIZED LEAST-SQUARES CRITERION
c
      implicit double precision (a-h,o-z)
      REAL*8 WGHT(NMAX),X(NMAX),y(NMAX)
      REAL*8 QTY(NMAX),V(NMAX,7)
      REAL*8 DIFF,PREV
      INTEGER NPOINT,I,NPM1
c
      NPM1=NPOINT -1
      V(1,4)=X(2)-X(1)
      if(v(1,4).eq.0.d0) then
                ierr=5
                return
      endif

      DO 11 I=2,NPM1
         V(I,4)=X(I+1) - X(I)
         if(v(I,4).eq.0.d0) then
                ierr=5
                return
         endif
         if(itp.eq.0) then 
            V(I,1)=WGHT(I-1)/V(I-1,4)
            V(I,2)=-WGHT(I)/V(I,4) - WGHT(I)/V(I-1,4)
            V(I,3)=WGHT(I+1)/V(I,4)
         else
            V(I,1)=1.d0/V(I-1,4)
            V(I,2)=-1.d0/V(I,4) - 1.0/V(I-1,4)
            V(I,3)=1.d0/V(I,4)
         endif
   11 continue
c
      V(NPOINT,1)=0.d0
      DO 12 I=2,NPM1
   12    V(I,5)=V(I,1)**2 + V(I,2)**2 + V(I,3)**2
      IF(NPM1 .LT. 3)GO TO 14
      DO 13 I=3,NPM1
   13    V(I-1,6)=V(I-1,2)*V(I,1) + V(I-1,3)*V(I,2)
   14 V(NPM1,6)=0.d0
      IF(NPM1 .LT. 4)GO TO 16
      DO 15 I=4,NPM1
   15    V(I-2,7)=V(I-2,3)*V(I,1)
   16 V(NPM1-1,7)=0.d0
      V(NPM1,7)=0.d0
c
C  CONSTRUCT Q-TRANSP. * Y IN QTY
      PREV=(Y(2) - Y(1))/V(1,4)
      DO 21 I=2,NPM1
         DIFF=(Y(I+1)-Y(I))/V(I,4)
         QTY(I)=DIFF - PREV
   21    PREV=DIFF
c
      RETURN
      END




      SUBROUTINE dCHOLD(P,V,QTY,NPOINT,U,QU,NMAX)
c CONSTRUCTS THE UPPER THREE DIAGONALS IN V(I,J),I=2,
C   NPOINT-1,J=1,3 OF THE MATRIX 6*(1-P)*Q-TRANSP*
C   (D**2)*Q + P*R . (HERE R IS THE MATRIX PROPORTIONAL
C   TO Q-TRANSP*KSUBN*Q , WHERE KSUBN IS THE MATRIX 
C   WITH ELEMENTS K(X(I),X(J)) AND K IS THE USUAL
C   KERNEL ASSOCIATED WITH CUBIC SPLINE SMOOTHING.)
C   THE CHOLESKY DECOMPOSITION OF THIS MATRIX IS COMPUTED
C   AND STORED IN V(.,1-3) AND THE EQUATION SYSTEM FOR 
C   THE QUADRATIC COEFFICIENTS OF THE SPLINE IN ITS 
C   PIECEWISE POLYNOMIAL REPRESENTATION IS SOLVED . THE 
C   SOLUTION IS RETURNED IN U.
c
      REAL*8 P,QTY(NMAX),QU(NMAX),U(NMAX),V(NMAX,7)
      REAL*8 SIX1MP,TWOP,RATIO,PREV
      INTEGER NPOINT,I,NPM1,NPM2
c
      NPM1=NPOINT - 1
C****  CONSTRUCT 6*(1-P)*Q-TRANSP.*(D**2)*Q + P*R
      SIX1MP=6.d0*(1.d0 - P)
      TWOP=2.d0*P
      DO 2 I=2,NPM1
         V(I,1)=SIX1MP*V(I,5) + TWOP*(V(I-1,4)+V(I,4))
         V(I,2)=SIX1MP*V(I,6) + P*V(I,4)
         V(I,3)=SIX1MP*V(I,7)
   2  continue
      NPM2=NPOINT - 2
      IF(NPM2 .GE. 2)GO TO 10
      U(1)=0.d0
      U(2)=QTY(2)/V(2,1)
      U(3)=0.d0
      GO TO 41
C  FACTORIZATION : FACTORIZE THE MATRIX AS L*B-INV*
C     L-TRANSP , WHERE B IS A DIAGONAL MATRIX AND L
C     IS UPPER TRIANGULAR.
   10 DO 20 I=2,NPM2
         RATIO=V(I,2)/V(I,1)
         V(I+1,1)=V(I+1,1) - RATIO*V(I,2)
         V(I+1,2)=V(I+1,2) - RATIO*V(I,3)
         V(I,2)=RATIO
         RATIO=V(I,3)/V(I,1)
         V(I+2,1)=V(I+2,1) - RATIO*V(I,3)
   20    V(I,3)=RATIO
C  FORWARD SUBSTITUTION
      U(1)=0.d0
      V(1,3)=0.d0
      U(2)=QTY(2)
      DO 30 I=2,NPM2
   30    U(I+1)=QTY(I+1) - V(I,2)*U(I) -V(I-1,3)*U(I-1)
C  BACK SUBSTITUTION
      U(NPOINT)=0.d0
      U(NPM1)=U(NPM1)/V(NPM1,1)
      DO 40 I=NPM2,2,-1
   40    U(I)=U(I)/V(I,1) - U(I+1)*V(I,2) - U(I+2)*V(I,3)
C  CONSTRUCT Q*U
   41 PREV=0.d0
      DO 50 I=2,NPOINT
         QU(I)=(U(I)-U(I-1))/V(I-1,4)
         QU(I-1)=QU(I) - PREV
   50    PREV=QU(I)
      QU(NPOINT)=-QU(NPOINT)
c
      RETURN
      END




      SUBROUTINE dLV(NPOINT,V,WGHT,SIX1MP,TR,LEV,NMAX)
c CONSTRUCTS THE UPPER THREE DIAGONALS OF (6*(1-P)*
C   Q-TRANSP*(D**2)*Q + P*R)-INV USING THE RECURSION
C   FORMULA IN HUTCHINSON,M.F. AND DEHOOG,F.R.(1985).
C   NUMER. MATH. 47,99-106, AND STORES THEM IN V(.,5-7).
C   THESE ARE USED IN POST AND PRE MULTIPLICATION BY
C   Q-TRANSP AND Q TO OBTAIN THE DIAGONAL ELEMENTS OF
C   THE HAT MATRIX WHICH ARE STORED IN THE VECTOR LEV.
C     THE TRACE OF THE HAT MATRIX IS RETURNED IN TR.
c
      REAL*8 V(NMAX,7),TR,W1,W2,W3,SIX1MP
      REAL*8 wght(NMAX)
      REAL*8 LEV(npoint)
      INTEGER NPM1,NPM2,NPM3,NPOINT
c
      NPM1=NPOINT - 1
      NPM2=NPOINT - 2
      NPM3=NPOINT - 3
C   RECURSION FOR DIAGONALS OF INVERSE MATRIX
      V(NPM1,5)=1/V(NPM1,1)
      V(NPM2,6)=-V(NPM2,2)*V(NPM1,5)
      V(NPM2,5)=(1/V(NPM2,1)) - V(NPM2,6)*V(NPM2,2)
      DO 10 I=NPM3,2,-1
         V(I,7)=-V(I,2)*V(I+1,6) - V(I,3)*V(I+2,5)
         V(I,6)=-V(I,2)*V(I+1,5) - V(I,3)*V(I+1,6)
         V(I,5)=(1/V(I,1))- V(I,2)*V(I,6) - V(I,3)*V(I,7)
   10 CONTINUE
C   POSTMULTIPLY BY (D**2)*Q-TRANSP AND PREMULTIPLY BY Q TO
C     OBTAIN DIAGONALS OF MATRIX PROPORTIONAL TO THE 
C     IDENTITY MINUS THE HAT MATRIX.
      W1=1.d0/V(1,4)
      W2= -1.d0/V(2,4) - 1.d0/V(1,4)
      W3=1.d0/V(2,4)
      V(1,1)=V(2,5)*W1
      V(2,1)=W2*V(2,5) + W3*V(2,6)
      V(2,2)=W2*V(2,6) + W3*V(3,5)
      LEV(1)=1.d0 - (WGHT(1)**2)*SIX1MP*W1*V(1,1)
      LEV(2)=1.d0 - (WGHT(2)**2)*SIX1MP*(W2*V(2,1) + W3*V(2,2))
      TR=LEV(1) + LEV(2)
      DO 20 I=4,NPM1
         W1=1.d0/V(I-2,4)
         W2= -1.d0/V(I-1,4) - 1.d0/V(I-2,4)
         W3=1.d0/V(I-1,4)
         V(I-1,1)=V(I-2,5)*W1 + V(I-2,6)*W2 + V(I-2,7)*W3
         V(I-1,2)=V(I-2,6)*W1 + V(I-1,5)*W2 + V(I-1,6)*W3
         V(I-1,3)=V(I-2,7)*W1 + V(I-1,6)*W2 + V(I,5)*W3
         LEV(I-1)=1.d0 - (WGHT(I-1)**2)*SIX1MP*(W1*V(I-1,1) 
     .             + W2*V(I-1,2) + W3*V(I-1,3))
         TR= TR + LEV(I-1)
   20 CONTINUE
      W1=1.d0/V(NPM2,4)
      W2= -1.d0/V(NPM1,4) - 1.d0/V(NPM2,4)
      W3=1.d0/V(NPM1,4)
      V(NPOINT,1)=V(NPM1,5)*W3
      V(NPM1,1)=V(NPM2,5)*W1 + V(NPM2,6)*W2
      V(NPM1,2)=V(NPM2,6)*W1 + V(NPM1,5)*W2
      LEV(NPM1)=1.d0 - (WGHT(NPM1)**2)*SIX1MP*(W1*V(NPM1,1)
     .             + W2*V(NPM1,2))
      LEV(NPOINT)=1.d0 - (WGHT(NPOINT)**2)*SIX1MP*W3*V(NPOINT,1)
      TR= TR + LEV(NPM1) + LEV(NPOINT)
      RETURN 
      END
