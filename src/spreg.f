      subroutine devssr(h,npoint,x,y,wght,sy,trace,diag,vlam,
     +                  ngrid,xg,yg,job,ideriv,ierr)
 
C COMPUTES A CUBIC SMOOTHING SPLINE FIT TO A SET OF DATA GIVEN
C NPOINT(=NUMBER OF DATA VALUES) AND LAMBDA(=VALUE OF
C THE SMOOTHING PARAMETER). THE CODE IS ADAPTED FROM A
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
c    h : natural log of lambda ( more convenient scale when passing a 
c                               real*4)
c        if h is passed with a value less than or equal -1000 no smoothing will 
c        be done and the spline will interploate the data points
c    npoint: number of observations
c    (x,y) : pairs of data points to be smoothed
c            x(i) are assumed to be increasing. Repeated 
c            observations at the same x are not handled by this routine.
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
c   job: in decimal job=ab  (a=igcv, b=igrid) 
c        a=0  evaluate spline at x values, return predicted values in sy
c        a=1  same as a=0 plus returning values of trace, diag and vlam
c        a=2  do no smoothing, interpolate the data
c        a=3  same as a=1 but use the passed value invlam argument as
c             a cost factor in the GCV function
c
c        b=0  do not evaluate the spline
c        b=1  evaluate the spline (ideriv=0) or the ith deriv (i=ideriv)
c             of the spline at the (unique) data points, xg, return yg
c        b=2  evaluate the spline (ideriv=0) or the ith deriv (i=ideriv)
c             of the spline at ngrid equally spaced points between x(1)
c             and x(npoints)    
c        b=3  evaluate the spline (ideriv=0) or the ith deriv (i=ideriv)
c             of the spline at ngrid points on grid passed in xg. NOTE1:
c             ** xg MUST BE SORTED! xg(1)>= x(1) and xg(ngrid)<=x(npoint))
c             NOTE 2: extrapolation of the estimate beyond the range of 
c                     the x's will just be a linear function.
c
c  ierr : on return ierr>0 indicates all is not well. 
c
 
      PARAMETER (mxM=20000,NMAX=mxM,mxfit=20000)
      implicit real*8 (a-h,o-z)
      REAL*8 h,trace,vlam
      REAL*8 wght(npoint),X(npoint),Y(npoint),sy(npoint),diag(npoint)
      REAL*8 xg(ngrid),yg(ngrid)
      REAL*8 A(NMAX,4),V(NMAX,7)
      REAL*8 P,SIXP,SIX1MP,cost
      REAL*8 ux(NMAX),uy(NMAX), uw(NMAX),ud(NMAX),utr
      integer idx(NMAX)
      integer npoint
      data eps / 1d-14/
c    cost =1.0 corrsponds to GCV
      cost=1.0
      iprint=0
      if (iprint .gt. 4) then
         write(*,*) '>>** sub devssr  input'
         write(*,*) 'h,npoint,trace,vlam,ngrid,job,ideriv,ierr'
         write(*,*) h,npoint,trace,vlam,ngrid,job,ideriv,ierr
         write(*,*) 'x: ',(x(kk),kk=1,npoint)
         write(*,*) 'y: ',(y(kk),kk=1,npoint)
         write(*,*) 'wght: ',(wght(kk),kk=1,npoint)
         write(*,*) 'xg: ',(xg(kk),kk=1,50)
      endif
c
c      write(*,*) ' evss npoint, job', npoint,job
      if( npoint.gt.NMAX) then
           ierr=1
           return
       endif
c
C**** construct vector consisting of unique observations. 
       jj=1
       ux(jj)= x(1)
       uy(jj)= y(1)
       uw(jj)= wght(1)
       idx(1)=jj

       do 10 k=2,npoint
          if( abs( ux(jj)-x(k)).lt.eps) then
c**** we have a repeat observation, update weight and the weighted 
c***** average at this point
            temp1=  1.0d0/uw(jj)**2
            temp2=  1.0d0/wght(k)**2
            temp3 = (temp1 + temp2)
            uy(jj)=  ( uy(jj)*temp1 + y(k)*temp2)/temp3
            uw(jj)= 1.0d0/dsqrt(temp3)
          else
            jj= jj+1
            ux(jj)= x(k)
            uy(jj)= y(k)
            uw(jj)= wght(k)
          endif
          idx(k)=jj
   10 continue
      nunq= jj

        igcv= job/10
        igrid= job- 10*(job/10)
        itp=0
        if(igcv.eq.2) itp=1

c   handle condition for interpolation      

      if(itp.eq.0) then 
          P=1.d0/(npoint*dexp(h)+ 1.d0)
      else
          P=1.d0
      endif
c      write(*,*) 'evss igcv igrid p', igcv,igrid,p

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
            sy(k)= a(jj,1)
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
c            yg(i)=a(i,1)
   65    continue
         ngrid=nunq
      endif

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

c**** first evaluate spline beyond the range of X values
         uxmin= ux(1)
         uxmax= ux(nunq)

         a1= a(1,1)
         an= a(nunq,1)

         b1= a(1,2)

         d= ux(nunq)- ux(nunq-1)
         ind= nunq-1
         bn= a(ind,2) + a(ind,3)*d + a(ind,4)*(d**2)/2.d0

         do 195 j=1,ngrid
          d= xg(j)-uxmin
         if( d.lt.0) then 
             
            if (ideriv .eq. 0) 
     -         yg(j)= a1 + b1*d 
            if (ideriv .eq. 1)
     -         yg(j)=  b1
            if (ideriv .eq. 2)
     -         yg(j)= 0.0
         endif
          d= xg(j)-uxmax
         if( d.gt.0) then 
             
            if (ideriv .eq. 0) 
     -         yg(j)= an + bn*d 
            if (ideriv .eq. 1)
     -         yg(j)=  bn
            if (ideriv .eq. 2)
     -         yg(j)= 0.0
         endif
  195    continue

         ind=1

         do 200 j=1,ngrid

c**** ignore grid points outside range
         if( xg(j).gt.uxmax) goto  300
         if( xg(j).lt.uxmin) goto 200

  210       continue
            if(ux(ind+1).lt.xg(j)) then
               ind= ind+1
               goto 210
            endif

c            a1=a(ind,1)
c            a2=a(ind,2)
c            b=a(ind,3)/2.d0
c            c=a(ind,4)/6.d0
c
            d= xg(j)-ux(ind)
c
c            write(*,*) 'err chk in devssr: ideriv=',ideriv
            if (ideriv .eq. 0) 
     -         yg(j)= a(ind,1) + a(ind,2)*d + a(ind,3)*(d**2)/2.d0
     -                + a(ind,4)*(d**3)/6.d0
            if (ideriv .eq. 1)
     -         yg(j)= a(ind,2) + a(ind,3)*d + a(ind,4)*(d**2)/2.d0
            if (ideriv .eq. 2)
     -         yg(j)= a(ind,3) + a(ind,4)*d
c
  200    continue
c**** branch to outside loop
  300    continue

      endif
c****end of evaluation block

      if((igcv.eq.1).or.( igcv.eq.3)) then
        if( igcv.eq.3) cost=vlam
c*****                       begin computing gcv and trace
C COMPUTE LVRAGE VALUES ,THE VARIANCE ESTIMATE 
C     SGHAT2 AND CONFIDENCE INTERVALS.
c
         call dLV(nunq,V,uw,SIX1MP,utr,ud,NMAX)

         rss=0.d0
         trace=0.d0
         vlam=0.d0

         do 100 i=1,nunq
            rss= rss + ((uy(i)-a(i,1))/uw(i))**2
            trace= trace +ud(i)
  100    continue

         do 110 k=1,npoint
            diag(k)= ud(idx(k))
  110    continue
         ctrace=  2+  cost*( trace-2)
         if( nunq-ctrace .gt. 0.d0) then
            vlam= rss*nunq/( nunq- ctrace)**2 
         else
            vlam=1e20
         endif

      endif
c************************ end computing diag and gcv
      if (iprint .gt. 4) then
         write(*,*) '>>** sub devssr  final results'
         write(*,*) 'h,npoint,trace,vlam,ngrid,job,ideriv,ierr'
         write(*,*) h,npoint,trace,vlam,ngrid,job,ideriv,ierr
         write(*,*) 'x: ',(x(kk),kk=1,npoint)
         write(*,*) 'y: ',(y(kk),kk=1,npoint)
         write(*,*) 'sy: ',(sy(kk),kk=1,npoint)
         write(*,*) 'diag: ',(diag(kk),kk=1,npoint)
         write(*,*) 'wght: ',(wght(kk),kk=1,npoint)
         write(*,*) 'xg: ',(xg(kk),kk=1,50)
         write(*,*) 'yg: ',(yg(kk),kk=1,50)
      endif
c

      return

      END


      subroutine hgcvc(n,x,y,wt,c,nstep,maxit,hopt,vopt,
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
         ierr=1
         return
      endif

      range= x(n)- x(1)

      xnobs=n
c**  c= cost value on GCV function there is a pole at trp=2+ nobs/c
      trp=Xnobs*.95d0/c
      hmin=   -0.755762d0 +  0.706693d0*log(xnobs)
     *  +  0.01884722d0*(log(xnobs)**2)  -4.918114d0*log(trp)
     *  +  0.0879931d0*(log(trp)**2)
 
      hmin= hmin + log(range)*3.d0-log(xnobs)

      trp=2.01d0
      hmax=   -0.755762d0 +  0.706693d0*log(xnobs)
     *  +  0.01884722d0*(log(xnobs)**2)  -4.918114d0*log(trp)
     *  +  0.0879931d0*(log(trp)**2)
 
      hmax= hmax + log(range)*3.d0-log(xnobs)

      hstep=(hmax-hmin)/(nstep-1)

      do 3 j=0,nstep-1
         h=hmin+j*hstep
         gcvh= gcvfc(h,n,x,y,wt,c,trace)
         fout(j+1,1)=h
         fout(j+1,2)=trace
         fout(j+1,3)=gcvh
         if ((gcvh .lt. gcvmin) .or. (j .eq. 0)) then
            hopt=h
            best=h
            gcvmin=gcvh
            trbest=trace
         endif
c         write(*,5000) h, trace,gcvh
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

      gcvhl=gcvfc(hl,n,x,y,wt,c,trace)
      gcvhr=gcvfc(hr,n,x,y,wt,c,trace)
      hlm=hl*tau+hr*tausq
      hrm=hl+hr-hlm
      gcvhlm=gcvfc(hlm,n,x,y,wt,c,trchlm)
      gcvhrm=gcvfc(hrm,n,x,y,wt,c,trchrm)

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
            gcvhrm=gcvfc(hrm,n,x,y,wt,c,trchrm) 
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
            gcvhlm=gcvfc(hlm,n,x,y, wt,c,trchlm)
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






