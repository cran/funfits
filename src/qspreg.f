           



c** finds value of h minimizing the  generalized cross-validation

      double precision function cvrf
     * (h,nobs,x,y,wt,sy,trace,diag,din,dout,ierr)  


      implicit double precision (a-h,o-z)    
      REAL*8 h,trace
      REAL*8 x(nobs),y(nobs),wt(nobs)
      REAL*8 sy(nobs),diag(nobs),dumm1(1),dumm2(1)
      real*8 din(10), dout(10)
      integer ierr, job(3),ideriv,ngrid
      data job/3,0,0/
      data ideriv,ngrid/0,1/
c      write(*,*) " call to rcss", h, nobs, ierr
      call rcss(h,nobs,x,y,wt,sy,trace,diag,cv,  
     +          ngrid,dumm1,dumm2,job,ideriv,din,dout,ierr)  
      nit= int( dout(1))
c      write(*, 100) h,  trace, cv, nit,ierr
 100   format( 3e12.4, i8,i5)
c
c
      trace=dout(3)
       cvrf= cv
      return
      end

       subroutine cvrcss(
c    arguments passed through to rcss
     + n,x,y,wt,sy,diag,din,dout,  
c
c  arguments needed for cv search
     +  nstep,maxit,hmin,hmax,hopt,vopt,
     +     tropt,mxstep,tabout,
     +                 ierr)
      implicit double precision (a-h,o-z)
      real*8  tabout(mxstep,4)
      real*8  x(n),y(n),wt(n),sy(n), diag(n), din(10), dout(10)
      real*8  hopt,vopt,tropt,hmin, hmax
      data tau,tausq/.6180339d0,.3819660d0/
**** coarse search in bandwidth h: this search feeds (hl,hr) to golden search
c compute range for coarse search   
      if (mxstep .lt. nstep) then
         ierr=10
         return
      endif
      hstep=(hmax-hmin)/(nstep-1)
      do 3 j=0,nstep-1
         h=hmax - j*hstep
         cvh= cvrf(h,n,x,y,wt,sy,trace,diag,din,dout,ierr)  


         tabout(j+1,1)=h
         tabout(j+1,2)=dout(3)
         tabout(j+1,3)=cvh
         tabout(j+1,4)=dout(1)
         if ((cvh .lt. cvmin) .or. (j .eq. 0)) then
            hopt=h
            best=h
            cvmin=cvh
            trbest=trace
         endif
c          write(*,5000) h, trace,cvh
 5000   format(5e12.4)
c
  3   continue
c
c      write(*,*) ' Crude search of optimal h completed'
c      write(*,*) 'best at : ',best, tropt, cvmin
c**** fast return if crude search minimum cv is hmin or hmax
           hopt= best
           vopt=cvmin
           tropt=trbest
      if( (best.le.hmin) .or. (best.ge.hmax) ) then
           ierr=-1
           return
      endif


c**** start values for golden search
      hl=best-hstep
      hr=best+hstep

*** Golden section search for min cv on interval (hl,hr)--maxit iterations
***   cv(h) must be quasiconvex on initial interval (hl,hr).
***   On return, h=abscissa of minimum, v=cv(h)=actual minimum achieved.
***   Interval of uncertainty is (hl,hr), hl <= hlm <= hrm <= hr.

      cvhl=cvrf(hl,n,x,y,wt,sy,trace,diag,din,dout,ierr)  

      cvhr=cvrf(hr,n,x,y,wt,sy,trace,diag,din,dout,ierr)  

      hlm=hl*tau+hr*tausq
      hrm=hl+hr-hlm
      cvhlm=cvrf(hlm,n,x,y,wt,sy,trchlm,diag,din,dout,ierr)  

      cvhrm=cvrf(hrm,n,x,y,wt,sy,trchrm,diag,din,dout,ierr)  


      do 5 it=1,maxit
         if( cvhlm .ge. cvhrm ) then
            if (cvhl .lt. cvhlm) then
               err= cvhl/cvhlm
c              write(*,10) it,err,hl,hlm
               ierr=2
               return
            endif
c
            hl=hlm
            cvhl=cvhlm
            hlm=hrm
            hrm=hrm+(hrm-hl)*tau
            cvhlm=cvhrm
            cvhrm=cvrf(hrm,n,x,y,wt,sy,trchrm,diag,din,dout,ierr)  

         else
            if (cvhr .lt. cvhrm) then
               err= cvhrm/cvhr
c              write(*,11) it,err,hr,hrm
               ierr=2
               return
            endif
c
            hr=hrm
            cvhr=cvhrm
            hrm=hlm
            hlm=hlm+(hlm-hr)*tau
            cvhrm=cvhlm
            cvhlm=cvrf(hlm,n,x,y,wt,sy,trchlm,diag,din,dout,ierr)  
            
         endif
5     continue

**** finished -- take the best h so far
      if( cvhlm .ge. cvhrm) then
         hopt=hrm
         vopt=cvhrm
         tropt=trchrm
      else
         hopt=hlm
         vopt=cvhlm
         tropt=trchlm
      endif
c
10    format(' cv(h) is NOT quasiconvex ',/,
     - ' iter,error,hl,hlm: ',I4,e15.5,2e15.5)
11    format(' cv(h) is NOT quasiconvex ',/,
     - ' iter,error,hr,hrm: ',I4,e15.5,2e15.5)
c
c
c     write(*,*) 'hopt', hopt
      return
      end






       subroutine rcss(h,npoint,x,y,wt,sy,trace,diag,cv,  
     +                  ngrid,xg,yg,job,ideriv,din,dout,ierr)  
c This a program to compute a robust univariate spline according to the
c model:
c   minimize   (1/n) sum(i=1,n)[rho( y(i)-f(x(i) )] + lambda*J(f)
c    over f
c definition of the rho function and its derivative are in rcsswt
c and rcssr
c
c One way of speeding convergence is to use the results from a previous
c estimate as the starting values for the another estimate. This is
c particulary appropriate when lambda or a parameter in the rho function
c is being varied. Moderate changes in lambda will often yield similar
c estimates. The way to take advantage of this is to pass the weights
c from the previous fit as teh starting values for the next estimate
c   
c   Arguments of rcss:   
c    h : natural log of lambda  
c                               
c        if h is passed with a value less than or equal -1000 no smoothing will   
c        be done and the spline will interploate the data points  
c    npoint: number of observations  
c    (x,y) : pairs of data points to be smoothed  
c            x(i) are assumed to be increasing. Repeated   
c            observations at the same x are not handled by this routine.  
c    sy : on return predicted values of f at x  
c    wt : weights used in the iterivatively reweighted least 
c           squares algorithm to compute robust spline. Vector
c           passed are used as the starting values. Subsequent
c           iterations compute the weights by a call to  the 
c           subroutine     rcsswt  
c
c          that is the linear approximation of teh estimator at 
c              convergence. 
c          trace= tr(A(lambda)) = " effective number of paramters"  
c   diag: diagonal elements of A(lambda) ( this is the most   
c         computationally intetnsive opertation in this subroutine.)  
c   cv: approximate cross-validation function 
c         using the linear approximation at convergence
c 
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
c   din:  Vector of input parameters 
c           din(1)= cost for cv
c           din(2)= offset for cv
c           din(3)= max number of iterations
c           din(4)= tolerance criterion for convergence
c
c           din(5)= C scale parameter in robust function (transition from
c                   quadratic to linear.
c           din(6)= alpha   1/2 slope of rho function for large, positive
c                   residuals   slope for residuals <0 : is  1/2 (1-alpha) 
c                  see comments in rcsswt for defintion of the  rho and psi 
c                  functions
c
c   job: in decimal job=(a,b,c)  (a=igcv, b=igrid)   
c        a=0  evaluate spline at x values, return predicted values in sy  
c        a=1  same as a=0 plus returning values of trace, diag and cv  
c        a=2  do no smoothing, interpolate the data  
c        a=3  same as a=1 but use the passed values in din(1) din(2)
c             for computing cv function
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
 
c Arguments of subroutine:
c    dout(1)= numit
c    dout(2)=tstop
c    dout(3) = trace
c    dout(4)= cv
c      numit: actual number of iterations for convergence
c      tstop: value of convergence criterion at finish.
c
c ierr: if ierr>0 something bad has happened
c       ierr>10 indicates problems in the call to the cubic spline
c       routine.
c  
   
      parameter(NMAX=20000)  
      implicit double precision (a-h,o-z)  
      REAL*8 h,trace,cv  
      REAL*8 wt(npoint),X(npoint),Y(npoint),sy(npoint),diag(npoint)  
      REAL*8 xg(ngrid),yg(ngrid)  
      REAL*8 din(10), dout(10),cost,  offset, dum1, dum2
      
      integer npoint,ngrid ,itj(3), job(3)
 
      if(npoint.gt.NMAX) then
c            write(*,*) 'nobs= ',npoint,' too large for rspl'
           ierr=1
           return
      endif

c Set up initial estimates for first iteration
c if sav(1)=1 uses the previous weights (even if these
c came from another call)
c

       maxit= int(din(3))
       tstop=din(4)
       ybar=0.0
       ysd=0.0       
       do 10 k=1,npoint
            diag(k) = y(k)
            ybar= ybar+ y(k)
            ysd= ysd+ y(k)**2
   10  continue 
        ybar= ybar/npoint
        ysd= sqrt( ysd/npoint - ybar**2)     

c Start iterating
          test=0.0

       itj(1)= 0
       itj(2)=0
       itj(3)=0
      do 500 it=1,maxit
c      write(*,*) 'iteration', it
      if( it.gt.1) then
        itj(3)=2
      endif
c fit a weighted least squares cubic smoothing spline
      call css(h,npoint,x,y,wt,sy,
     *  dum1,diag,dum2,ngrid,xg,yg,
     *  itj,ideriv,ierr)

c check for an error returned by spline routine
      if(ierr.gt.0) then
c         add 10 so these can be distinguished from errors in this routine
          ierr= 10 + ierr
c          write(*,*) " error from css ----------------------------"
          return
       endif

c compute convergence criterion
c The intent of this criterion is to find out when successive iterations
c produce changes that are small 
c     write(*,*) 'compute convergence criterion'


          do 25 i=1,npoint
               test=(diag(i)-sy(i))**2 + test 
               diag(i)= sy(i)
 25       continue

           test=sqrt(test/npoint)/ysd
c           write(*,*) "test value, stop value: ", test, stop
           if( test.lt.tstop  ) then
c            * exit loop *
               numit=it
               goto 1000
            endif
             
c
c    make up new set of weights
c
          call rcsswt( npoint, y,sy,wt, din(5))
c          write(*,*)"new weights"
 480       format( 5e15.6)
c           write(*,480) ( wt(k), k=1,npoint)
c    reinitialize test criterion for convergence
c
      test=0.0
500   continue

      numit= maxit
 1000 continue
c One last call if job code is not 0 
      if( (job(1).ne.0).or.(job(2).ne.0)) then
c
      call css(h,npoint,x,y,wt,sy,
     *  trace,diag,cv,ngrid,xg,yg,
     *  job,ideriv,ierr)

      ia= job(1)
      ig= job(2)
c      if(ig.gt.0)  then
c          write(*,*) " info from grid eval"
c          write(*,*) job, ngrid , xg(1), yg(1), xg(10), yg(10)
c      endif
c calculate cv value if asked for 
      if( (ia.eq.1) .or.( ia.eq.3) ) then 
             if(ia.eq.3) then 
             cost= din(1)
             offset= din(2)/npoint
             else
             cost=1
             offset= 0
             endif

             cv=0.0
             do 1500 k=1,npoint
c compute approx. cross-validated residual
                  resid= (y(k)- sy(k))/( 1- cost*(diag(k)+offset))

c plug cv residual into rho function, din(5) is the begining of parameter 
c vector for the rho function (scale and alpha)
                  cv= cv + rcssr(resid, din(5))    
 1500        continue
             cv= cv/npoint
      endif 
      endif

      dout(1)=numit
      dout(2)=test
      dout(3)=trace
      dout(4)=cv

      return
      end
 
      double precision function rcssr(r,par)
c
c     robust rho function:
c  This is a peicewise polynomial with knots at -C , 0 and C
c  the function is quadratic for -C<u<0 and 0<u<C
c  the function is linear for u<-C and u>C
c  rho is continuous for all u aqnd differentiable for all points 
c  except u=0 when a != 1/2 
c   
c
c    rho(u) =      2*a*u - a*c     for u>C
c                  a*u**2/C        for   0<u< C   
c                  (1-a)*u**2/C    for -C<u<0
c                  2*(1-a)*u - (1-a)*C  for u< -C
c
c        Note a= par(1), C= par(2)
      implicit real*8 (a-h, o-z)
      real*8 r, par(2),c,a
      c= par(1)     
      if( r.gt.0 ) then 
         a=par(2)
       else
         a =(1-par(2))
         r= -r
      endif 
      if( r.le.c) then
            rcssr= a*r*r/c
      else
           rcssr= 2*a*(r) - a*c
      endif
      return
      end
c**********
      subroutine rcsswt(n,y, sy, wt, par)
      implicit double precision (a-h, o-z)
      real*8 y(n), sy(n), wt(n),psi,a,am1,c
      real*8 par(2)
c
c   psi(u) is the derivative of rho(u) defined in rcssr above
c
c   It is composed of peicewise linear and peicewise constant segements
c   and will be continuous except at u for a!=.5. 
c
c
        a= par(2)
        am1 = (1-par(2))
        c= par(1)
        do 100 k=1, n
c   find scaled residual
               r= (y(k)- sy(k))/c
               if( (r.gt. 0)) then 
                         if( r.lt. 1) then 
                               psi= 2*a*r
                               
                         else
                              psi= 2*a
                         endif
               else 
                         if( r.gt.-1) then 
                               psi= 2*am1*r
                               
                         else
                              psi= -2*am1
                         endif

               endif
c
c note weights supplied to cubic spline routine follow the convention that
c they are in terms of standard deviations. ( The more common form is
c   as reciprocal variances
c
        wt(k) = dsqrt( 2*r/psi)
  100   continue
        return
        end
  
  
  
         



