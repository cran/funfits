c------- computes an additve model using a cubic smoothing spline and 
c        the backfitting algorithm
c 
      subroutine addreg(x,ldx,nvar,nobs,y,wght,
     +  lam,trace,sxy,dsxy,sy,din,dout,job,ierr)
      parameter(  NMVAR=20,NMAX=10000,NGCV=100*3)
      implicit double precision (a-h,o-z)
      real*8 x(ldx,nvar), y(nobs), wght(nobs), 
     +               sxy(ldx,nvar),dsxy(ldx,nvar)

      real*8 sy(nobs), din(*), dout(*), lam(1)
      real*8 trace(1)
      real*8 ttrace(NMVAR),ltemp(NMVAR)
      real*8 syold(NMAX), stemp(NMAX),work(NMAX),xg(1),yg(1)
      real*8 diag(NMAX), thold, vlam,gcvout(NGCV)
      integer ierr, ngcv,job, jobspl(3)
c set up some constants
c     error return code
      ierr =0
      jobspl(1)=0
      jobspl(2)=0
      jobspl(3)=0
      ideriv=0
      ngrid=0
c      tolerance for convergence criterion
      tol=din(1)

c     maximum number of backfits
      maxit= int( din(2))
c parameters for gcv search
      nitcv=int( din(3))
      nstep= int( din(4))
c parameters for GCV function

      cost= din(5)
c set up start values for traces if included 
      do 15 k= 1, nvar
        ttrace(k)= trace(k)
 15   continue

c      offset= din(6)

      test=0.0
c check size of arrays

      if(  
     *   (nstep.gt.NGCV).or.(nvar.gt.NMVAR).or.(nobs.gt.ldx)
     *      .or.(nobs.gt.NMAX).or.(ldx.gt.NMAX)) then
      ierr=100
      return
      endif
c     write(*,*) tol, maxit
c zero out arrays, otherwise the passed values are 
c     used for starting     
c

         
      if( job.eq.0) then
          do 100 k=1,nobs
             sy(k)=0.0
        
  100     continue
          do 200 j=1,nvar
            
            do 250 k=1,nobs
               syold(k)=0.0
               sxy(k,j)=0.0         
  250     continue
  200     continue
      
      else
         do 260 k=1,nobs
            syold(k)= sy(k)
  260    continue
      endif
c set up work value and find stand. dev. of the y's
       ss=0.0
       sum=0.0
      do 300 k=1,nobs      
      ss= ss+ y(k)**2
      sum=sum+ y(k)
  300 continue
      sum= sum/nobs

      vary= ss/nobs - ((sum)**2)
c begin outer backfit loop
      do 1000 ib=1,maxit

c begin inner loop through variables
         do 1005 id=1,nvar
c          subtract off other  predicted components from y besides variable id
              do 1100 k=1,nobs
                  work(k)= y(k)- sy(k) + sxy(k,id)
c                  write(*,*)  id,x(k,id), work(k)
 1100        continue
        if( lam(id).gt.1e-20) then
 
        ltemp(id)=  log(lam(id))
        endif
c**** test for estimate of smoothting parameter
       if( lam(id).eq.-1) then
c**** calculate offset
        tsum=0.0
        do  1105 kk=1,nvar
         tsum= ttrace(kk) + tsum

 1105    continue
c*** if this is first pass then coerce search not to saturate
c** on one variable. This is done by adjusting the offset

        if (ib.eq.1) then

        offset= nobs- (nobs/cost)/(nvar)  +1 
        else
        
        offset= (tsum - 2*nvar -ttrace(id))*cost + 2*(nvar-1)
        endif
c       offset=0.0
c       write(*,*)ib, id, offset,cost
        call gcvcss( nobs,x(1,id),work,wght,cost,offset,
     -       nstep,nitcv,hopt,vopt,
     -     tropt,nstep,gcvout,
     -                 igcv)
 1500   format( 5e15.4)
       ltemp(id)= hopt
       ttrace(id)= tropt
c   some error codes from gcvcss are not fatal ( -1  through 2)
       if( igcv. ge.10) then
c but > 10 is fatal ! 
       ierr= igcv+1000*id
       return 
c      write(*,*) "error in call to gcvcss"
       endif
       endif
       call css(ltemp(id),nobs,x(1,id),work,wght,stemp,thold,diag,vlam,  
     +                  ngrid,xg,yg,jobspl,ideriv,ierr) 
              if( ierr.ne.0) then
              ierr= ierr + 100
c            write(*,*)  ierr, " error in call to css"
              return 
              endif
c         update predicted values and the id component
c            write(*,*) " new values for component:",id
c            write(*,*) ltemp(id),ttrace(id)
              do 1200 k=1,nobs
                 sy(k)=  sy(k) - sxy(k,id) + stemp(k)
                 sxy(k,id)= stemp(k)
c                write(*,*) ib,id, x(k,id), sxy(k,id)
1200         continue
c   break out of loop if there is only one variable to fit
        if( nvar.eq.1) goto 2000

 1005   continue
c       end inner loop over variables

c compute convergence criterion

         test=0.0
c     write(*,*) "old and updated fitted values"
      do 1300 k=1,nobs

         test= ( syold(k)-sy(k))**2 + test
c       write(*,*) ib, syold(k), sy(k),y(k)
         syold(k)= sy(k)
 1300 continue
      test= ( test/ vary)/nobs
c      write(*,*)   " end of backfit=",ib, " R**2 test value=", test

c break if test is less than tolerance or if only fitting one variable
      if( nvar.eq.1) goto 2000
      if( test.le.tol) goto 2000
 1000 continue
c     end of backfitting loop

 2000 continue
      dout(1)= ib
      dout(2)= test
c find traces 
      jobspl(1)=1
      jobspl(2)=3
      jobspl(3)=0
      ideriv=1
      ngrid=nobs
      do 3000 id= 1,nvar
              do 3100 k=1,nobs
                  work(k)= y(k)- sy(k) + sxy(k,id)
 3100        continue
       call css(ltemp(id),nobs,x(1,id),work,wght,stemp,thold,diag,vlam,  
     +                  ngrid,x(1,id),syold,jobspl,ideriv,ierr) 

        do 3200 k=1, nobs
          dsxy(k,id)= syold(k)
 3200   continue

        dout(2+id)=exp(ltemp(id))
        trace(id)= thold

              if( ierr.ne.0) then
              ierr= ierr + 100
c             write(*,*)  ierr, " error in call to css"
              return 
              endif

 3000 continue
      return
      end
