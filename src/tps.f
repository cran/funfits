
c  in free format and will estimate a smoothing spline
c for the additvie model
c    
c    y= f(t) + Xb  + e
c
c the input file 
c should have the parameters dim, m, ncov, job
c  on the first line
c
c dim   = dimension of t.
c m = order of roughness penalty.  ( m-1 degree polynomials in null
c                                                 space)
c ncov = dimension of b ( number of covariates)  
c job   = 1  for generalized cross validation
c         101 evaluate the smoothing spline at the value of the smoothing 
c          parameter given on the next line of input
c
c            
c the  name of the data file should follow
c
c format of data in the file  should be :
c t(1),t(2),...,t(dim), x(1),...x(ncov), y,   weight
c                       [ optional    ]     [optional]
c    NB = nobs
c    NCTS= number of covariates plus dimension of null space
c    NTAB= related to number of points searched for GCV function
c
      parameter (NDATA=400,NB=400,NVAR=2,MAXD=4)
      parameter(MMAX=4,ZN=35 + NVAR,NPAR= ZN+NB)
c   35 =  (MMAX + MAXD -1 choose  MAXD) this is the 
c maximum number of polynomial terms in the null space
c old      parameter (LWA1=ZN*(2*NB+1) +NPAR*(NB+NPAR))
c old      parameter (LWA2= (NPAR-ZN)*(NPAR- 2*ZN +2 +NB) +NPAR+NB)
c old      parameter (NW= LWA1 + LWA2)
      parameter( NW= NB*(2+ NPAR) +NDATA )
      parameter (NWI= NB +2*NDATA-ZN, NB3=3*NB,NCOEF=NPAR)
      parameter (NTAB= 100)
      integer nobs,info,ntbl,ncov,lds,job,m,ld,dim,
     * lwa,iwork(NWI),liwa,ldtbl,iout(4)
      double precision s(NDATA,NVAR),lamlim(2),des(NDATA,MAXD),
     * y(NDATA),
     * adiag(NDATA),tbl(NTAB,3),des2(NDATA,MAXD),s2(NDATA,NVAR),
     * coef(NCOEF),auxtbl(3,3),svals(NDATA),dout(5),
     * work(NW),trA,diag,wt(NDATA)
      double precision savy(NDATA)
      double precision sighat
      character*40 filen
c
      ld=NDATA
      ldtbl=NTAB
      lds=NDATA
      lwa=NW
      liwa=NWI
c     write(*,*) NW,lwa,liwa,NB,NDATA
      ntbl = 80 
c**** decompositions will be saved in case this program is altered
c**** and dntpss is called more than once with same design points
      jobrd=1

      call gpar(dim,m, ncov,job,lamlim, filen)
c
c  set up pointers to returned coefficients
c
      if(dim.gt.MAXD) then
                write(*,*) 'number of variables (d) is too large'
                stop
            endif
      if( 2*m-dim .lt.1) then
          write(*,*) " 2m -d must be greater than 0"
          stop
       endif
      
      
      nd= 0
      nc= nd +dim
      ny=nc+ncov
      nn= ny +1
c 

      call gdata(dim,ncov,nobs,y,wt,ld,des,lds,s,job,filen)

c save data 
        call dcopy( nobs,y,1,savy,1)

      do  10 j=1,dim
          call dcopy(nobs,des(1,j),1,des2(1,j),1)
   10 continue

      do  20 j=1,ncov
          call dcopy(nobs,s(1,j),1,s2(1,j),1)
   20 continue
c
c
      if( mod(job,10).eq.2) then
        write(*,*) 'Warning: Replicated observtions will not be '
        write(*,*) 'weighted correctly'
      endif
      write(*,*) 'call to dntpss'
      call dntpss(des,ld,nobs,dim,m,s,lds,ncov,y,wt,ntbl,adiag,lamlim,
     * dout,iout,coef,svals,tbl,ldtbl,auxtbl,work,lwa,iwork,liwa,
     * job,jobrd,info)
      nunq= iout(4)

      if (info .ne. 0) write(*,*) 'dtpss info',info
      write(*,*)  'nobs= ',nobs,'  lamhat = ',dout(1)
      write(*,*) 'number of unique observations = ',nunq
      
      write(*,*) 'log10(nobs*lambda)= ',auxtbl(1,1)
      write(*,*)  'penlty = ',dout(2)

       if ( dout(4).gt. 1.0d-8) sighat= dsqrt( dout(3)/dout(4))
      write(*,*) 'sigma hat = ',sighat

      write(*,*) 'effective degrees of freedom for residuals = ',dout(4)
      write(*,*) 'total number of conventional parameters = ',iout(3)

      trA=dble(nobs)-dout(4)
      if (abs(trA-diag) .gt. 1.0d-8) write(*,*)
     *  'effective number of parameters including spline: ',trA
c
c
c***** pointer to coefficients of covariates   
      ncv= iout(3)-ncov

      if(ncov.ne.0) then 
            write(*,*) 'coefficients for covariates',
     *                       (k,coef(ncv+k),k=1,ncov)
      end if
      open(2,file='tps.spl')
      rewind 2
 
       if( ncov.eq.0) then 
       write(*,*) ' FILE: tps.spline contains  ',
     *            ' spline, residual '
               else
       write(*,*) ' FILE: tps.spline contains ',
     *                  ' spline, predicted, residual '
                endif

      do  400 l=1,nobs
           pf=y(l)
           if( ncov.ne.0)  then
               do 450 k=1,ncov
                    pf= pf - s2(l,k)*coef(ncv+k)
  450          continue
           endif
           res= savy(l)-y(l)
c
           if( ncov.eq.0) then 
c                write(2,500) (des2(l,k),k=1,dim),y(l),res
                write(2,500) y(l),res
           else
c                write(2,500) (des2(l,k),k=1,dim),pf,y(l),res
                write(2,500) pf,y(l),res
           end if

  400 continue
      close(2) 

  500 format(e20.10)
      open(2, file='tps.par')
      rewind 2
      write(*,*) ' FILE: tps.gcvf contains'
      write(*,*) ' nobs,dim,m,ncov1,(iout(k),k=1,4),log10(nobs*lambda)' 
      write(2,*)  nobs,dim,m,ncov1,(iout(k),k=1,4)
      write(2,500)  auxtbl(1,1)
      write(2,500) trA
      write(2,500) (dout(jj), jj=1,4)
      close(2)

      open(2,file='tps.gcv')
      rewind 2

      write(*,*) ' FILE: tps.gcvf contains',
     *             ' log10(n*lambda), gcvf(lambda) , tr(A(lambda)'
      it= (job - 1000*( job/1000))/100
      if( it.eq.1) then 
        nrows=1
      else
         nrows=ntbl
      endif
      do 600 l=1,nrows
         write(2,500) tbl(l,1),tbl(l,2),tbl(l,3)
  600 continue

      close(2)

c
c  save estimate
c
      call rest(0,des,ld,nobs,dim,m,ncov,iout,coef)
      stop
      end


      subroutine gpar(dim,m,ncov,job,lamlim,filen)
      integer dim,ncov,job
      doubleprecision lamlim(2)
      character*40 filen
c			  
      read(*,*) dim,m, ncov,job
      write(*,*) '  dimension= ',dim,  ' m=' , m, '# of covariates= ',ncov
      write(*,*) ' job= ',job

      it= (job - 1000*( job/1000))/100
      if( it.eq.1) then
         read(*,*) temp
          lamlim(1)=temp
          lamlim(2)=temp
          write(*,*) 'ln(nobs*lambda)  used for evaluation : ' ,temp
      else 
          write(*,*) ' lambda found by GCV '
      endif
      read(*,*) filen
      write(*,*) 'data file= ', filen
       return
      end
c
c
      subroutine gdata( dim,ncov,nobs,y,wt,ld,des,lds,s,job,filen)

      double precision s(lds,1),des(ld,1)
      double precision y(1),wt(1),buff(20)
      character*40 filen
      integer dim

      i=0

       iwght= job- 10*(job/10) 
      nd= 0
      nc= nd +dim
      ny=nc+ncov
      nn= ny +1
      if(iwght.eq.2) nn=nn+1


      open(2, file=filen)
      rewind(2)

   10 continue


           read(2,*,end= 45) (buff(k),k=1,nn)
          i= i+1
                if (i.eq.1) then
                write(*,*) ' first record: ', (buff(k),k=1,nn)
                endif
c          write(*,*) ' i=' , i, buff(1)
                 
          do 20 k=1,dim
                des(i,k)= buff(k+nd)
   20     continue
          if(ncov.gt.0) then
               do 25 k=1,ncov
                    s(i,k)= buff(k+nc)
   25          continue
          endif
          y(i)= buff(ny+1)
          if(iwght.eq.2) wt(i)=buff(ny+2)
      goto 10
c
   45 continue
      nobs= i
c

                write(*,*) ' last record: ', (buff(k),k=1,nn)

      write(*,*)   nobs, ' observations read from', filen
      close(2)
      return
      end
c
c



