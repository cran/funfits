c----------------- HP version, August 1995 ---------
        program nnreg
        implicit double precision(a-h,o-z)
        parameter(kmax=8,jdmax=16,nxmax=3000)
        parameter(npmax=1+kmax*(jdmax+2))
        dimension ydat(nxmax),xmat(nxmax,jdmax),
     *   theta(npmax),ctheta(npmax),
     *   tsav(200,npmax),rsav(200),
     *  yhat(nxmax)
        dimension h0(npmax,npmax), xsd(jdmax), xmean(jdmax)
        character*35 fname,outname
        data glow, ghigh, scale / -1.26, 1.26, 0.5 /,
     *   cgcv2, maxstep1, maxstep2 / 2.0, 50, 250 /,
     *   ibrent, fltol / 0, 0.0 /   
        common /size/ k,jd,nx,m,jforce
        common /xdata/ ydat,xmat
        external objfun

c----------- Data statement initializes many parameters
c----------- GCV cost=2, ibrent and fltol=0 so never used 

c----------- Set initial Jacobian to identity matrix
	call setidmat(h0,npmax)
 
c----------- Read in parameters and minimization routine options
        open(11,file='nnreg.par',status='old')
        read(11,10) fname
        read(11,10) outname
10      format(a35)
        read(11,*) nx,nxc
        if( (nx.gt.nxmax).or.(nxc.gt.jdmax)) then
          write(*,*) "dataset too large!"
          stop
        endif
        read(11,*) ngrid, ntries, npol, iprout, igreed

c----------- Initialize random number generator
        read(11,*) iseed
	iseed=abs(mod(iseed,125000)) 
        rn=rng(iseed)

c----------- if ngrid is negative set switch to read in start values
        if( ngrid.lt.0) then
          istart=1
          npol=1
          ntries=1
          ngrid=1
        else
           istart=0
        endif
        if( (npol.gt.200).or.(npol.gt.ngrid)) then
          write(*,*)"too many values to polish"
          stop
               endif

         if( (igreed.gt.0) .and.(iprout.gt.0)) then
           write(*,*) "Greedy algorithm can only output the best fit"
           write(*,*) " Not all the polsihed results"
           stop
          endif

        read(11,*) ftol1, ftol2, itmax1, itmax2
        iprint=itmax1+itmax2

        read(11,*) k1,k2
        close(11)
        if( ( 1+ (nxc+2)*k2).gt.npmax) then 
          write(*,*) "too many parameters in the model for array size"
           stop
          endif

        if( ( k1.gt.kmax).or.( k2.gt.kmax)) then
         write(*,*) " too many hidden units"
        stop
        endif
          m= nx
          jd=nxc

c----------- read data series into xdat; set nx=total series length
       call getxy(fname,xmat,ydat,nxmax,nx,nxc)

c----------- Set up the output file
        open(12,file=outname,status='unknown',access='append')
          
        write(12,*) "# first row of data X, Y:"
        write(12,*) "#",(xmat(1,k), k=1,nxc), ydat(1)	
        write(12,*) "#"," last row of data X, Y:"
        write(12,*) "#",(xmat(nx,k), k=1,nxc), ydat( nx)	
        write(12,*) "# number of observations:", nx
        write(12,*) "# number of columns of X", nxc
        if( istart.eq.0) then
        write(12,*) "# number of points for grid search", ngrid
        write(12,*)"#",  glow, ghigh, ngrid 
        else
        write(12,*)"# starting values read in "
        write(12,*)"# a grid search has been omitted !"
        endif
        if( igreed.gt.0) then
        write(12,*) "# greedy algorithm used to fit net"
        write(12,*) "# in steps of ", igreed
        endif
        write(12,*) "# d  k   RMS        BIC      GCV    par    iter"
        close(12)

c------------ Standardize values to predict: SD=1, mean=0 
              do 500 ic=1,nxc
                 call sdev(xmat(1,ic),1,nx,tempsd, tempm)
                  xmean(ic)=tempm
                  xsd(ic)=tempsd
                 do 501 jj=1,nx
                        xmat(jj,ic)= (xmat(jj,ic)-tempm)/tempsd
  501            continue
  500         continue
                 
                 call sdev(ydat,1,nx,ysd, ymean)
                 do 502 jj=1,nx
                    ydat(jj) = (ydat(jj)-ymean)/ysd
  502            continue


c----------- beginning of output for S function

c----------- open file of start values
            if( istart.eq.1) then
               open(11,file='nnreg.start',status='old')
            endif

c----------- START LOOP on # units
        do 1004 ik=k1,k2
           k=ik
c----- reset number of hidden units to one for the greedy algorithm
           if( (igreed.gt.0).and.( ik.gt.k1)) then
                k=igreed
           endif
c----------- compute number of parameters

          npar=1+k*(jd+2+2*jforce)
          if (2*npar.gt.m) goto 1004
          if (npar.gt.npmax) goto 1004

c----------- if istart=1 then read in start values 

          if( istart.eq.1) then
           read(11,*) ( tsav(1,jpar), jpar=1,npar)
          endif
          do 15 i=1,200
              rsav(i)=10000.
15        continue

c----------- START LOOP on repeated fits

c----------- if start=0 then do random search for starting values

           if( istart.eq.0) then 

c----------- to match lenns glow= -1.26 ghigh =1.26 ngrid = 256
            step = (ghigh- glow)/(ngrid-1)
            do 99 irep=1,ngrid

c----------- initialize parameter choices
             eps=scale*10**( step*irep + glow )

             call tinit(theta,npar,eps,ntries)
                
c----------- do a high-tolerance fit
	     
              call bfgsfm(objfun,npar,npmax,theta,h0,ftol1,fltol,
     *              itmax1,maxstep1,fnew,iter,iprint,inform,ibrent)

             xm=float(m)
             bic=1.419+dlog(fnew)+0.5*float(npar)*dlog(xm)/xm
             imax=idamax(npol,rsav,1)
             rmax=rsav(imax)
             if (fnew.lt.rmax) then
                rsav(imax)=fnew
                do 85 jpar=1,npar
                   tsav(imax,jpar)=theta(jpar)
85              continue

9000            format( 5e16.8)
           
             endif
99        continue
        endif               
c----------- END LOOP on repeated fits
          do 999 i=1,npol
            do 997 jpar=1,npar
                    theta(jpar)=tsav(i,jpar)
997         continue

c----- entering subroutine theta are starting values from grinding

            call bfgsfm(objfun,npar,npmax,theta,h0,ftol2,fltol,
     *               itmax2,maxstep2,fnew,iter,iprint,inform,ibrent)
c----------- convert the paramteres to canonical form to avoid 
c----------- identification problems

            call canpar(theta,k, nxc, ctheta)

c----------- save new parameter estimates. 
            do 998 jpar=1,npar
                    tsav(i,jpar)= ctheta(jpar)
998         continue

                xm=float(m)
                bic=1.419+dlog(fnew)+0.5*float(npar)*dlog(xm)/xm
		gcv=(fnew/(1.-(cgcv*float(npar)/xm)))**2
               if( (fnew.lt.ftest).or.(i.eq.1)) then 
                      ibest=i
                      ftest=fnew
                endif

c----------- Write out results 
           if( iprout.gt.0) then
             write(*,8900)  nxc, k, fnew
             write(*,9000) (xmean(jj),jj=1,nxc)
             write(*,9000) (xsd(jj), jj=1,nxc)
             write(*,9000) ymean, ysd
             write(*,9000) (ctheta(jj), jj=1,npar) 

           endif 

                open(12,file=outname,status='old',access='append')
                write(12,98) jd,k,fnew,bic,gcv,npar,iter,
     *                       inform
		close(12)
97          format(1x,2i3,2x,f7.3,2(2x,f9.6),2x,i5,i3)
98          format(1x,2i3,2x,f9.6,2(2x,f7.4),i3,i7,i3,2(2x,f7.3))
         
999         continue
            
            do 1002 jpar=1,npar
                    ctheta(jpar)=tsav(ibest,jpar)
1002         continue

c----- write out the fit with smallest rms out of polished results

            if( iprout.eq.0) then
            write(*,8900)  nxc, k, ftest
 8900       format( i8, i8, e15.8) 

            write(*,9000) (xmean(jj),jj=1,nxc)
            write(*,9000) (xsd(jj), jj=1,nxc)
            write(*,9000) ymean, ysd
            write(*,9000) (ctheta(jj), jj=1,npar) 
            endif
            
c----- calculate residual vector and substitute for y if this
c----- is the greedy algorithm
            if( igreed.gt.0) then
            call netev( npar,ctheta,k,xmat, nxmax, nx,nxc,yhat)
            do 9100 jj=1,nx
c----- find residuals on raw scale 
              ydat(jj)= (ydat(jj)- yhat(jj))*ysd
 9100       continue
c----- reset the scale for next fit  
                 call sdev(ydat,1,nx,ysd, ymean)
                 do 9200 jj=1,nx
                    ydat(jj) = (ydat(jj)-ymean)/ysd
 9200            continue

           endif            

 1004    continue
c----------- END LOOP on #units
       if( istart.eq.1) close(11)

       open(12,file=outname,status='old',access='append')
       write(12,*) '@end  '
         
       close(12)
       stop
       end











