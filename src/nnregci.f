c----------------- HP version, November 1995 -----------
        program nnregCI
        implicit double precision(a-h,o-z)
        parameter(kmax=8,jdmax=16,nxmax=3000)
        parameter(npmax=1+kmax*(jdmax+2),zero=0.d0)
	parameter(btol=1.d-4)
        dimension ydat(nxmax),xmat(nxmax,jdmax),
     *   theta(npmax),ctheta(npmax),g(npmax),dtheta(npmax),
     *   tsav(500,npmax),rsav(500),tgsav(500,npmax),grsav(500)
        dimension h0(npmax,npmax),xsd(jdmax),xmean(jdmax)
        character*35 fname,outname
        data glow, ghigh, scale / -1.26, 1.26, 0.5 /,
     *   cgcv2, maxstep1, maxstep2 / 2.0, 50, 250 /,
     *   ibrent, fltol / 0, 0.0 / 
        common /size/ k,jd,nx,m,jforce
        common /xdata/ ydat,xmat
        common /backup/ theta,dtheta,fcut,npar
        external objfun,func

c----------- Data statement initializes many parameters
c----------- GCV cost=2, ibrent and fltol=0 so never used 
	
c----------- Set initial Jacobian to identity matrix
	call setidmat(h0,npmax)

c----------- Read in parameters and minimization routine options
        open(11,file='nnci.par',status='old')
        read(11,10) fname
        read(11,10) outname
10      format(a35)
        read(11,*) nx,nxc
        if( (nx.gt.nxmax).or.(nxc.gt.jdmax)) then
          write(*,*) "dataset too large!"
          stop
        endif
        read(11,*) ngrid,ntries,npol,rms,cut1,cut2,nfits

c----------- Initialize random number generator
        read(11,*) iseed
	iseed=abs(mod(iseed,125000)) 
        rn=rng(iseed)

        read(11,*) ftol1, ftol2, itmax1, itmax2
        iprint=itmax1+itmax2

        read(11,*) k1
        close(11)
        if( ( 1+ (nxc+2)*k1).gt.npmax) then
          write(*,*) "too many parameters in the model for array size"
           stop
        endif

        if( k1.gt.kmax) then
          write(*,*) " too many hidden units"
          stop
        endif
        m= nx
        jd=nxc

	open(13,file='nnregCI.cut',status='unknown',access='append')
	rewind(13)
	write(13,*) '## icount1   icount2'
	write(13,808) itmax1,itmax2,scale,ftol2,maxstep1,maxstep2

	write(13,809) ngrid,npol,ntries
	close(13)
808	format(1x,2i10,2(2x,e12.6),2i10)
809	format(1x,3i10)

c------------ read data series into xdat; set nx=total series length
       call getxy(fname,xmat,ydat,nxmax,nx,nxc)

c---------- Set up the output file
        open(12,file=outname,status='unknown',access='append')
	rewind(12)
        write(12,*) "# first row of data X, Y:"
        write(12,*) "#",(xmat(1,k), k=1,nxc), ydat(1)
        write(12,*) "#"," last row of data X, Y:"
        write(12,*) "#",(xmat(nx,k), k=1,nxc), ydat( nx)
        write(12,*) "# number of observations:", nx
        write(12,*) "# number of columns of X", nxc
        write(12,*) "# number of points for grid search", ngrid
        write(12,*)"#",  glow, ghigh, ngrid
        write(12,*) "#", nx, nxc,ngrid,ntries
        write(12,*)"#", npol,rms,cut1,cut2,nfits
        write(12,*)"#",iseed,rn
        write(12,*)"#",ftol1,ftol2,itmax1,itmax2
        write(12,*)"#", k1, iprint
        write(12,*) "# d  k   RMS        BIC      GCV    par    iter"
        close(12)

c------------ Standardize values to predict: SD=1,mean=0
        do ic=1,nxc
             call sdev(xmat(1,ic),1,nx,tempsd,tempm)
             xmean(ic)=tempm
             xsd(ic)=tempsd
             do jj=1,nx
                xmat(jj,ic)= (xmat(jj,ic)-tempm)/tempsd
             end do
	end do

        call sdev(ydat,1,nx,ysd, ymean)
        do jj=1,nx
             ydat(jj) = (ydat(jj)-ymean)/ysd
	end do

c------------------ START LOOP on # units
        do 1004 k=k1,k1

c----------compute number of parameters
          npar=1+k*(jd+2+2*jforce)
          if (2*npar.gt.m) goto 1004
          if (npar.gt.npmax) goto 1004

          do i=1,500
              rsav(i)=10000.
          end do

c-------set counters to zero
         icount1=0
         icount2=0
	 inpol=0
         icountz=0

c------------------ START LOOP on repeated fits
c to match lenns glow= -1.26 ghigh =1.26 ngrid = 251
            step = (ghigh-glow)/(ngrid-1)

995         continue
	    do j=1,500
	    	grsav(j)=10000.
	    end do

            do 99 irep=1,ngrid
c----------------- initialize parameter choices
              eps=scale*10**(step*dfloat(irep) + glow)
              call tinit(theta,npar,eps,ntries)

c------------------ do a high-tolerance fit
	      call bfgsfm(objfun,npar,npmax,theta,h0,ftol1,fltol,
     *              itmax1,maxstep1,fnew,iter,iprint,inform,ibrent)

C------------ Save for polishing
              imax=idamax(npol,grsav,1)
              rmax=grsav(imax)
              if (fnew.lt.rmax) then
                grsav(imax)=fnew
                do jpar=1,npar
                   tgsav(imax,jpar)=theta(jpar)
                end do
	      endif

99          continue

C------------ Polish
           do 199 i=1,npol
             inpol = inpol +1
             do jpar=1,npar
                 theta(jpar)=tgsav(i,jpar)
	     end do
             call bfgsfm(objfun,npar,npmax,theta,h0,ftol2,fltol,
     *               itmax2,maxstep2,fnew,iter,iprint,inform,0)
	     call objfun(0,npar,theta,fnew,g,1)

C----------- Generate distribution for cut; cut1 with probability .2
C----------- and uniform between cut1 and rms with probability .8
	     if (rng(1) .le. .2) then
                fcut = cut1
	     else
                fcut = rms + rng(1)*(cut1-rms)
             endif
             
             if (fnew .gt. fcut) then
                do j=1,npol
                  do jpar=1,npar
                      theta(jpar)=tgsav(j,jpar)
	          end do
                  call bfgsfm(objfun,npar,npmax,theta,h0,ftol2,fltol,
     *               itmax2,maxstep2,fnew,iter,iprint,inform,0)
	          call objfun(0,npar,theta,fnew,g,1)
                  if (fnew .le. fcut) goto 875
                end do

8001	     format(1x,'failed to get rms below cut ',5i5,3(2x,f8.5))
	     open(13,file='nnregCI.cut',status='old',access='append')
	     write(13,8001) i,iter,inpol,icount1,icount2,
     *        fnew,cut1,cut2
	     close(13)
             goto 199
              endif

875           continue

C----------- Generate theta value with large enough rms to use root
C----------- finder, zeroin
                sum=0.d0
                tnorm=0.d0
                do j=1,npar
                     sum = sum + theta(j)*theta(j)
                end do
                tnorm = sqrt(sum)
		do j=1,npar
		   dtheta(j)=rng(1)*tnorm
		end do
             alpha =0.005
                do j=1,15
                   alpha=2.d0*alpha
                   y=func(alpha)
                   if (y.gt.zero) goto 876
                end do
876             continue

		f1=func(zero)
		f2=func(alpha)
		if ( (f1.lt.zero).and.(f2.gt.zero) ) then
                  icountz =icountz +1
                  zalpha=zeroin(zero,alpha,func,btol)
		  do j=1,npar
		    theta(j)=theta(j)+zalpha*dtheta(j)
		  end do

c----------- Check accuracy of root finder (not necessary)
	          call objfun(0,npar,theta,fnew,g,1)
                  fdiff=fcut-fnew
                  fnew=fcut
                  rsav(icountz)=fnew
                  do jpar=1,npar
                    tsav(icountz,jpar)=theta(jpar)
                  end do
		
C----------- Counts in bin1 and bin2
                  if ((fnew.gt.cut2).and.(fnew.le.cut1)) then
                     icount1=icount1+1
                  endif

        	  if (fnew.le.cut2) then
	             icount2=icount2+1
                  endif

             endif

9001	     format(1x,'polished ',5i5,4(2x,f8.5))
	     open(13,file='nnregCI.cut',status='old',access='append')
     	     write(13,9001) i,iter,inpol,icount1,icount2,
     *        fnew,fdiff,cut1,cut2
	     close(13)

c----------- convert the parameters to canonical form
            do jpar=1,npar
               theta(jpar)=tsav(i,jpar)
	    end do
            call canpar(theta,k,nxc,ctheta)
            do jpar=1,npar
               tsav(i,jpar)= ctheta(jpar)
	    end do

            fnew=rsav(i)
            xm=float(m)
            bic=1.419+dlog(fnew)+0.5*float(npar)*dlog(xm)/xm
	    gcv=(fnew/(1.-(cgcv*float(npar)/xm)))**2

C----------- Output for S function
            write(*,8900)  nxc, k, fnew
            write(*,9000) (xmean(jj),jj=1,nxc)
            write(*,9000) (xsd(jj), jj=1,nxc)
            write(*,9000) ymean, ysd
            write(*,9000) (ctheta(jj), jj=1,npar)
8900        format(1x, i8, i8, e15.8)
9000         format(5e16.8)

c------------------ Write out results
            open(12,file=outname,status='old',access='append')
            write(12,98) jd,k,fnew,bic,gcv,npar,iter,inform
	    close(12)
97          format(1x,2i3,2x,f7.3,2(2x,f9.6),2x,i5,i3)
98          format(1x,2i3,2x,f9.6,2(2x,f7.4),i3,i7,i3,2(2x,f7.3))

             if ((icount1+icount2).eq.nfits) goto 996
199        continue
             if ((icount1+icount2).lt.nfits) goto 995

996     continue
c---------------------- END LOOP on repeated fits

1004    continue
c---------------------- END LOOP on #units

       open(12,file=outname,status='old',access='append')
       write(12,*)"#count is ",icount1,icount2

       close(12)
       stop
       end


       function func(alpha)
       implicit double precision(a-h,o-z)
       parameter(kmax=8,jdmax=16,nxmax=3000)
       parameter(npmax=1+kmax*(jdmax+2))
       dimension theta(npmax),dtheta(npmax),g(npmax),
     *	    atheta(npmax)
       common /backup/ theta,dtheta,fcut,npar
       do j=1,npar
          atheta(j)=theta(j)+alpha*dtheta(j)
       end do
       call objfun(0,npar,atheta,fnew,g,1)
       func=(fnew-fcut)/fcut
       return
       end

      function zeroin(ax,bx,f,tol)
      double precision zeroin,ax,bx,f,tol
      external func
c
c      a zero of the function  f(x)  is computed in the interval ax,bx .
c
c  input..
c
c  ax     left endpoint of initial interval
c  bx     right endpoint of initial interval
c  f      function subprogram which evaluates f(x) for any x in
c         the interval  ax,bx
c  tol    desired length of the interval of uncertainty of the
c         final result ( .ge. 0.0d0)
c
c
c  output..
c
c  zeroin abcissa approximating a zero of  f  in the interval ax,bx
c
c
c      it is assumed  that   f(ax)   and   f(bx)   have  opposite  signs
c  without  a  check.  zeroin  returns a zero  x  in the given interval
c  ax,bx  to within a tolerance  4*macheps*abs(x) + tol, where macheps
c  is the relative machine precision.
c      this function subprogram is a slightly  modified  translation  of
c  the algol 60 procedure  zero  given in  richard brent, algorithms for
c  minimization without derivatives, prentice - hall, inc. (1973).
c
c
      double precision  a,b,c,d,e,eps,fa,fb,fc,tol1,xm,p,q,r,s
      double precision  dabs,dsign
c
c  compute eps, the relative machine precision
c
      eps = 1.0d0
   10 eps = eps/2.0d0
      tol1 = 1.0d0 + eps
      if (tol1 .gt. 1.0d0) go to 10
c
c initialization
c
      a = ax
      b = bx
      fa = f(a)
      fb = f(b)
c
c begin step
c
   20 c = a
      fc = fa
      d = b - a
      e = d
   30 if (dabs(fc) .ge. dabs(fb)) go to 40
      a = b
      b = c
      c = a
      fa = fb
      fb = fc
      fc = fa
c
c convergence test
c
   40 tol1 = 2.0d0*eps*dabs(b) + 0.5d0*tol
      xm = .5*(c - b)
      if (dabs(xm) .le. tol1) go to 90
      if (fb .eq. 0.0d0) go to 90
c
c is bisection necessary
c
      if (dabs(e) .lt. tol1) go to 70
      if (dabs(fa) .le. dabs(fb)) go to 70
c
c is quadratic interpolation possible
c
      if (a .ne. c) go to 50
c
c linear interpolation
c
      s = fb/fa
      p = 2.0d0*xm*s
      q = 1.0d0 - s
      go to 60
c
c inverse quadratic interpolation
c
   50 q = fa/fc
      r = fb/fc
      s = fb/fa
      p = s*(2.0d0*xm*q*(q - r) - (b - a)*(r - 1.0d0))
      q = (q - 1.0d0)*(r - 1.0d0)*(s - 1.0d0)
c
c adjust signs
c
   60 if (p .gt. 0.0d0) q = -q
      p = dabs(p)
c
c is interpolation acceptable
c
      if ((2.0d0*p) .ge. (3.0d0*xm*q - dabs(tol1*q))) go to 70
      if (p .ge. dabs(0.5d0*e*q)) go to 70
      e = d
      d = p/q
      go to 80
c
c bisection
c
   70 d = xm
      e = d
c
c complete step
c
   80 a = b
      fa = fb
      if (dabs(d) .gt. tol1) b = b + d
      if (dabs(d) .le. tol1) b = b + dsign(tol1, xm)
      fb = f(b)
      if ((fb*(fc/dabs(fc))) .gt. 0.0d0) go to 20
      go to 30
c
c done
c
   90 zeroin = b
      return
      end

C From FMM (Forsythe, George E., Malcolm, Michael A., and 
C Moler, Cleve B) on netlib 
