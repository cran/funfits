        subroutine bfgsfm(fcn,nvar,np,p,h0,ftol,fltol,maxiter,
     *                  maxstep,fnew,iter,iprint,inform,ibrent)

        implicit real*8(a-h,o-z)
        character*1 trans
        parameter(nvmax=100,eps=1.d-16,one=1.0,zero=0.0,trans='N')
        dimension p(np),pnew(nvmax),h0(np,np),g(nvmax),
     *            gnew(nvmax),h(nvmax,nvmax),xi(nvmax),work(nvmax),
     *            dx(nvmax),shy(nvmax),shys(nvmax,nvmax),
     *            pc(nvmax),gc(nvmax)
        external fcn
        sftol=sqrt(ftol)
        cftol=ftol**0.35
        ymin=-0.001
        inform=0
        jset=1

        call fcn(2,nvar,p,fp,g,0)
        do 10 i=1,nvar
            do 5 j=1,nvar
                h(i,j)=h0(i,j)
5           continue
10      continue

c---------------- set xi = -h*g
        call dgemv(trans,nvar,nvar,-1.d0,h,nvmax,g,1,zero,xi,1)

        nfail=0
        jhct=0
        do 99 iter=1,maxiter
12          call bhhhstp(fcn,nvar,p,fp,xi,g,maxstep,s,fnew,iters,info)

c---------------- check for failed step; reset hessian & search direction
            if (fnew.gt.fp) then
                if (jset.eq.1) then
                    inform=1
                    fnew=fp
                    goto 999
                endif
                nfail=nfail+1
                jset=1
                do 15 i=1,nvar
                    do 14 j=1,nvar
                        h(i,j)=h0(i,j)
14                  continue
15              continue

                call dgemv(trans,nvar,nvar,-1.d0,h,nvmax,g,1,zero,xi,1)
                goto 12
            else

c---------------- set pnew=p+s*xi and compute new gradient
            call dcopy(nvar,p,1,pnew,1)
            call daxpy(nvar,s,xi,1,pnew,1)
            call fcn(1,nvar,pnew,fjunk,gnew,0)

            if (ibrent.eq.1) then
c-------------- try to bracket minimum
            ax=0.
            bx=s
            fb=fnew
	    step=s
            do 16 iterb=1,5
                cx=bx+step
                call dcopy(nvar,p,1,pc,1)
                call daxpy(nvar,cx,xi,1,pc,1)
                call fcn(0,nvar,pc,fc,gc,0)
                if (fc.gt.fb) then
                    goto 17
                else
                    ax=bx
                    fb=fc
                    bx=cx
		    step=2.*step
                endif
16          continue
17          continue

c------------- If bracketing succeeded, use Brent's method to find minimum
            if (fc.gt.fb) then
               fnew=ffmin(fcn,nvar,p,xi,ax,bx,cx,fb,fltol,xmin)
               call dcopy(nvar,p,1,pnew,1)
               call daxpy(nvar,xmin,xi,1,pnew,1)
	       jbr=jbr+1

c--------------If bracketing failed, use lowest value found
            else
                s=cx
                fnew=fc
                call dcopy(nvar,pc,1,pnew,1)
           endif

          endif

c------------------ print out results every iprint iteration
            itermod=mod(iter,iprint)
            if (itermod.eq.0) then
                imax=idamax(nvar,gnew,1)
                gnorm=ddot(nvar,gnew,1,gnew,1)
                write(*,21)iter,fnew,(fp-fnew)/abs(fp),gnorm,
     *                      nfail,jhct,jbr
21              format(1x,i8,2x,f12.8,2x,e13.6,2x,f12.8,2x,3i4)
                nfail=0
                jhct=0
		jbr=0
            endif

c----------------- Check if convergence criteria satisfied
            imax=idamax(nvar,gnew,1)
            gmax=abs(gnew(imax))
            if(gmax.lt.eps) then
                goto 999
            else
                if ( (fp-fnew).lt.ftol*(1.+fnew)) then
                    if(gmax.lt.cftol*(1.+abs(fnew))) then
                        imax=idamax(nvar,pnew,1)
                        pmax=abs(pnew(imax))
                        call dcopy(nvar,pnew,1,work,1)
                        call daxpy(nvar,-1.d0,p,1,work,1)
                        imax=idamax(nvar,work,1)
                        dpmax=abs(work(imax))
                        if (dpmax.lt.sftol*(1.+pmax)) then
                            goto 999
                        endif
                    endif
                endif
            endif

c--------------------- If not, prepare for the next iteration
c        set dx=pnew-p
            call dcopy(nvar,pnew,1,dx,1)
            call daxpy(nvar,-1.d0,p,1,dx,1)

c        set shy=dx-h*gnew+h*g
            call dcopy(nvar,dx,1,shy,1)
            call dgemv(trans,nvar,nvar,-1.d0,h,nvmax,gnew,1,1.d0,shy,1)
            call dgemv(trans,nvar,nvar,1.d0,h,nvmax,g,1,1.d0,shy,1)

c        set shys=(shy*dx')+(shy*dx')'
            do 30 i=1,nvar
                do 25 j=1,nvar
                    shys(i,j)=shy(i)*dx(j) + shy(j)*dx(i)
25              continue
30          continue

c        set yshy=gnew'*shy-g'*shy
            y1=ddot(nvar,gnew,1,shy,1)
            y2=ddot(nvar,g,1,shy,1)
            yshy=y1-y2

c        set ys=gnew'*dx-g'*dx
            y1=ddot(nvar,gnew,1,dx,1)
            y2=ddot(nvar,g,1,dx,1)
            ys=y1-y2

c       update h if ys>0, otherwise reset h
            if (ys.gt.0.) then
                jset=0
                call dger(nvar,nvar,-yshy/ys,dx,1,dx,1,shys,nvmax)
                do 40 i=1,nvar
                    do 35 j=1,nvar
                        h(i,j)=h(i,j)+shys(i,j)/ys
35                  continue
40              continue
            else
                jset=1
                do 50 i=1,nvar
                    do 45 j=1,nvar
                        h(i,j)=h0(i,j)
45                  continue
50              continue
                jhct=jhct+1
            endif

c       update f,g,p, and search direction
            call dcopy(nvar,gnew,1,g,1)
            call dcopy(nvar,pnew,1,p,1)
            fp=fnew
            call dgemv(trans,nvar,nvar,-1.d0,h,nvmax,g,1,0.d0,xi,1)

c------------------- If next step "points up", reset Hessian
            y1=ddot(nvar,xi,1,g,1)
            y2=ddot(nvar,xi,1,xi,1)
            if (y1.ge.ymin*y2) then
                jset=1
                do 60 i=1,nvar
                    do 55 j=1,nvar
                        h(i,j)=h0(i,j)
55                  continue
60              continue
                call dgemv(trans,nvar,nvar,-1.d0,h,nvmax,g,1,0.d0,xi,1)
                jhct=jhct+1
            endif

        endif

99      continue

999     return
        end

        subroutine bhhhstp(fcn,nvar,x0,f0,d,g,itermax,s,fs,
     *                      iters,info)
        implicit real*8(a-h,o-z)
        parameter(delta=0.25,pwr=0.618,one=1.0,nvmax=100)
        dimension x0(nvar),d(nvar),g(nvar),work(nvmax)

c--------------------- Initializations
        info=1
        iter=0
        iup=0
        idown=0
        factor=2.
        dg=ddot(nvar,d,1,g,1)
        downfact=dg*delta
        upfact=dg*(1.-delta)
        alam=1.0
        if (itermax.eq.0) itermax=25

        do 10 iter=1,itermax

c----------------- Set g= alam*d + x0, vof=fcn(g)
            call dcopy(nvar,x0,1,work,1)
            call daxpy(nvar,alam,d,1,work,1)
            call fcn(0,nvar,work,vof,g,1)

c------------------ vof above cone, hence reduce step size
            if ((vof-f0).gt.(downfact*alam)) then
                idown=1
                if (iup.eq.1) then
                    factor=factor**pwr
                    iup=0
                endif
                alam=alam/factor

c------------------ vof below cone, hence increase step size
            else if ((vof-f0).lt.(upfact*alam)) then
                iup=1
                if (idown.eq.1) then
                    factor=factor**pwr
                    idown=0
                endif
                alam=alam*factor

c------------------ vof within cone, hence stop iterations
            else
                info=0
                goto 20
            endif
10      continue

20      continue

c------------------- compute returned values
        iters=iter
        s=alam
        fs=vof
        return
        end

C     function ffmin(ax,bx,cx,f,tol)
      function ffmin(fcn,nvar,x0,dx,AX,BX,CX,fb,TOL,XMIN)
      implicit real*8(a-h,o-z)
      parameter(nvmax=100,itmax=100)
      dimension x0(nvar),dx(nvar),work(nvmax),g(nvmax)
c
c      an approximation  x  to the point where  f  attains a minimum  on
c  the interval  (ax,bx)  is determined.
c
c
c  input..
c
c  ax    left endpoint of initial interval
c  bx    right endpoint of initial interval
c  f     function subprogram which evaluates  f(x)  for any  x
c        in the interval  (ax,bx)
c  tol   desired length of the interval of uncertainty of the final
c        result ( .ge. 0.0d0)
c
c
c  output..
c
c  fmin  abcissa approximating the point where  f  attains a minimum
c
c
c      the method used is a combination of  golden  section  search  and
c  successive parabolic interpolation.  convergence is never much slower
c  than  that  for  a  fibonacci search.  if  f  has a continuous second
c  derivative which is positive at the minimum (which is not  at  ax  or
c  bx),  then  convergence  is  superlinear, and usually of the order of
c  about  1.324....
c      the function  f  is never evaluated at two points closer together
c  than  eps*abs(fmin) + (tol/3), where eps is  approximately the square
c  root  of  the  relative  machine  precision.   if   f   is a unimodal
c  function and the computed values of   f   are  always  unimodal  when
c  separated by at least  eps*abs(x) + (tol/3), then  fmin  approximates
c  the abcissa of the global minimum of  f  on the interval  ax,bx  with
c  an error less than  3*eps*abs(fmin) + tol.  if   f   is not unimodal,
c  then fmin may approximate a local, but perhaps non-global, minimum to
c  the same accuracy.
c      this function subprogram is a slightly modified  version  of  the
c  algol  60 procedure  localmin  given in richard brent, algorithms for
c  minimization without derivatives, prentice - hall, inc. (1973).
c
C
C       Modified by S Ellner to use relative error tolerance, tol*(1.+dabs(x)),
C	instead of absolute error tolerance; and to assume a machine
C       epsilon of 10d-16. July 1992
C
C       This version is adapted for calling by bfgsi, and uses the
C       same conventions concerning the objective function.
c
c
c  c is the squared inverse of the golden ratio
c
      c = 0.5d0*(3. - dsqrt(5.0d0))
c
c  eps is approximately the square root of the relative machine
c  precision.
c
	eps=1.d-8
c
c  initialization
c
      a = min(ax,cx)
      b = max(ax,cx)
      v = bx
      w = v
      x = v
      e = 0.0d0
      fx = fb
      fv = fx
      fw = fx
c
c  main loop starts here
c
      do 86 iter=1,itmax
      xm = 0.5d0*(a + b)
      tol1 = eps + tol*(1.+dabs(x))
      tol2 = 2.0d0*tol1
c
c  check stopping criterion
c
      if (dabs(x - xm) .le. (tol2 - 0.5d0*(b - a))) go to 90
c
c is golden-section necessary
c
      if (dabs(e) .le. tol1) go to 40
c
c  fit parabola
c
      r = (x - w)*(fx - fv)
      q = (x - v)*(fx - fw)
      p = (x - v)*q - (x - w)*r
      q = 2.0d00*(q - r)
      if (q .gt. 0.0d0) p = -p
      q =  dabs(q)
      r = e
      e = d
c
c  is parabola acceptable?
c
   30 if (dabs(p) .ge. dabs(0.5d0*q*r)) go to 40
      if (p .le. q*(a - x)) go to 40
      if (p .ge. q*(b - x)) go to 40
c
c  a parabolic interpolation step
c
      d = p/q
      u = x + d
c
c  f must not be evaluated too close to ax or bx
c
      if ((u - a) .lt. tol2) d = dsign(tol1, xm - x)
      if ((b - u) .lt. tol2) d = dsign(tol1, xm - x)
      go to 50
c
c  a golden-section step
c
   40 if (x .ge. xm) then
          e = a - x
      else
          e = b - x
      endif
      d = c*e
c
c  f must not be evaluated too close to x
c
   50 if (dabs(d) .ge. tol1) then
          u = x + d
      else
          u = x + dsign(tol1, d)
      endif
      call dcopy(nvar,x0,1,work,1)
      call daxpy(nvar,u,dx,1,work,1)
      call fcn(0,nvar,work,fu,g,1)
C     fu = f(u)


c
c  update  a, b, v, w, and x
c
      if (fu .gt. fx) go to 60
      if (u .ge. x) then
          a = x
      else
          b=x
      endif
      v = w
      fv = fw
      w = x
      fw = fx
      x = u
      fx = fu
      go to 86
   60 if (u .lt. x) then
          a = u
      else
          b = u
      endif
      if (fu .le. fw) go to 70
      if (w .eq. x) go to 70
      if (fu .le. fv) go to 80
      if (v .eq. x) go to 80
      if (v .eq. w) go to 80
      go to 86
   70 v = w
      fv = fw
      w = u
      fw = fu
      go to 86
   80 v = u
      fv = fu
      go to 86
c
c  end of main loop
c
86    CONTINUE
90    xmin = x
      ffmin=fx
      return
      end
