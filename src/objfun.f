c -------------------------------------------------------------------------
c     Single-layer net, "quadratic" squasher
c         Inputs: X       (n-d) x d data matrix    (d lags, n data points)
c                 gama    K x d matrix             (K=# of hidden units)
c                 beta    K+1 x 1 column vector
c                 mu      K x 1 row vector
c         Output: gk=   beta * G(X*gama' .+ mu')
c -------------------------------------------------------------------------
        subroutine net(xout,xin,beta,amu,gama,nx,jd,k)
        implicit real*8(a-h,o-z)
        parameter(kmax=8,jdmax=16,nxmax=3000,
     *     hf=0.5,one=1.d0,two=2.d0)
        dimension xout(nxmax),xin(nxmax,jdmax),uin(nxmax,kmax)
        dimension gama(kmax,jdmax),beta(kmax+1),amu(kmax)
	do 10 it=1,nx
		xout(it)=beta(1)
10	continue
	do 50 jun=1,k
		u0=amu(jun)
		do 20 it=1,nx
			uin(it,jun)=u0
20 		continue

		do 40 lag=1,jd
			gjunlag=gama(jun,lag)
			do 30 it=1,nx
			   uin(it,jun)=uin(it,jun)+gjunlag*xin(it,lag)
30			continue
40		continue
50	continue
	do 100 jun=1,k
		do 90 it=1,nx
			uinij=uin(it,jun)
			auin=abs(uinij)
			xout(it)=xout(it)+beta(1+jun)*uinij
     *   		*(one+hf*auin)/(two+auin+hf*auin**2)
90		continue
100	continue
        return
        end
C=============================================================
        subroutine objfun(mode,npar,theta,f,g,nstate)
c
c           npar  = # of variables   (input integer)
c           theta = evaluation point (input vector of length nvar)
c           f     = function value   (output scalar)
c           g     = gradient value   (output vector of length nvar)
c           mode  ==  0 then compute f, leave g unchanged
c                 ==  1 then compute g, leave f unchanged
c                 ==  2 then compute both f and g
c          nstate = not used; included for compatibility with NPSOL

        implicit real*8(a-h,o-z)
        parameter(kmax=8,jdmax=16,nxmax=3000)
        parameter(npmax=1+kmax*(jdmax+2))
        parameter(zero=0.d0,one=1.d0,hf=0.5d0,two=2.d0)
        common /size/ k,jd,nx,m,jforce
        common /xdata/ xt,xlag
        dimension gama(kmax,jdmax)
        dimension xt(nxmax),xout(nxmax),xlag(nxmax,jdmax)
        dimension theta(npar),g(npar)

        m0=k+1
        j0=2*k+1
	jdf=jd+2*jforce

c ------  extract gamma matrix from theta vector
        do 30 iunit=1,k
            loc0=j0+(iunit-1)*jdf
            do 20 lag=1,jdf
                gama(iunit,lag)=theta(loc0+lag)
20          continue
30      continue

c ---------- compute predicted values (xout)
        call net(xout,xlag,theta(1),theta(m0+1),gama,nx,jdf,k)

c----------- compute objective function f (mean square prediction error)
        sse=zero
        if (mode.eq.0) then
            do 50 i=1,m
                sse=sse+(xout(i)-xt(i))**2
50          continue
            f=sqrt(sse/dfloat(m))
        else

c ---------- compute objective function & its gradient vector
	    do 55 j=1,npar
		g(j)=zero
55	    continue
            do 90 i=1,m
	      fi=xout(i)-xt(i)
              sse=sse+fi**2
	      g(1)=g(1)+fi
              do 80 j=1,k
                uin=theta(m0+j)
                do 60 lag=1,jdf
                    uin=uin+gama(j,lag)*xlag(i,lag)
60              continue
                auin=abs(uin)
                dij=(two+auin +hf*(auin**2))
	        g(j+1)=g(j+1)+fi*uin*(one + hf*auin)/dij
                gpij=two*(one+auin)/(dij**2)
		g(m0+j)=g(m0+j)+fi*theta(1+j)*gpij
                loc0=j0+(j-1)*jdf
                do 70 lag=1,jdf
                   g(loc0+lag)=g(loc0+lag)
     &                          +fi*theta(1+j)*xlag(i,lag)*gpij
70              continue
80            continue
90         continue
       
           fac=one/sqrt(sse*dfloat(m))
           do 120 j=1,npar
               g(j)=fac*g(j)
120        continue
           if (mode.eq.2) f=sqrt(sse/dfloat(m))
        endif
        return
        end
