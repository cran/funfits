
C=============================================================
        subroutine  netev(npar,theta, k,xmat,nxmax,nx,nxc, yhat)
c
c           npar  = # of variables   (input integer)
c           theta = evaluation point (input vector of length npar)
c           f     = function value   (output scalar)
c           g     = gradient value   (output vector of length npar)
c           mode  ==  0 then compute f, leave g unchanged
c                 ==  1 then compute g, leave f unchanged
c                 ==  2 then compute both f and g
c           nstate is not used (optional flag for initialization)

        implicit double precision(a-h,o-z)
        parameter(kmax=8,jdmax=16)
        parameter(npmax=1+kmax*(jdmax+2))
        parameter(zero=0.d0,one=1.d0,hf=0.5d0,two=2.d0)
        dimension gama(kmax,jdmax)
        dimension  yhat(nxmax),xmat(nxmax,jdmax)
        dimension theta(npar)
        jd= nxc
        m0=k+1
        j0=2*k+1
        do 30 iunit=1,k
            loc0=j0+(iunit-1)*jd
            do 20 lag=1,jd
                gama(iunit,lag)=theta(loc0+lag)
20          continue
30      continue
        call net(yhat,xmat,theta(1),theta(m0+1),gama,nx,jd,k)

        return
        end 
