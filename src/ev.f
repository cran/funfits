c
c    NB = nobs

c    NCTS= number of covariates plus dimension of null space
c    NTAB= related to number of points searched for GCV function
c
      parameter (NB=400,NVAR=4,NCTS=50)
      parameter (NW= NB*(NCTS+ 2+ NB)+NB,NWI=3*NB -NCTS)
      parameter (NB3=3*NB,NCOEF=NB+NCTS)
      parameter (NTAB=200)
      implicit double precision (a-g,o-z)
      integer nobs,info,ncov1,m,lddes,dim,
     * lwa,iwork(NWI),liwa,iout(4)
      double precision des(NB,NVAR),coef(NCOEF), work(NW)
      double precision des2(1,NVAR),f(1), ptab(100,NVAR)
c
      lddes=NB
      lddes2=1
      lwa=NW
      liwa=NWI
      ntbl = NTAB
      ldptab=100
c
c    get  previous estimation results from tps.ev
c
      call rest(1,des,lddes,nobs,dim,m,ncov1,iout,coef)

      nuobs=iout(4)
c
c
          if( nuobs.gt.lddes) then
           write(*,*) "Design matrix too large for arrays "
           write(*,*) " change NB in ev.f "
           stop
           endif
         ndes2= 1
         ib=0 

   5     read ( *,*, end=10) (des2(1,k),k=1,dim)
      call       ddevf(des2,lddes2,ndes2,dim,m,des,lddes,nuobs,
     *           ncov1,ncov2,coef,iout(2), f,work,lwa,iwork,info,ib)

c     evaluate spline not predicted values with covariates 
c
c             write(*,*) (des2(1,k),k=1,dim), f(1)
              write(*,100) f(1)
         goto 5
   10    continue

100   format(e20.14)

      stop
      end







