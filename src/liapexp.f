      subroutine liapexp(nlags,ldder,der,n,lexpsv,lexpqr,lexp11,state)

c Purpose: create the Jacobian and update the product
c       of Jacobians with the current vector of partial 
c       derivatives; using this product the
c       Lyapunov exponents are estimated
c
c On Entry:
c   nlags                number rows in jac
c   der(nc,n)            matrix of partial derivatives
c   ldder                leading dimension (#rows) of der in calling routine
c   
c On Exit:
c   lexpsv(nlags)          arrays of ordered estimated Lyapunov exponents, from
c   lexpqr(nlags)             SVD and QR decompositions, respectively
c   lexp11             magnitude of (1,1) entry
c Subprograms Called Directly:
c       Linpack - dsvdc,dqrdc
c Linked with Double-precision BLAS

      parameter(LMAX=20,LDW=6*LMAX)
      INTEGER nlags,info, state

      double precision jac(lmax,lmax),
     * s(lmax,lmax),t(lmax,lmax),
     * v(lmax,lmax),u(lmax,lmax),e(LMAX),epsln,
     * der(ldder,n),lexpsv(nlags),cnorm,w(LMAX),norm,work(LDW),
     * lexpqr(nlags),lexp11,qraux(LMAX),one,zero
       integer jpvt(LMAX)
 	parameter(one=1.d0, zero=0.d0,epsln=1.0d-15)

      if (LMAX.lt.ldder) then
	  write(3,*) ' Not enough storage in liapexp'
	  stop
      endif

c---------------- Initialize S=identity matrix
      call setidmat(s,lmax)

      CNORM=zero
c--------------------- Loop over kk=1,2,.. n to compute T_n
      do 300 kk=1,n

      if( state.eq.1) then
c--------------------- Put derivatives at time kk into Jac
      ind=1
      DO 130 I=1,nlags
	  DO 120 J=1,nlags                                          
	      jac(I,J)=der(ind,kk)                                      
	      ind=ind+1                                     
120           CONTINUE                                        
130       CONTINUE                                           
                                           
c--------------------- Update T(kk-1) by multiplying by Jac(kk)

      DO 13 I=1,NLAGS
	  DO 12 J=1,NLAGS                                     
	      sum=zero

	      DO 11 K=1,NLAGS                               
		  sum=jac(I,K)*S(K,J)+sum              
11            CONTINUE
              T(i,j)=sum
12        CONTINUE
13    CONTINUE
      else

c
c---- special multiplication for reconstructed state space form
c 
      DO 23 I=1,NLAGS
	      sum=zero

	      DO 24 K=1,NLAGS                               
		  sum=der(K,kk)*S(K,I)+sum              
24            CONTINUE
              T(1,I)=sum
23        CONTINUE
         if( nlags.gt.1) then
         do 25 I=2,nlags
           do 26 J=1,nlags
            T(I,J)= S(I-1,J)
 26        continue
 25      continue
         endif

      endif


c------------------ Normalize to Max(abs(T_i,j)) == 1
      NORM=zero
      DO 15 I=1,NLAGS                                             
	  DO 14 J=1,NLAGS                                        
	      NORM=MAX(NORM,ABS(T(I,J)))                        
14            CONTINUE                                         
15        CONTINUE                                            

      DO 17 I=1,NLAGS                                                 
	  DO 16 J=1,NLAGS                                            
	      T(I,J)=T(I,J)/NORM
c             copy T to S for next multiplication
	      S(I,J)=T(I,J)
16            CONTINUE                                             
17        CONTINUE                                                

c---------------- Update Cnorm=cumulative log of renormalizations
      CNORM=LOG(NORM)+CNORM                                     
    
300   continue
c------------- end loop over products





c     calculate (1,1) entry
      lexp11 = (LOG(ABS(S(1,1)))+CNORM)/dfloat(n)	

c------------- Estimate of LE's based on SV decomposition
      call dsvdc(t,LMAX,nlags,nlags,W,E,U,LMAX,V,LMAX,         
     *           work,00,info)

      if( info.gt.0) then
	   write(3,*) 'info from svddc', info
      endif

      do 210 k=1, nlags
	if( abs(w(k)).gt.epsln) then
	   lexpsv(k)= (LOG(ABS(w(k)))+CNORM)/dfloat(n)
	else
	   lexpsv(k)=(log(epsln) +cnorm)/dfloat(n)
       endif
210   continue

c---------- Estimate of LE's based on QR decomposition
      call dqrdc(S,LMAX,NLAGS,NLAGS,QRAUX,JPVT,WORK,0)

      do 220 k=1,nlags
c--- minus signs added so taht sort is reversed order
	if(abs(S(k,k)).gt.epsln) then
	   lexpqr(k)= -(LOG(ABS(S(k,k)))+CNORM)/dfloat(n)
	else
	   lexpqr(k)= -(log(epsln)+cnorm)/dfloat(n)
	endif
220   continue

      if(nlags.gt.1) then
	call hsort(nlags,lexpqr)
        endif
c
c-- fix up the signs
c
	do 230 k=1,nlags
	  lexpqr(k) = lexpqr(k)* -1

230     continue


      RETURN                                                          
      END

      subroutine hsort(n,k)
C  HEAPSORT ALGORITHM FOR SORTING ON VECTOR OF KEYS K OF LENGTH N    
C  J F MONAHAN        TRANSCRIBED FROM KNUTH, VOL 2, PP 146-7.       
      double precision K(1),KK                                                 
      INTEGER R                                                      
      IF(N.LE.1) RETURN                                              
      L=N/2+1                                                        
      R=N                                                            
  2   IF(L.GT.1) GO TO 1                                             
      KK=K(R)                                                        
      K(R)=K(1)                                                      
      R=R-1                                                          
      IF(R.EQ.1) GO TO 9                                             
      GO TO 3                                                        
  1   L=L-1                                                          
      KK=K(L)                                                        
  3   J=L                                                            
  4   I=J                                                            
      J=2*J                                                          
      IF(J-R) 5,6,8                                                  
  5   IF(K(J).LT.K(J+1)) J=J+1                                       
  6   IF(KK.GT.K(J)) GO TO 8                                         
  7   K(I)=K(J)                                                      
      GO TO 4                                                        
  8   K(I)=KK                                                        
      GO TO 2                                                        
  9   K(1)=KK                                                        
      RETURN                                                         
      END   





      subroutine setidmat( a,m)
      real*8 a(m,m)
      do 10 j=1, m
      do 20 k= 1,m
      a(j,k) =0.0
   20 continue
   10 continue
      do 30 j=1,m
          a(j,j)=1.0
   30 continue
      return
      end






