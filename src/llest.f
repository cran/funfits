c**** subroutine to loop over stretches of size  nprod
      subroutine
     *  llest(nlags, ldder, der, n,nprod,nll,llesv, lleqr,lle11,
     *   state)
      integer state
      double precision der(ldder,n),llesv(ldder,n),
     *     lleqr(ldder,n),lle11(n)
       nll= n-nprod+1

       do 10 k=1,nll
      
      call liapexp( nlags, ldder, der(1,k),
     *             nprod,llesv(1,k), lleqr(1,k), lle11(k), state)  


 10     continue
      return
      end







