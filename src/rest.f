       subroutine rest( mode, des,lddes,nobs,dim ,m,ncov1, iout,coef)
       integer mode,dim,lddes,nobs,m,iout(4),ncov1
       real*8 des(lddes,1),coef(1)
c
       open( 2,file='tps.ev')
       rewind 2


c
       if( mode.eq.1) then 
       read( 2,*) nobs,dim,m,ncov1,(iout(k),k=1,4)
       if( iout(4).gt.lddes) then
        return
       endif
       read(2,*) ( (des(i,k),k=1,dim),i=1,iout(4))
       read(2,*)  ( coef(i),i=1,iout(2) )
       else
       write(*,*) 'FILE: tps.ev intermediate file for ev '
       write( 2,*) nobs,dim,m,ncov1,(iout(k),k=1,4)
       n= iout(4)
        write(2,100) ( (des(i,k),k=1,dim),i=1,iout(4))
       write(2,100) ( coef(i),i=1,iout(2) )
       endif
c
       close(2)
 100   format(e20.14)
c
       return
       end
