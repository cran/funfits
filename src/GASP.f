
c**** subroutine to fill in the omega ( or K) matrix for 
c**** GASP type covariance
c**** K_ij= exp(-( sum_k( |x1_ik- x2_jk|**par(k)) )
C**** 
C*** in the code k is really ic
c
       subroutine gaspbs( nd,x1,n1, x2,n2, par, k)
       implicit double precision (a-h,o-z)
       integer nd,n1,n2,ic
       
       real*8 par(nd),x1(n1,nd), x2(n2,nd), k(n1,n2)
c **** loop through columns of output matrix K
c*** outer most loop over columns of x1 and x2 should reduce paging 

       do 5 ic= 1, nd
       do 10 j =1,n2
            xtemp= x2(j,ic)
            do  15 i= 1, n1
c
c** pth power differences

               k(i,j)=  (abs(x1(i,ic)- xtemp))**par(ic) + k(i,j)
 15             continue
 10    continue
 5      continue
c**** at this point k( i,j) is the 
c***** sum of differences each raised to the par(ic) power
c***
c*** now evaluate like radial basis functions
         nbig= n1*n2
c***** Now evalute the radial basis functions with the
c      distances. gaspf will just loop through the matrix 
c as stacked column vectors. 

         call gaspfn( nbig,k(1,1),par)

       return
       end
C** evaluates radial basis functions 
c**** K_ij= radfun( distance( x1_i, x2_j))
c
       subroutine multgb( nd,x1,n1, x2,n2, par, c,h,work)
       implicit double precision (a-h,o-z)
       integer nd,n1,n2,ic
       
       real*8 par(nd),x1(n1,nd), x2(n2,nd), c(n2), h(n1),sum
       real*8 work( n1), ddot

c****** work aray must be dimensioned to size n2
c **** loop through columns of output matrix K
c*** outer most loop over columns of x1 and x2 should reduce paging 

       do 5 ir= 1, n1

c
 
c evaluate all basis functions at  x1(j,.)       
       do 10 j =1,n2
c
c  zero out sum accumulator
c
         sum=0.0
      do 15  ic=1,nd
c
c** accumulate squared differences
c 

            sum= sum+ (dabs(x1(ir,ic)- x2(j,ic)))**par(ic)

 15             continue
        work(j)=sum
 10    continue

C**** evaluate squared distances  with basis functions. 

          call gaspfn( n2,work(1),par)
c
c***** now the dot product you have all been waiting for!
c
          h(ir)= ddot( n2, work(1), 1, c(1),1)
 5      continue

       return
       end

       subroutine gaspfn(n,d2, par)
       real*8 d2(n), par(1)
       integer n

         do 5 k =1,n
         d2(k)=  exp(-1*d2(k))
   5     continue

        return
        end
