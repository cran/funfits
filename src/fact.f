      integer function fact(i)
      integer i
c
c Purpose: quick factorial function for the bspline routine
c	returns zero for negative i.
c
c On Entry:
c   i			a non-negative integer
c On Exit:
c   fact		i factorial
c
c $Header: /usr/local/cvsroot/funfits/src/fact.f,v 1.1.1.1 1998/05/24 21:50:08 agebhard Exp $
c
      integer j
      fact = 0
      if (i .ge. 0) fact = 1
      if (i .le. 1) return
      do 10 j = 2,i
	 fact = fact*j
   10	continue
      return
      end
