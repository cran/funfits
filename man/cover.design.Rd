\name{cover.design}
\title{
Generate a coverage design.
}
\usage{
cover.design(R, nd, nruns = 1, nn = T, num.nn = 100,
		fixed = NULL, scale.type = "unscaled", 
		R.center, R.scale,
		P = -20, Q = 20, start = NULL, 
		DIST = NULL, return.grid = T)
}
\arguments{
\item{R}{
An Nxd matrix consisting of the N candidate points describing a
d-dimensional experimental region.  The region need not be regular.
}
\item{nd}{
An integer giving the number of points in the generated design. If
fixed points are used in the design, then nd is the number of new
points.
}
\item{nruns}{
An integer value giving the number of optimal designs to compute, each
starting with a different set of random points. Upon convergence of
the algorithm, the optimum design may be different for different
starting configurations of points. The design with the minimum
coverage criterion is retained. The default is nruns=1.
}
\item{nn}{
Logical value specifying whether or not to consider only nearest
neighbors in the swapping algorithm. When nn=F, then the swapping
algorithm will consider all points in the candidate space. When nn=T,
then the swapping algorithm will consider only the num.nn closest
points for possible swapping. The default is to use nearest neighbors
only (nn=T).
}
\item{num.nn}{
An integer the specifies the number of nearest neighbors to be
considered in the swapping algorithm. For large candidate sets, we
recommend num.nn >= 100. Values of num.nn >= 200 will probably give
results equivalent to full exploration of the space. num.nn=150 is
often a good compromise between speed of convergence and goodness of
the design. The default is 100.
}
\item{fixed}{
A matrix or vector specifying points to be forced into the
experimental design. If fixed is a matrix, it gives coordinates of the
fixed points in the design. In this case fixed must be a subset of
R. If fixed is a vector, then fixed gives the row numbers from the
candidate matrix R that identify the fixed points. The number of
points to be generated, nd, is in addition to the number of points
specified by fixed.
}
\item{scale.type}{
A character string that tells how to scale the candidate matrix, R,
before calculating distances. The default is "unscaled" in which case
no scaling is done.  This is appropriate when all of the columns of R
are commensurate; for example, when R gives x and y in spatial
coordinates. When the columns of R are not commensurate, then it is
generally thought that an appropriate choice of scaling will provide a
better design. This would be the case, for example, for a typical
process optimization. Other choices for scale.type are "range", which
scales all columns of R to the range (0,1), "unit.sd", which scales
all columns of R to have 0 mean and unit standard deviation, and
"user", which allows a user specified scaling. We typically use the
option "range" for physical experiments.
}
\item{R.center}{
A vector giving the centering values if scale.type='user'.
}
\item{R.scale}{
A vector giving the scale values if scale.type='user'.
}
\item{P}{
A scalar value specifiying a parameter of the criterion to be
optimized.  It affects how the distance from a point x to a set of
design points D is calculated (see discussion below).  P=1 gives
average distance. P=-1 gives harmonic mean distance.  P=-Inf would
give minimum distance (not available as a value). We typically use
P=-20 (the default) as an approximation to minimum distance.  Values
of P>=0 give trivial designs (e.g., all points equal to the
centroid). As P gets large (negatively) points will tend to be more
spread out.
}
\item{Q}{
A scalar value specifying a parameter of the criterion to be
optimized.  It affects how distances from all points not in the design
to points in the design are averaged.  When Q=1 (the default), simple
averaging of the distances is employed. As long as P is large
(negatively), Q is generally thought to have more affect on the ease
of optimization than on the characteristics of the design. Q=Inf (not
available as a value) in combination with P=-Inf would give a
classical minimax design. We often choose P=-20, Q=20 as an
approximation to the minimax design.  Values of Q<=0 give trivial
designs.
}
\item{start}{
A matrix or vector giving the initial design from which to start
optimization. If start is a matrix, it gives the coordinates of the
design points. In this case start must be a subset of R. If start is a
vector, then start gives the row numbers of the initial design. The
default is to use a random starting design.
}
\item{DIST}{
An S-PLUS function to be used for calculating the distances between
design points. The default is to use euclidian distance.
}
\item{return.grid}{
Logical value that tells whether or not to return the candidate matrix
as an attribute of the computed design. The default is return.grid=T.
The candidate matrix is used by plot.spatial.design if it is available.
}}
\value{
Returns a design object of class "spatial.design", which inherits from
"data.frame" and "matrix". The design object has the following attributes:

\item{best.id}{
Row numbers of the final design from the original candidate matrix, R.
}
\item{fixed}{
Row numbers of the fixed points from the original candidate matrix, R.
}
\item{opt.crit}{
Value of the optimality criterion for the final design.
}
\item{start.design}{
Row numbers of the starting design from the original candidate matrix, R.
}
\item{start.crit}{
Value of the optimality criterion for the starting design.
}
\item{history}{
The swapping history and corresponding values of the optimality 
criterion for the best design.
}
\item{other.designs}{
The designs other than the best design generated when nruns>1.
}
\item{other.crit}{
The optimality criteria for the other designs when nruns>1.
}
\item{DIST}{
The distance function used in calculating the design criterion.
}
\item{grid}{
The matrix R is returned if the argument return.grid=T.
}
\item{transform}{
The type of transformation used in scaling the data and the values
of the centering and scaling constants.
}
\item{call}{
The calling sequence.
}}
\description{
For nd design points di in the set D and nc candidate points ci in the set C, 
the coverage criteria is defined as:


  
   M(D,C) = \{SUM_\{ci in C\} 
            [sum_\{di in D\} (dist(di,ci)^P)]^(Q/P) \}^(1/Q)        [1]


Where P < 0, and Q > 0 are parameters.  The algorithm used in
"cover.design" to find the set of nd points in C that minimize this
criterion is an iterative swapping algorithm which will be described
briefly. The resulting design is referred to as a "coverage design"
from among the class of space-filling designs.


ALGORITHM:


An initial set of nd points is chosen randomly if no starting 
configuration is provided. The nc x nd distance matrix between the 
points in C and the points in D is computed, and raised to the power P.
The "row sums" of this matrix are computed.  Denote these as rs_\{i\} and
the vector of row sums as rs.  Using rs, M(D,C) is computed as:  


  (sum_\{i=1,nc\} (rs_\{i\})^(Q/P))^(1/Q)                               [2]


And note that if point d_\{i\} is "swapped" for point c_\{j\}, one must
only recompute 1 column of the original distance matrix, and 1 row.
The row elements not in the ith column will be the same for all j and
so only need computed when the first swapping occurs for each d_\{i\}.
Denote the sum of these off-i elements as "newrow(i)".  The index is
i here since this is the same for all rows j=1,\dots{}nc.


Thus, for each swap, the row sums vector is updated as 


  rs(new) = rs(old) - column(i,old) + column(i,new)


And the jth element of rs(new) is replaced by:


  rs(new)[j] = column(i,new)[k] + newrow(i)


Finally, M(D,C) is computed for this swap of the ith design point 
for the jth candidate point using [2].  The point in C that when 
swapped produces the minimum value of M(D,C) replaces d_\{i\}.  
This is done for all nd points in the design, and is iterated until 
M(D,C) does not change.


When the nearest neighbor option is selected, then the points
considered for swapping are limited to the num.nn nearest neighbors
of the current design point.


STABILITY 


The algorithm described above is guaranteed to converge.  However,
upon convergence, the solution is sensitive to the initial configuration
of points.  Thus, it is recommended that multiple optimizations be 
done (i.e. set nruns > 1 ).  Also, the quality of the solution depends
on the density of the points on the region.  At the same time, for
large regions (e.g. > 30 x 30 grids), optimization can be computationally
prohibitive using "cover.design" unless the nearest neighbor option
is employed.




}
\references{


Johnson, M.E., Moore, L.M., and Ylvisaker, D. (1990).  Minimax and
maximin distance designs. Journal of Statistical Planning and
Inference 26, 131-148.


SAS/QC Software. Volume 2: Usage and Reference. Version 6. First
Edition (1995). "Proc Optex". SAS Institute Inc. SAS Campus Drive,
Cary, NC 27513.


}
\seealso{
spread, expand.grid, plot.spatial.design
}
\examples{
# We wish to cover an 11 by 11 region with 9 points and take the
# best design out of 5 optimizations.
# This is a small candidate space, but we take num.nn=10 in order
# to speed computation for the example.
# There is no need to scale this example as x and y are commensurate.
grid11 <- expand.grid(list(x=1:11,y=1:11))
cover.des <- cover.design(grid11,9,nruns=5,num.nn=10)


## look at the design history
summary(cover.des)


# plot the region, and the optimum design:
plot(cover.des)


# display the change in the coverage criterion after each step of the
# optimization:
plot(x=attr(cover.des,'history')[,"step"],
	y=attr(cover.des,'history')[,"new.crit"])


}
\keyword{FUNFITS}
% Converted by Sd2Rd version 0.2-a3.


