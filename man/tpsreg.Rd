\name{tpsreg}
\title{
Thin plate spline regression
~function to do ???
}
\usage{
tpsreg(x, y, spar, m=2, clean=T)
}
\arguments{
\item{x}{
 A matrix of independent variables that are the arguments
}
\item{y}{
Vector of dependent variables
}
\item{spar}{
Value of the smoothing parameter. 
}
\item{m}{
Order of spline surface. 
}
\item{clean}{
Remove temporary files from the fitting process. 
}}
\value{
A list of class tpsreg. This includes the predicted surface in
fitted.value and the residuals in residual. The results of the grid
search to minimize the Generalized Cross Validation function is
returned in gcv.grid.  
}
\section{Side Effects}{
The computations are done by writing the data and job parameters to
temporary UNIX files and executing a stand alone FORTRAN program. This
strategy is preferred due to the memory requirements and the complexity
of the thin plate spline algorithms. The bulk of the computations are
done by a set of subroutines for thin plate spline problems: GCVPACK.
}
\description{
A thin plate spline is result of minimizing the residual sum of
squares subject to the constraint that the function have a certain
level of smoothness (or roughness penalty). Here smoothness is
quantified by the integral of squared m^th order derivatives. For one
dimension and m=2 the roughness penalty is the integrated square of
the second derivative of the function. For two dimensions the
roughness penalty is the integral of the square of
        Dxx(f) + 2Dxy(f) +Dyy(f)
 (where Duv denotes the second partial derivative with respect to u
and v.) Besides controlling the order of the derivatives, the value of
m also determines the base polynomial that will be fit to the data.
The degree of this polynomial will be (m-1).


The smoothing parameter controls the amount that the data is
smoothed. In the usual form this is denoted by lambda, the Lagrange
multiplier of the minimization problem. Although this is an awkward
scale, lambda =0 corresponds to no smoothness constraints and the data
is interpolated.  lambda=infinity corresponds to just fitting the
polynomial base model by ordinary least squares.  GCVPACK and this
function use a more convenient scale for the smoothing parameter which
in terms of lambda is log10(n*lambda). Thus interpolation and
smoothing correspond to the extremes -infinity and +infinity.  The
preferred is the effective number of parameters associated with the
fitted surface. This scale is a complicated but monotone transformation
of the smoothing parameter and these values are reported in gcv.grid
and the eff.df components.


 


 




}
\references{
See Additive Models by Hastie and Tibshriani.
}
\seealso{
predict.tpsreg, plot.tpsreg, summary.tpsreg


}
\examples{
Fitting a surface to ozone measurements. 
tpsreg(ozone$x, ozone$y) -> hold
plot(hold) # residual plots and a plot of the GCV function verses the 
           # effective number of parameters


 
}
\keyword{splines,}
% Converted by Sd2Rd version 0.2-a3.
