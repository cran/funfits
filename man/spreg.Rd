\name{spreg}
\title{
Smoothing spline regression
}
\usage{
spreg(x, y, lambda, xgrid, weight=rep(1, length(x)), derivative=0, 
Adiag=T, cost=1)
}
\arguments{
\item{x}{
Vector of x values 
}
\item{y}{
Vector of y values
}
\item{lambda}{
Smoothing parameter. If omitted this is estimated by GCV.
}
\item{xgrid}{
Vector of points to evaluate the estimated curve. Default is unique
sorted x's.
}
\item{weight}{
A vector that is proportional to the standard deviation of the errors.
}
\item{derivative}{
If equal to 1 or 2 returns the estimated first or second derivative of
the estimate
}
\item{Adiag}{
If true will compute leverage values for the estimate
}
\item{cost}{
Cost value to be used in the GCV criterion.
}}
\value{
A list of class spreg. The values of the GCV function and the
effective number of parameters are tabulated in the component gcv.grid.
The component predicted is a two column matrix that contains the
values from xgrid (or sorted unique x's) and the estimated curve at
these points. 
}
\description{
A smoothing spline is a locally weighted average of the data y's based
on the relative locations of the x values. Formally the estimate is
the curve that minimizes the criterion:
    (1/n) sum(k=1,n) ( Y_k - f( X_k))**2  + lambda* R(f)
where R(f) is the integral of the squared second derivative of f over
the range of the X values. The solution is a piecewise cubic
polynomial with the join points at the unique set of X values. The
polynomial segments are constructed so that the entire curve has
continuous first and second derivatives and the second and third
derivatives are zero at the boundaries.  The smoothing parameter has the range
[0,infinity]. Lambda equal to  zero gives a cubic spline interpolation
of  the data. As lambda diverges to infinity ( e.g lambda =1e20) the 
estimate will converge to the straight line estimated by least squares.


    The values of the estimated
function at the data points can be expressed in the matrix form:


 
    predicted.values= A(lambda)Y
where A is an nXn symmetric matrix that does NOT depend on Y.
The diagonal elements are the leverage values for the estimate and the
sum of these  (trace(A(lambda)) can be interpreted as the effective
number of parameters that are used to define the spline function. 


}
\references{
Additive Models by Hastie and Tibishirani
}
\seealso{
 predict.spreg, splint, tpsreg
}
\examples{
spreg( auto.paint$thick, auto.paint$DOI)-> out
plot(out)
lines(out$predicted)
}
\keyword{smoothing}
% Converted by Sd2Rd version 0.2-a3.
