\name{sreg}
\title{
Smoothing spline regression
}
\usage{
sreg(x, y, lam=NA, offset=0, wt=rep(1, length(x)), cost=1, nstep.cv=50, maxit.cv=10, xgrid=sort(unique(x)), deriv=0, find.trA=T)
}
\arguments{
\item{x}{
Vector of x values
}
\item{y}{
Vector of y values
}
\item{lam}{
Smoothing parameter. If omitted this is estimated by GCV.
}
\item{offset}{
GCV is RSS/((1-(tr(A)-offset)*cost + offset)/n)^2, so that the degrees of
freedom can be adjusted with the offset.
}
\item{wt}{
A vector that is proportional to the standard deviation of the errors.
}
\item{cost}{
Cost value to be used in the GCV criterion.
}
\item{nstep.cv}{
Number of grid points for minimum GCV search
}
\item{maxit.cv}{
Maximum number of iterations for Golden Section search of optimum
}
\item{xgrid}{
Vector of points to evaluate the estimated curve. Default is unique
sorted x's.
}
\item{deriv}{
If equal to 1 or 2 returns the estimated first or second derivative of
the estimate. 
}
\item{find.trA}{
Calculate the trace of A
}}
\value{
Returns a list of class sreg. This includes the predicted values and 
residuals. The lambda and the effective number of parameters for the fit are
also returned. The results of the grid search to minimize GCV are returned
in cv.grid.
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
derivatives are zero at the boundaries.  The smoothing has the range
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
splint, tps
}
\examples{
sreg(ozone$x[,1],ozone$y)-> fit # fit sreg to ozone at longitude values
plot(fit)                       # plot sreg fit
}
\keyword{smoothing}
% Converted by Sd2Rd version 0.2-a3.






