\name{qspreg}
\title{
Quantile spline regression 
}
\usage{
qspreg(x, y, lam=NA, maxit=50, maxit.cv=10, tol=0.0001, cost=1,
 offset=0, sc=sqrt(var(y)) * 1e-07, alpha=0.5, wt=rep(1, length(x)),
 nstep.cv=50, xgrid=sort(unique(x)), deriv=0, hmin=-35, hmax=-15)
}
\arguments{
\item{x}{
Vector of independent variables
}
\item{y}{
Vector dependent variables
}
\item{lam}{
Values of the smoothing parameter. If omitted is found by GCV based on the 
the quantile criterion
}
\item{maxit}{
Maximum number of iterations used to estimate each quantile spline.
}
\item{maxit.cv}{
Maximum number of iterations to find GCV minimum.
}
\item{tol}{
Tolerance for convergence when computing quantile spline.
}
\item{cost}{
Cost value used in the GCV criterion. Cost=1 is the usual GCV 
denominator.
}
\item{offset}{
Constant added to the effective degrees of freedom in the GCV function. 
}
\item{sc}{
Scale factor for quantile function. Default is a scale on the order of 
machine precision. Scales on the order of the residuals will result is a 
robust regression fit using the Huber weight function. 
}
\item{alpha}{
Quantile to be estimated. Default is find the median.
}
\item{wt}{
Weight vector default is constant values. Passing nonconstant weights is a
pretty strange thing to do. 
}
\item{xgrid}{
Grid of x values to evaluate the estimated quantile function. Default 
is the unique sorted values of x.
}
\item{derivative}{
Specifies whether the function itself of derivatives should be evaluated 
at xgrid. 
}
\item{hmin}{
Minimum value of log( lambda) used for GCV grid search.
}
\item{hmax}{
Maximum value of log( lambda) used for GCV grid search.
}}
\value{
Object of class qspreg with many arguments similar to sreg object. 
}
\description{
This is an experimental function to find the smoothing parameter for a 
quantile spline using a more appropriate criterion than mean squared 
error prediction. 
The quantile spline is found by an iterative algorithm using weighted 
least squares cubic splines. At convergence the estimate will also be a 
weighted natural  cubic spline but the weights will depend on the estimate.
Of course these weights are crafted so that the resulting spline is an 
estimate of the alpha quantile instead of the mean. 


 
}
\references{
Nychka,D. Oconnell, M. (1996)  "
}
\seealso{
sreg 
}
\examples{


}
\keyword{quantile}
% Converted by Sd2Rd version 0.2-a3.
