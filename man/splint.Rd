\name{splint}
\title{
Cubic spline interpolation
}
\usage{
splint(x, y, xgrid, derivative=0)
}
\arguments{
\item{x}{
The x values the define the curve or a two column matrix of 
x and y values.
}
\item{y}{
The y values that are paired with the x's.
}
\item{xgrid}{
The grid to evaluate the fitted cubic interpolating curve.
}
\item{derivative}{
Indicates whether the function or a a first or second derivative
should be evaluated.
}}
\value{
A vector consisting of the spline evaluated at the grid values.
}
\description{
Fits a piecewise cubic polynomial to the x and y values assuming that
the second and third derivatives are zero at the range of the x
values. For this reason extrapolation outside the range of the x
values will be a linear function.
}
\references{
See Additive Models by Hastie and Tibshriani.
}
\seealso{
spreg, sreg, tps 
}
\examples{
splint(ozone$x[,1],ozone$y,seq(min(ozone$x[,1]),max(ozone$x[,1]),,40)) -> fit 
# evaluate cubic spline to ozone at longitude values on grid of 
# longitude values
}
\keyword{splines,}
% Converted by Sd2Rd version 0.2-a3.


