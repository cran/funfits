\name{make.Amatrix.krig}
\title{
Computes matrix that relates predicted values to the vector of
dependent variables (Y)
}
\usage{
make.Amatrix.krig(out, x0, lambda)
}
\arguments{
\item{out}{
The output object from using krig to fit a data set
}
\item{x0}{
The points where the function is to be predicted
}
\item{lambda}{
Value of the smoothing parameter. Default is value in the krig object.
}}
\value{
A matrix such that when multiplied times Y gives the predicted values
of the spline at the points x0


}
\references{
FUNFITS manual
}
\seealso{
krig, make.Amatrix, predict.se.krig
}
\examples{
# Compute the A matrix or "hat" matrix for a thin plate spline
# check that this gives the same predicted values


krig( ozone$x, ozone$y, exp.cov)-> krig.out
make.Amatrix( krig.out, ozone$x)-> A
A\%*\%ozone$y -> test.fitted.values
# now compare this to predict( krig.out) or krig.out$fitted.values
#                    they are all  the same!


}
\keyword{FUNFITS,}
% Converted by Sd2Rd version 0.2-a3.






