\name{make.Amatrix.tps}
\title{
Computes matrix that relates predicted values to the vector of
dependent variables (Y)
}
\usage{
make.Amatrix.tps(out, x0, lambda)
}
\arguments{
\item{out}{
The output object from using tps to fit a data set
}
\item{x0}{
The points where the function is to be predicted
}
\item{lambda}{
Value of the smoothing parameter. Default is value in the tps object.
}}
\value{
A matrix such that when multiplied times Y gives the predicted values
of the spline at the points x0


}
\references{
FUNFITS manual
}
\seealso{
tps, make.Amatrix, predict.se.tps
}
\examples{
# Compute the A matrix or "hat" matrix for a thin plate spline
# check that this gives the same predicted values


tps( ozone$x, ozone$y)-> tps.out
make.Amatrix( tps.out, ozone$x)-> A
A\%*\%ozone$y -> test.fitted.values
# now compare this to predict( tps.out) or tps.out$fitted.values
#                    they should be the same


}
\keyword{FUNFITS,}
% Converted by Sd2Rd version 0.2-a3.
