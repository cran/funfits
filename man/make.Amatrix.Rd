\name{make.Amatrix}
\title{
Matrix relating predicted values to the dependent (Y) values
}
\usage{
make.Amatrix(object, \dots)
Output object from fitting a data set using a FUNFITS regression method.
Currently this is supported only for the tps and krig functions.
}
\arguments{
\item{object}{


}
\item{...}{
Additional arguments that indicate the value of the smoothing
parameter to use and the X values where predictions should be made. 
}}
\value{
A matrix where the number of rows is equal to the number of predicted points
and the number of columns is equal to the length of the Y vector. 
}
\description{
For the tps and krig functions the A matrix is constructed based on
the representation of the estimate as a generalized ridge regression.
The matrix expressions are explained in the FUNFITS manual. For linear
regression the matrix that gives predicted values is often referred to
as the "hat" matrix and is useful for regression diagnostics.  For
smoothing problems the effective number of parameters in the fit is
usually taken to be the trace of the A matrix. Note that while the A
matrix is usually constructed to predict the estimated curve at the
data point this S function  does not have such restrictions. This is
possible because any value of the estimated curve will be a linear
function of Y. 
}
\references{
FUNFITS manual
}
\seealso{
predict.se.krig, predict.se.tps
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
\keyword{FUNFITS}
% Converted by Sd2Rd version 0.2-a3.
