\name{rtps}
\title{
Robust thin plate spline regression
}
\usage{
rtps(\dots, max.iter=20, acc=10 * .Machine$single.eps^0.5, finish=F, 
conv.method="lambda")
}
\arguments{
\item{...}{
Arguments for tps function. Matrix of independent variables x and vector of
dependent variables Y are required.
}
\item{max.iter}{
Maximum number of iterations to convergence.
}
\item{acc}{
Accuracy of convergence criteria.
}
\item{finish}{
Finish the iterations with two iterations of the bisquare.
}
\item{conv.method}{
Method of convergence, lambda or weights.
}}
\value{
tps object.
}
\references{
See Nonparametric Regression and Generalized Linear Models
by Green and Silverman.


See Additive Models by Hastie and Tibshirani.
}
\seealso{
tps, summary.tps
}
\examples{
#2-d example
rtps(ozone$x, ozone$y) -> fit # fitting a surface to ozone measurements,
using iterative reweighted least squares.
plot(fit) # plots fit and residuals
}
\keyword{splines,}
% Converted by Sd2Rd version 0.2-a3.
