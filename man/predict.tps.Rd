\name{predict.tps}
\title{
Evaluation of thin plate spline surface.
}
\usage{
predict.tps(out, x, lambda, derivative=0, model=NA)
}
\arguments{
\item{out}{
Fitted tps object.
}
\item{x}{
Matrix of x values on which to evaluate the thin plate spline surface.
If omitted, out$x will be used.
}
\item{lambda}{
Smoothing parameter. If omitted, out$lambda will be used.
}
\item{derivative}{
Derivative of function is returned if derivative=1.
Generic argument that may be used to pass a different lambda.
}}
\value{
Vector of predicted responses. If derivative=1 a vector of
derivatives or a matrix of partial derivatives is returned.
}
\seealso{
tps, predict.surface
}
\examples{
tps(ozone$x,ozone$y) -> fit                 # tps fit
cbind(seq(87,89,,10),seq(40,42,,10)) -> x   # new x matrix
predict(fit,x) -> out                       # evaluate fit at x
}
\keyword{splines,}
% Converted by Sd2Rd version 0.2-a3.
