re.BG
\name{predict.krig}
\title{
Evaluation of kriging surface
}
\usage{
predict.krig(out, x, lambda=NA, model=NA)
}
\arguments{
\item{out}{
Fit object from the kriging function.
}
\item{x}{
Matrix of x values on which to evaluate the kriging surface.
If omitted, the data x values, i.e. out$x  will be used.
}
\item{lambda}{
Smoothing parameter. If omitted, out$lambda will be used.
}
\item{model}{
Generic argument that may be used to pass a different lambda. 
}}
\value{
Vector of predicted responses.
}
\seealso{
krig, predict.surface
}
\examples{
krig(ozone$x,ozone$y,exp.cov) -> fit        # krig fit
cbind(seq(87,89,,10),seq(40,42,,10)) -> x   # new x matrix
predict(fit,x) -> out                       # evaluate fit at x
}
\keyword{splines,}
% Converted by Sd2Rd version 0.2-a3.
