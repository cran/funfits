\name{predict.se.tps}
\title{
Standard errors of predictions for thin plate spline regression
}
\usage{
predict.se.tps(out, x)
}
\arguments{
\item{out}{
A fitted tps object.
}
\item{x}{
Matrix of x values on which to calculate the standard errors of predictions
of the thin plate spline regression. If omitted, the out$x will be used.
}}
\value{
A vector of standard errors for the predicted values of the thin plate
spline regression.
}
\description{
The predictions are represented as linear combinations of the
dependent variables using the function make.Amatrix. Given the
estimate for sigma2 and the assumption of uncorrelated errors, the
prediction variance associated with the estimate at a particular point 
is sigma2 multiplied by the sum of the squared
linear weights. Note that no adjustment is made for the bias in the
spline. This can be substantial and so these standard errors only
represent part of the uncertainty in the estimate. For more details
refer tot he FUNFITS manual.  


 
}
\seealso{
tps, predict.tps
}
\examples{
tps(ozone$x,ozone$y) -> fit                 # tps fit
predict.se.tps(fit)                         # std errors of predictions
cbind(seq(87,89,,10),seq(40,42,,10)) -> x   # new x matrix
predict.se.tps(fit,x) -> out                # std errors of predictions


tps(BD[,1:4],BD$lnya,scale.type="range") -> fit # fitting a DNA strand
# displacement surface to various buffer compositions
predict.se.tps(fit) -> out # std erros of predictions
}
\keyword{splines,}
% Converted by Sd2Rd version 0.2-a3.
