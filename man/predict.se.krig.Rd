\name{predict.se.krig}
\title{
Standard errors of predictions for kriging estimate
}
\usage{
predict.se.krig(out, x, se=F, cov.fun, rho, sigma2, stationary=T)
}
\arguments{
\item{out}{
A fitted krig object.
}
\item{x}{
Matrix of x values on which to calculate the standard errors of predictions
of the thin plate spline regression. If omitted, the out$x will be used.
}
\item{se}{
~Describe se here
}
\item{cov.fun}{
Covariance function for the random surface. If omitted then the
function from the krig object ( out$cov.function) is used. 
}
\item{rho}{
Parameter that multiplies cov.function. If omitted the value out$rho, estimated
from the data is used. 
}
\item{sigma2}{
Variance of the measurement error. If omitted the value out$sigma2, estimated
from the data is used
}
\item{stationary}{
If true the covariance function is assumed to be stationary and more
efficient computations can be made. 
}}
\value{
A vector of standard errors for the predicted values of the kriging fit.
}
\description{
The predictions are represented as a linear combination of the
dependent variable, Y. Based on this representation the conditional
variance is the same as the expected value of (P(x) + Z(X) - LY)^2.
where P(x)+Z(x) is the value of the surface at x and LY is the linear
combination that estimates this point. Finding this expected value is
straight forward given the unbiasedness of LY for P(x) and the
covariance for Z and Y. In these calculations is assumed that the
covariance parameters are fixed. This is an approximation since in
most cases they have been estimated from the data. 
Note that the linear commbination is based on the covariance function
from the krig object. The covariance for the random surface may be
this same fucntion or another if it is explicitly passed as an
argument.
See the FUNFITS manual for more details. 
}
\seealso{
krig, predict.krig
}
\examples{
krig(ozone$x,ozone$y,exp.cov) -> fit        # krig fit
predict.se.krig(fit)                        # std errors of predictions
cbind(seq(87,89,,10),seq(40,42,,10)) -> x   # new x matrix
predict.se.krig(fit,x) -> out               # std errors of predictions
}
\keyword{splines,}
% Converted by Sd2Rd version 0.2-a3.
