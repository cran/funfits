\name{summary.tps}
\title{
Summary for thin plate spline regression
}
\usage{
summary.tps(x, \dots)
}
\arguments{
\item{x}{
A tps object.
}}
\value{
Gives a summary of the thin plate spline regression. The components
include the function call, number of observations, effective degrees 
of freedom, residual degrees of freedom, root mean squared error, 
R-squared and adjusted R-squared, log10(lambda), cost, GCV minimum, 
power and a summary of the residuals.  
}
\description{
This function is a method for the generic function summary for class tps.
}
\seealso{
tps, summary
}
\examples{
tps(ozone$x,ozone$y) -> fit  # tps fit
summary(fit)                 # summary of fit
}
\keyword{splines,}
% Converted by Sd2Rd version 0.2-a3.
