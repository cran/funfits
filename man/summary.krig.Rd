\name{summary.krig}
\title{
Summary for krigging fit
}
\usage{
summary.krig(x, \dots)
}
\arguments{
\item{x}{
A krig object.
}}
\value{
Gives a summary of the krigging oject. The components
include the function call, number of observations, effective degrees 
of freedom, residual degrees of freedom, root mean squared error, 
R-squared and adjusted R-squared, log10(lambda), cost, GCV minimum and 
a summary of the residuals. 
}
\description{
This function is a method for the generic function summary for class krig.
}
\seealso{
krig, summary
}
\examples{
krig(ozone$x, ozone$y, exp.cov) -> fit  # krig fit
summary(fit)                            # summary of fit
}
\keyword{splines,}
% Converted by Sd2Rd version 0.2-a3.
