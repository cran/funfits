\name{summary.nnreg}
\title{
Summary for neural net regression
}
\usage{
summary.nnreg(out, noprint=F)
}
\arguments{
\item{out}{
A nnreg object
}
\item{noprint}{
Print the summary
}}
\value{
Gives a summary of the neural net regression. The components
include the function call, number of hidden units, number of parameters,
degree of freedom of residuals, root mean squared error and GCV.
}
\description{
This function is a method for the generic function summary for class nnreg.
}
\seealso{
nnreg, summary
}
\examples{
nnreg(ozone$x,ozone$y,1,2) -> fit # fitting a surface to ozone 
# measurements, from 1 to 2 hidden units
summary(fit) # summary of fit
}
\keyword{neural}
% Converted by Sd2Rd version 0.2-a3.
