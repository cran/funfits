\name{plot.sreg}
\title{
Plots smoothing spline regression object
}
\usage{
plot.sreg(out, pch=".", main=NA)
}
\arguments{
\item{out}{
A sreg object
}
\item{pch}{
Ploting symbol for data points
}
\item{main}{
Title of the plot. Default is the function call.
}}
\description{
This function creates three summary plots of the sreg object. The first plot
is a plot of the sreg fit to the data. The second is predicted values vs
residuals, and the third is effective number of parameters vs estimated 
average squared prediction errors.
}
\seealso{
sreg, plot
}
\examples{
sreg(ozone$x[,1],ozone$y)-> fit # fit sreg to ozone at longitude values
plot(fit)                       # plot sreg fit
}
\keyword{smoothing}
% Converted by Sd2Rd version 0.2-a3.
