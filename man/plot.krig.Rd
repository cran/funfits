\name{plot.krig}
\title{
Diagnoistic and summary plots of the kriging object
}
\usage{
plot.krig(out, main=NA, digits=4, graphics.reset=T, \dots)
}
\arguments{
\item{out}{
A krig object.
}
\item{main}{
Title of the plot. Default is the function call.
}
\item{digits}{
Number of significant digits for the RMSE label.
}
\item{graphics.reset}{
Reset to original graphics parameters after function plotting.
}}
\description{
This function creates four summary plots of the krig object. The first
plot is a plot of the krig fit to the data. The second is predicted
values vs residuals. The third is effective number of parameters
vs GCV and the fourth is lambda vs GCV.
}
\seealso{
krig, summary.krig
}
\examples{
krig(ozone$x, ozone$y,exp.cov) -> fit # fitting a surface to ozone 
# measurements
plot(fit) # plots fit and residual plots
}
\keyword{splines,}
% Converted by Sd2Rd version 0.2-a3.
