\name{plot.tps}
\title{
Diagnoistic and summary plots for a thin plate spline regression object
}
\usage{
plot.tps(out, main=NA, digits=4, graphics.reset=T, \dots)
}
\arguments{
\item{out}{
A tps object.
}
\item{main}{
Title of the plot. Default is the function call.
}
\item{digits}{
Number of significant digits for the RMSE label.
}
\item{graphics.reset}{
Reset to original graphics parameters after function plotting.
}
\item{...}{
Any plotting arguments.
}}
\description{
This function creates four summary plots of the tps object. The first
plot is a plot of the tps fit to the data. The second is predicted
values vs residuals. The third is effective number of parameters
vs GCV and the fourth is lambda vs GCV.
}
\seealso{
tps, summary.tps
}
\examples{
tps(ozone$x, ozone$y) -> fit # fitting a surface to ozone measurements
plot(fit) # plots fit and residual plots


tps(BD[,1:4],BD$lnya,scale.type="range") -> fit # fitting a DNA strand
# displacement amplification surface to various buffer compositions
plot(fit) # plots fit and residual plots
}
\keyword{splines,}
% Converted by Sd2Rd version 0.2-a3.
