\name{print.tps}
\title{
Print thin plate spline fit results
}
\usage{
print.tps(x, \dots)
}
\arguments{
\item{x}{
Object from tps fitting.
}}
\value{
Selected summary results from tps.
}
\seealso{
print, summary.tps
}
\examples{
tps(ozone$x,ozone$y) -> fit # tps fit
print(fit) # print a summary
fit # this will too
}
\keyword{splines,}
% Converted by Sd2Rd version 0.2-a3.
