\name{print.krig}
\title{
Print kriging fit results.
}
\usage{
print.krig(obj)
}
\arguments{
\item{obj}{
Oject from krig fitting.
}}
\value{
Selected summary results from krig.
}
\seealso{
print, summary.krig
}
\examples{
krig(ozone$x,ozone$y,exp.cov) -> fit # krig fit
print(fit) # print the summary
fit # this will too
}
\keyword{splines,}
% Converted by Sd2Rd version 0.2-a3.
