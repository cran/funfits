\name{ozone}
\title{
Data set of ozone measurements at 20 Chicago monitoring stations.
}
\section{SUMMARY}{
The \code{ozone} data is a list of two components, x and y. The x
component is longitude and latitude position of each of the 20
Chicago monitoring stations. The y component is the average
daily ozone values over the time period 6/3/87-8/30/87.
}
\section{"DATA}{
This data set is a list containing the following components:
}
\arguments{
\item{x}{
Longitude-latitude postion of monitoring stations.
}
\item{y}{
Average daily ozone values over 1987 summer.
}}
\section{SOURCE}{
AIR, the EPA air quality data base.
}
\seealso{
tps, krig, nnreg
}
\examples{
tps(ozone$x, ozone$y) -> fit # fitting a surface to ozone measurements.
plot(fit) # plots fit and residuals
}
\keyword{FUNFITS}
% Converted by Sd2Rd version 0.2-a3.
