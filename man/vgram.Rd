\name{vgram}
\title{
Variogram
}
\usage{
vgram(loc, y, lon.lat=F)
}
\arguments{
\item{loc}{
Matrix where each row is the coordinates of an observed point of the field
}
\item{y}{
Value of the field at locations
}
\item{lon.lat}{
If true, locations are assumed to be longitudes and latitudes and
distances found are great circle distances. Default is false.
}}
\value{
Object of class vgram.  This is list with components d the pairwise
distances
, vgram the value of the variogram and loc the locations.
}
\references{
See any standard reference on spatial statistics. For example
Cressie Spatial Statistics 
}
\examples{
# compute variogram for the Chicago ozone field
vgram( ozone$x, ozone$y, lon.lat=T)-> out
plot( out$d, out$vgram, xlab="Distance", ylab="Squared differences")


}
\keyword{spatial}
% Converted by Sd2Rd version 0.2-a3.






