\name{surface.addreg}
\title{
Plots a surface and contours
}
\usage{
surface.addreg(obj, grid.list=NA, extrap=F, graphics.reset=T, \dots)
}
\arguments{
\item{obj}{
A addreg object.
}
\item{grid.list}{
A list with as many components as variables describing the surface.
All components should have a single value except the two that give the
grid points for evaluation. If the matrix or data frame has column names, 
these must appear in the grid list.
}
\item{extrap}{
Extraplolation beyond the range of the data.
}
\item{graphics.reset}{
Reset to original graphics parameters after function plotting.
}
\item{...}{
Any plotting arguments.
}}
\value{

}
\description{
Two plots are made, a surface plot and a contour plot.
}
\seealso{
addreg, surface
}
\examples{
# Additive model fit to the minitriathon data. A prediction of
# run times based on the swim and bike performance.
# The smoothing parameters found by GCV with cost=2.
cbind( minitri$swim, minitri$bike)-> x
addreg( x,minitri$run,cost=2)-> out
surface(out) # plots surface and contours of fit
}
\keyword{splines,}
% Converted by Sd2Rd version 0.2-a3.
