\name{surface.krig}
\title{
Plots a surface and contours
}
\usage{
surface.krig(obj, grid.list=NA, extrap=F, graphics.reset=T, \dots)
}
\arguments{
\item{obj}{
A krig object.
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
Any plotting options.
}}
\description{
Two plots are made, a surface plot and a contour plot.
}
\seealso{
krig, surface
}
\examples{
krig(ozone$x,ozone$y) -> fit # krig fit
surface(fit) #plots surface and contours of krig fit
}
\keyword{splines,}
% Converted by Sd2Rd version 0.2-a3.
