\name{surface.tps}
\title{
Plots a surface and contours
}
\usage{
surface.tps(obj, grid.list=NA, extrap=F, graphics.reset=T, \dots)
}
\arguments{
\item{obj}{
A tps object.
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
Any ploting options.
}}
\description{
Two plots are made, a surface plot and a contour plot.
}
\seealso{
tps, surface
}
\examples{
tps(ozone$x,ozone$y) -> fit # tps fit
surface(fit) #plots surface and contours of tps fit


tps(BD[,1:4],BD$lnya,scale.type="range") -> fit # fitting a DNA strand
# displacement amplification surface to various buffer compositions
surface(fit,grid.list=list(MgCl2='x',dNTP='y',KCl=35,KPO4=20))) # plots 
# surface on MgCl2-dNTP grid at specific KCL and KPO4 levels
}
\keyword{splines,}
% Converted by Sd2Rd version 0.2-a3.
