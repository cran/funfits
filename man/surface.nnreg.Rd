s.BG
\name{surface.nnreg}
\title{
Plots a surface and contours
}
\usage{
surface.nnreg(obj, grid.list=NA, extrap=F, \dots)
}
\arguments{
\item{obj}{
A nnreg object
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
\item{...}{
Any ploting options.
}}
\description{
Two plots are made, a surface plot and a contour plot.
}
\seealso{
nnreg, surface
}
\examples{
nnreg(ozone$x,ozone$y,1,3) -> fit # fitting a surface to ozone
surface(fit) # plots surface and contours of nnreg fit


nnreg(as.matrix(BD[,1:4]),BD$lnya,1,1) -> fit # fitting a DNA strand
# displacement amplification surface to various buffer compositions
surface(fit,grid.list=list(MgCl2='x',dNTP='y',KCl=35,KPO4=20)) # plots 
# surface on MgCl2-dNTP grid at specific KCL and KPO4 levels


}
\keyword{neural}
% Converted by Sd2Rd version 0.2-a3.
