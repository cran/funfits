\name{surface}
\title{
Plots a surface
}
\usage{
surface(x, gridlist=NA, extrap=F, \dots)
}
\arguments{
\item{x}{
An Splus object
}
\item{gridlist}{
A list with as many components as variables describing the surface.
All components should have a single value except the two that give the
grid points for evaluation. If the matrix or data frame has column names, 
these must appear in the grid list.
}
\item{extrap}{
Extraplolation beyond the range of the data.
}
\item{...}{
Method for specific class of data  may have additional arguments.
}}
\value{

}
\seealso{
surface.tps, surface.krig, surface.nnreg, plot.surface, grid.list
}
\examples{
tps(ozone$x,ozone$y) -> fit #tps fit to ozone data
surface(fit) #plots the surface fit to the ozone data
}
\keyword{FUNFITS}
% Converted by Sd2Rd version 0.2-a3.
