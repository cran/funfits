\name{as.surface}
\title{
Creates an object for plotting using contour, persp or image
}
\usage{
as.surface(grid.list, z, order.variables="xy")
}
\arguments{
\item{grid.list}{
A description of the grid used to evaluate the function. This can
either be in the form of a grid.list the grid of points produced
by make.surface.grid


}
\item{z}{
The value of the function evaluated at the gridded points.


}
\item{order.variables}{
Either "xy" or "yx" specifies how the x and y variables used to
evaluate the function are matched with the x and y grids in the surface
object.


}}
\value{
A list of class surface. This object is a modest generalization of the
list input format for the S functions contour, image or persp. 
Besides the x, y
and z components there are also components xlab, ylab, and main that are
used by the surface function to add axis labels and a title.  


}
\seealso{
grid.list, make.surface.grid, surface, contour
}
\examples{
# Make a perspective of the surface Z= X**2 -Y**2
# Do this by evaluating quadratic function on a 25 X 25 grid


 
grid.l<- list( X= seq( -2,2,,25), Y= seq( -2,2,,25))
make.surface.grid( grid.l)-> xg
# xg is a 625X2 matrix that has all pairs of X and Y grid values
xg[,1]**2 - xg[,2]**2 -> z
# now fold z in the matrix format needed for persp
as.surface( grid.l, z)-> out.p
persp( out.p)
# also try  plot( out.p) to see the default plot for a surface object


}
\keyword{surface,}
% Converted by Sd2Rd version 0.2-a3.
