\name{predict.surface}
\title{
Evaluates a function on a surface
}
\usage{
predict.surface(out, grid.list=NA, extrap=F, chull.mask)
}
\arguments{
\item{out}{
An object from fitting a function to data.
}
\item{grid.list}{
A list with as many components as variables describing the surface.
All components should have a single value except the two that give the
grid points for evaluation. If the matrix or data frame has column names, 
these must appear in the grid list.
}
\item{chull.mask}{
To evaluate on a convex hull, NA's are assigned to outside the convex hull
of the data
}
\item{extrap}{
Extraplolation beyond the range of the data.
}}
\value{
The usual list components for making contour and perspective plots along 
with the component containing the grid list.
}
\seealso{
tps, krig, nnreg, predict, persp, contour, interp, grid.list
}
\examples{
tps( BD[,1:4], BD$lnya)-> fit  # fit surface to data
list( KCl="x", MgCl2=mean(BD[,2]),KPO4="y", dNTP=mean(BD[,4]))-> grid
# make grid 
predict.surface(fit,grid) -> out.p  # evalute on a grid on two
# variables holding two fixed 
surface(out.p) # surface and contour plot 
}
\keyword{FUNFITS}
% Converted by Sd2Rd version 0.2-a3.












