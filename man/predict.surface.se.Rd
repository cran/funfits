 .BG
\name{predict.surface.se}
\title{
Standard error predictions of surface
}
\usage{
predict.surface.se(out, grid.list, chull.mask, extrap=F, \dots)
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
Exaplolation beyond the range of the data.
}
\item{...}{
Other plotting parameters
}}
\value{
The usual list for making contour and perspective plots along with the
component containing the grid list.
}
\seealso{
predict.surface, tps, krig, nnreg, persp, contour, grid.list
}
\examples{
tps( BD[,1:4], BD$lnya)-> fit  # fit surface to data
list( KCl="x", MgCl2=mean(BD[,2]),KPO4="y", dNTP=mean(BD[,4]))-> grid
# make grid 
predict.surface.se(fit,grid) -> out.p  # evalute standard erro surface
# on a grid
surface(out.p) # perspective and contour plot of standard error surface 
}
\keyword{standard}
% Converted by Sd2Rd version 0.2-a3.
