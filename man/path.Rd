\name{path}
\title{
Finds path of steepest ascent
}
\usage{
path(fit, r, x.center, x.scale, maximize=T, scale=0.0001, start, trace.z=F)
}
\arguments{
\item{fit}{
Object from a fitting procedure. A predict function for the fit must be
available.
}
\item{r}{
Vector of increasing values for the radii of the hyperspheres.
}
\item{x.center}{
Center of the hypersheres, in the original coordinates of the data.
}
\item{x.scale}{
Scale of the hyperspheres, in the original coordinates of the data.
}
\item{maximize}{
Default is to search for the maximum.
}
\item{scale}{
Transformation scale of data.
}
\item{start}{
Starting value for the search, the default is the middle of the region.
}
\item{trace.z}{
}}
\value{
Array with the components of the values of the data, the radius of the
hypersphere of the search and the predicted value of the surface.
}
\description{
The object from a fitting procedure must have a predict function.
Path calls the function nlminb which finds a local minimum of a 
smooth nonlinear function subject to bounded-constrained paramenters.
}
\seealso{
optim, nlminb, predict.tps, predict.krig, predict.nnreg
}
\examples{
tps(BD[,1:4],BD$lnya,scale.type="range") -> fit # fitting a DNA strand
# displacement amplification surface to various buffer compositions
path(fit,c(.01,.05)) -> fit.path # path from center of data to
# hyperspheres with radii .01 and .05
}
\keyword{FUNFITS}
% Converted by Sd2Rd version 0.2-a3.
