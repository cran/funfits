\name{surface.surface}
\title{
Plots a surface
}
\usage{
surface.surface(obj, lab=NA, type="b", zlab, xlab, ylab, \dots)
}
\arguments{
\item{obj}{
An object from a fit, for example a tps, krig, or nnreg object.
Object can also be of class surface from a predict.surface function. 
}
\item{lab}{
Axes labels
}
\item{type}{
A type="b" for both a perspective and contour plot, type "p"
for a perspective plot and type="c" for a contour plot. 
}
\item{zlab}{
z-axes label
}
\item{xlab}{
x-axes label
}
\item{ylab}{
y-axes labels
}
\item{...}{
Other graphical parameters
}}
\description{
A surface plot is generated.
}
\seealso{
surface, surface.tps, predict.surface
}
\examples{
tps( BD[,1:4], BD$lnya)-> fit  # fit surface to data
list( KCl="x", MgCl2=mean(BD[,2]),KPO4="y", dNTP=mean(BD[,4]))-> grid
# make grid 
predict.surface(fit,grid) -> out.p  # evalute on a grid on two
# variables holding two fixed 
surface(out.p) # surface plot and contour plot 
}
\keyword{FUNFITS}
% Converted by Sd2Rd version 0.2-a3.






