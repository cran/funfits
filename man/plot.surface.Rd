\name{plot.surface}
\title{
Plots a surface
}
\usage{
plot.surface(obj, lab=NA, type="p", zlab, xlab, ylab, graphics.reset = T, \dots)
}
\arguments{
\item{obj}{
A surface object at the minimum a list with components x,y and z
in the same form as are contour and perp functions.
Also an object from a fit, for example a tps, krig, or nnreg object.
}
\item{lab}{
Axes labels
}
\item{type}{
A type="p" for a perspective plot, type="c" for a contour plot 
and type="b" for both.
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
\item{graphics.reset}{
Reset to original graphics parameters after function plotting.
}
\item{...}{
Other graphical parameters
}}
\description{
A surface plot is generated.
}
\seealso{
surface, predict.surface, as.surface
}
\examples{
tps( BD[,1:4], BD$lnya)-> fit  # fit surface to data
list( KCl="x", MgCl2=mean(BD[,2]),KPO4="y", dNTP=mean(BD[,4]))-> grid
# make grid list
predict.surface(fit,grid) -> out.p  # evalute on a grid on two
# variables holding two fixed at their mean levels
plot.surface(out.p) # surface plot 
}
\keyword{FUNFITS}
% Converted by Sd2Rd version 0.2-a3.






