\name{persp.init,}
\title{
Sets up axes for perspective plotting using other functions from persp 
family.
}
\usage{
persp.init(xlim=c(-1, 1), ylim=c(-1, 1), zlim=c(-1, 1), eye=c(-6, -8, 5), \dots)
}
\arguments{
\item{xlim}{
Range of x coordinate.
}
\item{ylim}{
Range of y coordinate.
}
\item{zlim}{
Range of z coordinate.
}
\item{eye}{
Veiwing location. 
}
\item{...}{
Other graphical parameters that are passed on to persp. 
}}
\section{Side Effects}{
A data object is created for this session that has the information for 
creating the projections for subsequent plotting. 


  
}
\description{
The idea is to set up the plotting region for a given perspective and 
then use the other persp graphics primitives to build up a plot. 
This implementation is awkward in that SPLUS only sets up the perspective
information if a call to persp is made. persp.init must do this but it 
uses a trivial surface that will  match the axes. 
One can also use these functions to add to the usual SPLUS wire frame 
plot if the projection information it saved. An example of this is given 
below. 


}
\seealso{
perspp, persp, surface
}
\examples{


# add text and a symbol to a wire frame plot
# 


grid<- make.surface.grid( list( seq(-1,1,,15), seq(-1,1,,15)))
z<- grid[,1]**2 - grid[,2]**2   # saddle function
as.surface( grid, z)-> out.p


persp( out.p)-> proj.info   # save projection info


persp.points( 0,0,0, pch="+",cex=1.5,  persp.out=proj.info)
persp.text( 0,0,0, "   location of saddle point", cex=2,adj=0, 
persp.out=proj.info)






# plot a helix


theta<- seq( 0,50,,200)
x<- sin( theta)
y<- cos(theta) 
z<- theta/25




persp.init( range(x), range( y), range(z))
persp.lines( x,y,z)






#
# Draw a sphere in perspective. 
#
        u <- seq(0, 2 * pi,  , 30)
        v <- seq( - pi, pi,  , 30)      #
# make make up a parametric function for the sphere
#
function(uv)
{
        temp <- cbind(sin(uv[, 1]) * cos(uv[, 2]), cos(uv[, 1]) * 
cos(uv[, 2]), 
                sin(uv[, 2]) )
        temp
}


        out <- as.tiles(list(u, v), fun.ball)        #
# draw 3-d axes and setup perspective
#
        persp.init(xlim = c(-1, 1), ylim = c(-1, 1), zlim = c(0, 1))    #
# add tiles drawing from back to front
        persp.tiles(out)


 
}
\keyword{funfits}
