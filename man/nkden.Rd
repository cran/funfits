\name{nkden}
\title{
Normal kernel density estimate
}
\usage{
nkden(data, bandwidth, n.points, grid)


}
\arguments{
\item{data}{
A vector or matrix of oberservations. Rows are considerd to be
independent random samples from a continous distribution. 
}
\item{bandwidth}{
The bandwidht for the kernels estimates in the scale of standard
deviation for the normal density.
}
\item{n.points}{
Number of equally spaced points to evaluate a univariate density
estimate. 
}
\item{grid}{
A vector or matrix of values to evaluate the estimate.The defautl is
to use the data.
}}
\value{
A list where x is the points used for evaluation, y the density
estimates at these points and bandwidths and h the vector of
bandwidths. If more than one bandwidth is given then the estimates are
arranged as columns in the matrix y. 
}
\seealso{
nkreg, nkden.cv, ksmooth
}
\examples{


# univariate estimate with several bandwidths


nkden( minitri$swim, c(2.0,4.0,6.0),n.points=150)-> look 
matplot( look$x, look$y, type="l") # plot all of them togther




# a bivariate estimate


nkden( minitri[,1:2], 2.0) -> look2
# create gridded surface from values at data points.


interp( look2$x[,1], look2$x[,2], look2$y) -> surface
persp(surface, xlab="swim", ylab="bike")


}
\keyword{kernel}
% Converted by Sd2Rd version 0.2-a3.
