\name{nkreg}
\title{
Normal kernel regression estimate


}
\usage{
nkreg(data.x, data.y, bandwidth, n.points, grid)


}
\arguments{
\item{data.x}{
Vector or matrix of independent variables. 
}
\item{data.y}{
Vector of the dependent variable
}
\item{bandwidth}{
Vector of bandwidths
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
nkden, nkden.cv, ksmooth
}
\examples{
nkreg( auto.paint$thick, auto.paint$DOI, 5.0)-> look
plot( auto.paint$thick, auto.paint$DOI)
lines( look$x, look$y)
plot( look$y, auto.paint$DOI- look$y, xlab="Predicted", ylab="Residuals")


}
\keyword{kernel}
% Converted by Sd2Rd version 0.2-a3.
