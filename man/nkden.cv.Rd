\name{nkden.cv}
\title{
Cross- validation function for a normal kernel density estimate
}
\usage{
nkden.cv(data, h)
}
\arguments{
\item{data}{
A vector or matrix of observations. It is assumed that the rows are independent
}
\item{h}{
A vector of bandwidths to evaluate the Cross-validation function. The
default is to search on a grid scaled tot he ragne of the data. 
}}
\value{
A list with the bandwidth in h, the CV function in CV.f and the value of
h at the smallest value of the CV function in h.CV.
}
\description{
This function can be used to pick a data-based estimate of the
bandwidth for a normal kernel density estimate. Values of h that make
the CV function small are to be prefered. In some cases the CV
function should be ignored beacuse it has a minimum at h=0. The CV
function should always be plotted before the minimizing bandwidht is
used in a density estimate. 


}
\examples{
nkden.cv( minitri$swim)-> look
plot( look$CV.f , look$h, type="l", xlab="h", ylab=" LS CV function")
nkden( minitri$swim, look$h.CV, n.points=150)-> out
plot( out$x, out$y, type="l", xlab=" swim times", " ylab=" pdf")


}
\keyword{FUNFITS}
% Converted by Sd2Rd version 0.2-a3.
