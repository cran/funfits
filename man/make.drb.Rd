\name{make.drb}
\title{
Evaluates the Demmler-Reinsch basis fucntions 
}
\usage{
make.drb(out, x1=out$x)
}
\arguments{
\item{out}{
tps output object 
}
\item{x1}{
matrixc of points to evaluate the basis functions. 
}}
\value{
A matrx where the columns index different basis fucntions and rows
index the points in x1.


}
\references{
See the FUNFITS manual for details on teh relationship between this
basis and the usual one. 


}
\seealso{
tps
}
\examples{
# First create a tps object
tps( ozone$x, ozone$y)-> out
make.drb( out)-> hold
# to calculate the u vector important for tps calulations
 t(hold)\%*\% out$y -> u


}
\keyword{FUNFITS,}
% Converted by Sd2Rd version 0.2-a3.


