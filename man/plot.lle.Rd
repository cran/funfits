\name{plot.lle}
\title{
Plots local Lyapunov exponents
}
\usage{
plot.lle(out)
}
\arguments{
\item{out}{
An object from the lle function.
}}
\value{
Plots a series of boxplots of LLEs  calculated
from the vector of nprod specified in the lle function. The
horizontal line is the global LE.
}
\references{
S. Ellner, D.W. Nychka, and A.R. Gallant. 1992.  LENNS,  a
program  to  estimate  the  dominant  Lyapunov exponent of
noisy nonlinear systems from time series  data.  Institute
of  Statistics  Mimeo Series #2235, Statistics Department,
North Carolina State University, Raleigh, NC 27695-8203.
}
\seealso{
make.lle, lle
}
\examples{
make.lags(rossler$state.vector[1:200,1],c(1,2,3)) -> data  # create
# 3-d time delay vector model of the x variable of rossler system.
nnreg(data$x,data$y,5,5) -> fit # fit time series model using nnreg.
lle(fit) -> rossler.lle  # LLEs of Rossler data
plot(rossler.lle)  # plot LLEs
}
\keyword{FUNFITS}
% Converted by Sd2Rd version 0.2-a3.
