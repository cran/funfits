\name{make.lle}
\title{
Calculates global and local Lyapunov exponents
}
\usage{
make.lle(jac, nprod, statevector=F, verbose=T)
}
\arguments{
\item{jac}{
Jacobian matrix, the value of the partial derivatives evaluated at each time.
}
\item{nprod}{
Number of products of Jacobians multiplied together to calculate LLE. 
If nprod is -1, the global LE is calculated.
}
\item{statevector}{
If false, a time-delay reconstruction model is assumed and a Jacobian matrix
n by d is expected, where n is the length of the time series and d is the
dimension of the state space. 
If true, a state space vector model is assumed and a Jacobian matrix n by d^2 
is expected.
}
\item{verbose}{
Statements are printed to the screen.
}}
\value{
Two column matrix whose columns correspond the the calculation using a
singular value decomposition of the Jacobian matrix products or the
magnitude of the (1,1) entry of the Jacobian matrix products. The rows of
the matrix are the values of the LLEs at the time or the global LE.
}
\description{
}
\references{
S. Ellner, D.W. Nychka, and A.R. Gallant. 1992.  LENNS,  a
program  to  estimate  the  dominant  Lyapunov exponent of
noisy nonlinear systems from time series  data.  Institute
of  Statistics  Mimeo Series #2235, Statistics Department,
North Carolina State University, Raleigh, NC 27695-8203.
}
\seealso{
lle
}
\examples{
make.lags(rossler.state[1:200,1],c(1,2,3)) -> data  # create
# 3-d time delay vector model of the x variable of rossler system.
nnreg(data$x,data$y,5,5) -> fit # fit time series model using nnreg.
predict(fit,fit$x,derivative=1) -> jac # calculate Jacobian matrix
make.lle(jac,-1) -> rossler.le  # LE of Rossler data
}
\keyword{FUNFITS}
% Converted by Sd2Rd version 0.2-a3.
