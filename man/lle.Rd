\name{lle}
\title{
Calculates local Lyapunov exponents for plotting.
}
\usage{
lle(jac, model=1, nprod=c(5, 10, 20, 40, 80), skip, statevector=F, 
lags=NA)
}
\arguments{
\item{jac}{
Jacobian matrix or a nnreg fit.
}
\item{model}{
Model number of fit used to calculate Jacobians.
}
\item{nprod}{
Vector of LLE products of Jacobians.
}
\item{skip}{
Columns of Jacobian matrix to skip in calculating LLEs.
For example, skip the columns associated with forcing functions.
}
\item{statevector}{
If false, a time-delay reconstruction model is assumed and a Jacobian matrix
n by d is expected, where n is the length of the time series and d is the
dimension of the state space. 
If true, a state space vector model is assumed and a Jacobian matrix n by d^2 
is expected.
}
\item{lags}{
Lagged time values used in the Jacobian matrix.
}}
\value{


\item{local}{
Matrix of LLEs with columns corresponding to the LLEs of the nprod values.
}
\item{nprod}{
Vector of LLE products of Jacobians.
}
\item{glb}{
Global Lyapunov exponent.
}
\item{model}{
Model number used to calculate Jacobians.
}}
\references{
S. Ellner, D.W. Nychka, and A.R. Gallant. 1992.  LENNS,  a
program  to  estimate  the  dominant  Lyapunov exponent of
noisy nonlinear systems from time series  data.  Institute
of  Statistics  Mimeo Series #2235, Statistics Department,
North Carolina State University, Raleigh, NC 27695-8203.
}
\seealso{
make.lle
}
\examples{
make.lags(rossler.state[1:200,1],c(1,2,3)) -> data.r  # create
# 3-d time delay vector model of the x variable of rossler system.
nnreg(data.r$x,data.r$y,5,5) -> fit # fit time series model using nnreg.
jac<- predict(fit, derivative=1) 
lle(jac) -> rossler.lle  # LLEs of Rossler data
summary(lle)
plot(rossler.lle)  # plot LLEs

# here is an easier way
nlar( rossler[1:200], lags=1:3, method="nnreg", k1=5)-> ou
lle( out) -> rossler.lle

}
\keyword{FUNFITS}
% Converted by Sd2Rd version 0.2-a3.
