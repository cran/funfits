\name{nlar}
\title{
Nonlinear autoregressive model
}
\usage{
nlar(Y, lags, cov=NA, method="nnreg", \dots)


}
\arguments{
\item{Y}{
The time series
}
\item{lags}{
A vector that specifies which lags of Y to use in  the autoregressive function
}
\item{cov}{
A vector or matrix of covariates as long as the Y series these are
additional variables that will be used in the regression function
}
\item{method}{
Name of S function to fit the nonparametric model e.g. nnreg tps
addreg
}
\item{...}{
Optional argument that as passed through to the regression method
}}
\value{
An object of class nlar


}
\description{
his function fits a model of the form:
 Y_t =  f( Y_(t-l1),\dots{},Y_(t-ld),S_t) + e_t
Where e_t is assumed to mean zero, uncorrelated errors. Such a form is
useful for testing whether a system is chaotic.


}
\references{
FUNFITS manual
}
\seealso{
lle, predict.nlar
}
\examples{
# Fit the rossler series. A toy dynamical system that is chaotic
# Use a neural network with 4 hidden units based on lags 1, 2 and 3 of
the series. 
nlar( rossler,lags=c(1,2,3), method="nnreg",k1=4)-> out
summary(out)
plot( out)
lle( out) # calculate local and global Lyapunov exponents
}
\keyword{dynamical}
% Converted by Sd2Rd version 0.2-a3.
