\name{nnregCI}
\title{
Finds a confidence set of parameters for a neural net fit. 
}
\usage{
nnregCI(fit, model=fit$best.model, ngrind=250, ntries=100, npol=20, 
clevel=0.95, cut1=NA, cut2=NA, nfits=500, tol1=1e-06, tol2=1e-09, 
itmax1=250, itmax2=10000, fdata, fout="nnci.out", seed)
}
\arguments{
\item{fit}{
A nnreg object.
}
\item{model}{
Model number used in finding joint parameter confidence set. Default is the 
best model based on GCV(2).
}
\item{ngrind}{
Number of coarse optimizations.
}
\item{ntries}{
Number of random starting values for each coarse optimization.
}
\item{npol}{
Number of coarse fits improved, i.e polish, using smaller minimization
tolerance.
}
\item{clevel}{
Confidence level used in finding joint parameter confidence set. Default is the
0.95 level.
}
\item{cut1}{
RMSE value corresponding to the clevel confidence level.
}
\item{cut2}{
RMSE value corresponding to 80\% of the RMSE value corresponding to the clevel
confidence level.
}
\item{nfits}{
Number of fits (parameter sets) found in the confidence set. Maximum is 500.
}
\item{tol1}{
Minimization tolerance for coarse optimizations.
}
\item{tol2}{
Minimization tolerance for polish optimizations.
}
\item{itmax1}{
Maximum number of iterations performed in the minimization routine for
coarse optimizations.
}
\item{itmax2}{
Maximum number of iterations performed in the minimization routine for
polish optimizations.
}
\item{fdata}{
Temporary UNIX file name for the data.
}
\item{fout}{
Temporary UNIX file name for the output.
}
\item{seed}{
Seed used in generating the random parameter starts.
}}
\value{
Object of class nnreg. The component model is a list of the parameters
for each fitted model. Each component model is of class netfit. 

\item{model}{
Component model of class netfit. Includes a list of the dimension of
the x matrix, the number of hidden units used in the model, the mean
of each column of the x matrix, the mean of the y values, the standard
deviation of each column of the x matrix, the standard deviation of
the y values, the number of parameters in the model and the parameters
of model. 
}
\item{summary}{
Partial Fortan program output. Summary of the nnreg fit. Includes a 
summary of the specified number of fitted values.
}
\item{call}{
Call to the function.
}
\item{x}{
Matrix of independent variables.
}
\item{y}{
Vector of dependent variables.
}
\item{n}{
Number of observations or length of y.
}
\item{nfits}{
Number of fits (parameter sets) found in the confidence set.
}
\item{seed}{
Seed used in generating the random parameter starts.
}}
\section{Side Effects}{
This function does the bulk of the computation using a stand-alone
FORTRAN program running in the UNIX shell. This operation is
transparent to the user. For large problems the input files can be
setup using this function and the fitting program can be run
separately in the background.
}
\description{
The joint parameter confidence set for a neural net fit is all the 
neural net parameter sets (theta) such that 
S(theta) <= S(theta^hat)*[1+(p/n-p)*F(p,n-p,alpha)],
where S(theta) is the residual sum of squares, theta^hat is the least- 
squares estimate of theta, p is the number of parameters of the model 
and n is the number of data points. For the F distribution, alpha is the
probability level.


The program finds parameter sets which satisfy the above inequality.
The value of cut1 is RMSE(theta^hat)*sqrt([1+(p/n-p)*F(p,n-p,alpha)]).
The value of cut2 is .8*cut1. Approximately 20\% of the fits will have a 
RMSE of cut1 and the remaining 80\% will be uniform between RMSE(theta^hat)
and cut1. This distribution of parameter sets is to make sure that the 
parameter sets cover the confidence region. The actual value of cut2 is 
used only as a check for the covering of the confidence region. The returned 
component summary has a count of the fits between cut1 and cut2 and also 
below cut2.


Parameters of the model are estimated by nonlinear least squares. The
parameter space has a large number of local minimum so the strategy is
to generate "many" parameter sets at random and iterate these starts
with a minimization algorithm. The two function parameters ntries and
ngrid are used in generating the many starting parameter sets for
nonlinear least squares. Ngrind is the number of cubes growing
geometrically over a range of magnitude of parameters. Ntries is the
number of parameter sets generated at random by a uniform distribution
in each cube. The best parameter set ( out the Ntries ) in each cube
is used as the start of a coarse optimization.  Npol of these coarse
fits are selected for further refinement by a minimization with
smaller tolerance.


The target RMS for a fit is generated as described above. The parameter 
sets for the confidence sets are generated in the polishing stage and in
groups of the optional argument npol. The file nnregCI.cut contains 
information about the polished fits. The 7th column is target RMSE value 
the 8th column is the difference between target RMSE and the root finder's 
RMSE. The 9th column is the value of cut1 and the 10th column is the value 
of cut2.  


}
\references{
B.A. Bailey, S. Ellner, D.W. Nychka. 1996. Chaos with Confidence: Asymptotics 
and Applications of Local Lyapunov Exponents. Proceedings on Nonlinear Dynamics
and Time Series, Building a Bridge Between the Natural and Statistical
Sciences. Fields Institute Communications.


S. Ellner, D.W. Nychka, and A.R. Gallant. 1992. LENNS, a program to estimate
the dominant Lyapunov exponent of noisy nonlinear systems from time series 
data. Institute of Statistics Mimeo Series #2235, Statistics Department,
North Carolina State University, Raleigh, NC 27695-8203.


D.W. Nychka, S. Ellner, D. McCaffrey, and A.R. Gallant. 1992. Finding Chaos
in Noisy Systems. J. R. Statist. Soc. B 54:399-426. 
}
\seealso{
predict.nnreg, predict.netfit, plot.nnreg, summary.nnreg, print.nnreg
}
\examples{
nnreg(ozone$x,ozone$y,1,2) -> fit # fitting a surface to ozone 
# measurements, from 1 to 2 hidden units


nnregCI(fit) -> fit.ci # finds 500 fits in the .95 confidence set based
# on the best model from the above fit


}
\keyword{neural}
% Converted by Sd2Rd version 0.2-a3.
