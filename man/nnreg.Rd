\name{nnreg}
\title{
Fits a surface based on a neural network
}
\usage{
nnreg(x, y, k1, k2, start, ngrind=250, ntries=100, npol=20,
glow=-1.26, ghigh=1.26, scale=0.5, fdata, derivative=F,
fout="nnreg.out",run=T, just.setup=F, just.read=F,
fitted.values=F, all.fits=F, greedy=F, seed)
}
\arguments{
\item{x}{
Matrix of independent variables.
}
\item{y}{
Vector of dependent variables.
}
\item{k1}{
Lower limit for K, where K is the number of hidden units.
}
\item{k2}{
Upper limit for K.
}
\item{start}{
Starting values for parameters.
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
\item{glow}{
Lower limit for grid of initial parameter values.
}
\item{ghigh}{
Upper limit for grid of initial parameter values.
}
\item{scale}{
Scale factor for grid of initial parameter values.
}
\item{fdata}{
Temporary UNIX file name for the data.
}
\item{derivative}{
Return the derivative evaluated at the data points.
}
\item{fout}{
Temporary UNIX file name for the output.
}
\item{run}{
Runs the fitting program.
}
\item{just.setup}{
Sets up the input files but does not run the fitting program.
}
\item{just.read}{
Does no fitting, just reads in the results from a previous fit.
}
\item{fitted.values}{
Computes fitted values and residuals.
}
\item{all.fits}{
Includes all polished fits in the output file not just the best one.
}
\item{greedy}{
A value of zero fits the full model by nonlinear least squares. A
positive value uses the greedy algorithm to fit hidden units in chunks
of size greedy, sequentially adding hidden units fit to the residuals
of the previous fit.
}
\item{seed}{
Seed used in generating the random parameter starts.
}}
\value{
Object of class nnreg. The component model is a list of the parameters
for each fitted model. Columns of the components: residuals and
predicted values, correspond to the different fitted models. Each
component model is of class netfit. The best model number as judged by
minimizing the GCV statistic is also returned.

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
summary of the polished values.
}
\item{fitted.values}{
Predicted values from the fit.
}
\item{residuals}{
Residuals from the fit.
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
Number of different model specifications.
}
\item{lags}{
Time lags used in the x matrix, if a time series model.
}
\item{seed}{
Seed used in generating the random parameter starts.
}
\item{best.model}{
Number of the best model based on GCV with cost=2.
}}
\section{Side Effects}{
This function does the bulk of the computation using a stand-alone
FORTRAN program running in the UNIX shell. This operation is
transparent to the user. For large problems the input files can be
setup using this function and the fitting program can be run
separately in the background.
}
\description{
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


}
\references{
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
plot(fit) # plots fit and residuals


nnreg(as.matrix(BD[,1:4]),BD[,5],2,4) -> fit # fitting DNA strand
# displacement amplification surface to various buffer compositions
plot(fit) # plots fit and residuals
}
\keyword{neural}
% Converted by Sd2Rd version 0.2-a3.
