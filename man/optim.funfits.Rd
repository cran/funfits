\name{optim.funfits}
\title{
Finds optimum within a region
}
\usage{
optim.funfits(fit, start, maximize=T, lower, upper, \dots)
}
\arguments{
\item{fit}{
Object from a fitting procedure. A predict function for the fit 
must be available.
}
\item{start}{
Starting value for the search, the default is the middle of the region.
}
\item{maximize}{
Default is to search for the maximum.
}
\item{lower}{
Lower bound of parameters for optimization, default is the minimum of the data.
}
\item{upper}{
Upper bound of parameters for optimization, default is the maximum of the data.
}
\item{...}{
Optional arguments.
}}
\value{
Returns a list with the followin values:

\item{parameters}{
Final value of parameters at over which the optimization takes place.
}
\item{objective}{
Final value of the objective function.
}
\item{message}{
Statement of the reason for termination.
}
\item{grad.norm}{
Final norm of the objective gradient.
}
\item{iterations}{
Total number of iterations before terminiation.
}
\item{f.evals}{
Total number of residual evaluations before the termination.
}
\item{g.evals}{
Total number of jacobian evaluations before the termination.
}
\item{hessian}{
Final value of the Hessian matrix, only if hessian is supplied intially.
}
\item{scale}{
Final value of scale vector. 
}}
\description{
The object from a fitting procedure must have a predict function.
Optim calls the function nlminb which finds a local minimum of a 
smooth nonlinear function subject to bounded-constrained paramenters.
}
\seealso{
nlminb, predict.tps, predict.krig, predict.nnreg
}
\examples{
tps(BD[,1:4],BD$lnya,scale.type="range") -> fit # fitting a DNA strand
# displacement amplification surface to various buffer compositions
surface(fit) # plots fitted surface and contours
optim(fit) # find surface optimum
}
\keyword{FUNFITS}
% Converted by Sd2Rd version 0.2-a3.
