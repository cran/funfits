\name{addreg}
\title{
Additive nonparametric regression
}
\usage{
addreg(x, y, lam, nback=20, tol=1e-05, start, cost=1)
}
\arguments{
\item{x}{
Matrix of independent variables
}
\item{y}{
Vector of independent variables
}
\item{lam}{
Vector of length P with values of the smoothing parameters. If missing
then the parameters are found by generalized cross-validation.
}
\item{nback}{
Maximum number of iterations for backfitting. 
}
\item{tol}{
Tolerance value for judging convergence of backfitting.
}
\item{start}{
Matrix of starting values for the function estimates
}
\item{cost}{
Cost for each degree of freedom in the GCV criterion.
}}
\value{
A list of class addreg. Components fitted.values and residuals
contain the results of the fit. predicted.comp is matrix whose columns
contain the estimates for the individual functions. The component
trace has the effective number of parameters associated with each
smoother and
lambda has the corresponding smoothing parameters. Convergence
information is in the component converge. 


}
\description{
Estimates the univariate functions using backfitting. If lambda is
omitted in the call then at each step the smoothing parameter is
estimated by minimizing the GCV criterion:
          GCV(lambda)=(1/n)* RSS/(n- p*2 + cost*(sum( J_i))**2


Here J_i = tr(A_i) -2 where  A_i is the smoothing matrix for variable
i.


If p=1 and cost =1 this is the usual GCV function. A cost > 1 has the
effect of being more conservative in choosing estimates that deviate
from a linear function.  Note that the effective number of parameters
in the total model can be approximated by sum( tr( A_i)).  
}
\seealso{
plot.addreg, gam, predict.addreg
}
\examples{
# Additive model fit to the minitriathon data. A prediction of 
# run times based on the swim and bike performance.
# The smoothing parameters found by GCV with cost=2. 
cbind( minitri$swim, minitri$bike)-> x
addreg( x,minitri$run,cost=2)-> out
plot( out) # summary plots of fits
list( c(10,25), c(35,55))-> grid.list
predict.surface(out, grid.list)-> out.p 
persp( out.p) # perspective plot of fitted surface
}
\keyword{additive}
% Converted by Sd2Rd version 0.2-a3.










