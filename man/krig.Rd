\name{krig}
\title{
Kriging surface estimate 
}
\usage{
krig(x, Y, cov.function, lambda=NA, cost=1, knots, 
weights=rep(1, length(Y)), m=2, return.matrices=T, 
nstep.cv=50, scale.type="user", x.center=rep(0, ncol(x)), 
x.scale=rep(1, ncol(x)), rho=NA, sigma2=NA, \dots)
}
\arguments{
\item{x}{
Matrix of independent variables.
}
\item{Y}{
Vector of dependent variables.
}
\item{cov.function}{
Covariance function for data in the form of an S-PLUS function (see exp.cov). 
}
\item{lambda}{
Smoothing parameter that is the ratio of the error variance (sigma**2)
to the scale parameter of the 
covariance function. If omitted this is estimated by GCV.
}
\item{cost}{
Cost value used in GCV criterion. Corresponds to a penalty for 
increased number of parameters.
}
\item{knots}{
Subset of data used in the fit.
}
\item{weights}{
Weights are proportional to the reciprocal variance of the measurement 
error. The default is no weighting i.e. vector of unit weights.
}
\item{m}{
A polynomial function of degree (m-1) will be 
included in the model as the drift (or spatial trend) component.
}
\item{return.matrices}{
Matrices from the decompositions are returned. The default is T. 
}
\item{nstep.cv}{
Number of grid points for minimum GCV search.
}
\item{scale.type}{
The independent variables and knots are scaled to the specified scale.type.
By default the scale type is "unit.sd", whereby the data is scaled to
have mean 0 and standard deviation 1. Scale type of "range" scales the
data to the interval (0,1) by forming (x-min(x))/range(x) for each x.
Scale type of "user" allows specification of an x.center and x.scale by the
user. The default for "user" is mean 0 and standard deviation 1. Scale
type of "unscaled" does not scale the data.
}
\item{x.center}{
Centering values are subtracted from each column of the x matrix.
}
\item{x.scale}{
Scale values that divided into each column after centering.
}
\item{rho}{
Scale factor for covariance.
}
\item{sigma2}{
Variance of e.
}
\item{...}{
Optional arguments. Theta can be specified. If the cov.parameters are 
specified this list is assumed to be arguments to the covariance 
function.
}}
\value{
A list of class krig. This includes the predicted surface of
fitted.values and the residuals. The results of the grid
search to minimize the generalized cross validation function is
returned in gcv.grid.

\item{call}{
Call to the function
}
\item{y}{
Vector of dependent variables.
}
\item{x}{
Matrix of independent variables.
}
\item{weights}{
Vector of weights.
}
\item{knots}{
Locations used to define the basis functions. 
}
\item{transform}{
List of components used in centering and scaling data.
}
\item{np}{
Total number of parameters in the model.
}
\item{nt}{
Number of parameters in the null space.
}
\item{matrices}{
List of matrices from the decompositions (D, G, u, X, qr.T).
}
\item{gcv.grid}{
Matrix of values used in the GCV grid search. The first column
is the grid of lambda values used in the search, the second column 
is the trace of the A matrix, the third column is the GCV values and
the fourth column is the estimated variance.
}
\item{cost}{
Cost value used in GCV criterion.
}
\item{m}{
Order of the polynomial space: highest degree polynomial is (m-1).
}
\item{eff.df}{
Effective degrees of freedom of the model.
}
\item{fitted.values}{
Predicted values from the fit.
}
\item{residuals}{
Residuals from the fit.
}
\item{lambda}{
Value of the smoothing parameter used in the fit.
}
\item{yname}{
Name of the response.
}
\item{cov.function}{
Covariance function of the model.
}
\item{beta}{
Estimated coefficients in the ridge regression format
}
\item{d}{
Esimated coefficients for the polynomial basis functions that span the
null space
}
\item{fitted.values.null}{
Fitted values for just the polynomial part of the estimate
}
\item{trace}{
Effective number of parameters in model.
}
\item{c}{
Estimated coefficients for the basis functions derived from the
covariance.
}
\item{coefficients}{
Same as the beta vector.
}
\item{just.solve}{
Logical describing if the data has been interpolated using the basis 
functions. 
}
\item{shat}{
Estimated standard deviation of the measurement error (nugget effect).
}
\item{sigma2}{
Estimated variance of the measurement error (shat**2).
}
\item{rho}{
Scale factor for covariance.  COV(h(x),h(x')) = rho*cov.function(x,x')
}
\item{mean.var}{
Normalization of the covariance function used to find rho.
}
\item{best.model}{
Vector containing the value of lambda, the estimated variance of the 
measurement error and the scale factor for covariance used in the fit.
}}
\description{
The kriging model is Y(x)= P(x) + Z(x) + e where Y is the dependent
variable observed at location x, P is a low order polynomial, Z is a
mean zero, Gaussian field with covariance function
K and e is assumed to be independent normal
errors. The estimated surface is the best linear unbiased estimate (BLUE) 
of P(x) + Z(x) given the
observed data. For this estimate K, is taken to be rho*cov.function
and the errors have variance sigma^2. If these parameters are omitted
in the call, then they are estimated in the following way. If lambda is
given, then sigma2 is estimated from the residual sum of squares
divided by the degrees of freedom associated with the residuals. 
Rho is found as the difference between the sums of squares of the 
predicted values having subtracted off the polynomial part and sigma2. 


WARNING: The covariance functions often have a nonlinear parameter that 
controls the strength of the correlations as a function of separation, 
usually refered to as the range parameter. This parameter must be 
specified in the call to krig and will not be estimated. 


}
\references{
See "Additive Models" by Hastie and Tibshirani, "Spatial Statistics" by   
Cressie and the FUNFITS manual.
}
\seealso{
summary.krig, predict.krig, predict.se.krig, plot.krig, surface.krig
}
\examples{
#2-d example
krig(ozone$x, ozone$y, exp.cov) -> fit # fitting a surface to ozone 
# measurements.
plot(fit) # plots fit and residuals
# data using a Gaussian covariance
# first make up covariance function
test.cov <- function(x1,x2){exp(-(rdist(x1,x2)/.5)**2)}
krig(flame$x, flame$y, test.cov) -> fit.flame
surface(fit.flame)
}
\keyword{nonparametric}
% Converted by Sd2Rd version 0.2-a3.










