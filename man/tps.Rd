\name{tps}
\title{
Thin plate spline regression
}
\usage{
tps(x, y, lambda=NA, df=NA, cost=1, knots, weights=rep(1, length(y)), m, 
power,scale.type="unit.sd", x.center, x.scale, return.matrices=T, 
nstep.cv=80, method="GCV", rmse=NA, link.matrix=NA, verbose=F, 
subset=NULL, tol=0.0001, print.warning=T)
}
\arguments{
\item{x}{
variables. 
}
\item{Y}{
Vector of dependent variables.
}
\item{lambda}{
Smoothing parameter. If omitted this is estimated by GCV. Lambda=0 gives
an interpolating model.
}
\item{df}{
Specifies the effective degrees of freedom associated with the spline 
estimate. This parameter is an alternative to specifying lambda directly. 
}
\item{cost}{
increased number of parameters. 
}
\item{knots}{
Subset of data used in the fit.
}
\item{weights}{
Vector - default is no weighting i.e. vector of unit weights
(Weights are in units of reciprocal variance.)
}
\item{m}{
Order of spline surface, default is 2 corresponding to a linear (m-1=1) 
base polynomial model. If power is specified (m-1) will be the degree
of the polynomial null space. 
}
\item{power}{
Power used for the norm in the radial basis functions 
the default is 2*m-d and this will result in true thin-plate splines.  
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
Value subtracted from each column of the x matrix.
}
\item{x.scale}{
Value divided into each column for scaling.
}
\item{return.matrices}{
Matrices from the decompositions are returned.
}
\item{nstep.cv}{
Number of grid points for initial GCV grid search.
}
\item{method}{
The method for determining the smoothing parameter to evaluate the 
spline estimate. Choices are "GCV" (the default), "RMSE" and "pure error". 
The method may also be  set implicitly by the values of other 
arguments 
}
\item{rmse}{
~Describe rmse here
}
\item{link.matrix}{
A matrix that relates the function evaluated at the x values to the 
mean of y. This option is used when the mean of the observed data is a 
linear combination of the function evaluated at the x values.  
}
\item{verbose}{
If true prints all kinds of intermediate calculations. This is mainly for 
trouble shooting
}
\item{subset}{
A logical vector indicating the subset of data to use for fitting. 
}
\item{tol}{
Tolerance for convergence of the golden section and bisection searches 
in the GCV function and the df.to.lambda function. 


}
\item{print.warning}{


}}
\value{
A list of class tps. This includes the predicted surface of
fitted.values and the residuals. The results of the grid
search minimizing the generalized cross validation function is
returned in gcv.grid.

\item{call}{
Call to the function
}
\item{x}{
Matrix of independent variables.
}
\item{y}{
Vector of dependent variables.
}
\item{form}{
Logical denoting that the form of the model. Default is a thin plate spline.
}
\item{cost}{
Cost value used in GCV criterion. 
}
\item{m}{
Order of spline surface.
}
\item{trace}{
Effective number of parameters in model.
}
\item{trA2}{
trace of the square of the smoothing matrix, tr(A(lambda)**2)
}
\item{yname}{
Name of the response.
}
\item{weights}{
Vector of weights.
}
\item{knots}{
Subset of data used.
}
\item{transform}{
List of components used in scaling data.
}
\item{power}{
2*m-d unless specified explicitly in the call.
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
\item{beta}{
All coefficients i.e. corresponding to polynomial and radial basis functions.
}
\item{d}{
Parameters of the polynomial null space model. The powers for these
terms are accessible from the attribute matrix returned by make.tmatrix.
}
\item{c}{
Parameters for the radial basis terms. 
}
\item{coefficients}{
Same as beta. 
}
\item{just.solve}{
Logical indicating lambda=0 i.e. an interpolating function is returned.
}
\item{lambda.est}{
A matrix giving the values  of lambda found for different methods and the
corresponding estimates of sigma and the GCV function.
}
\item{shat}{
Estimated standard deviation of the errors using the specified value of 
lambda
}
\item{shat.pure.error}{
Estimate of standard deviation of the errors using replicated points.
the value is NA if there are no replicates
}
\item{GCV}{
Value of Generalized Cross Validation at lambda
}
\item{rmse}{
Root mean squared error used as a target for choosing lambda.
}
\item{q2}{
}
\item{press}{
}}
\value{
best.model 
Value of lambda used in the fit.
}
\description{
A thin plate spline is result of minimizing the residual sum of
squares subject to a constraint that the function have a certain
level of smoothness (or roughness penalty). Roughness is
quantified by the integral of squared m^th order derivatives. For one
dimension and m=2 the roughness penalty is the integrated square of
the second derivative of the function. For two dimensions the
roughness penalty is the integral of 
        (Dxx(f))^2 + 2(Dxy(f))^2 + (Dyy(f))^2
(where Duv denotes the second partial derivative with respect to u
and v.) Besides controlling the order of the derivatives, the value of
m also determines the base polynomial that is fit to the data.
The degree of this polynomial is (m-1).


The smoothing parameter controls the amount that the data is
smoothed. In the usual form this is denoted by lambda, the Lagrange
multiplier of the minimization problem. Although this is an awkward
scale, lambda =0 corresponds to no smoothness constraints and the data
is interpolated.  lambda=infinity corresponds to just fitting the
polynomial base model by ordinary least squares. 
}
\references{
See "Nonparametric Regression and Generalized Linear Models" 
by Green and Silverman.


See "Additive Models" by Hastie and Tibshirani.
}
\seealso{
summary.tps, predict.tps, predict.se.tps, plot.tps, surface.tps
}
\examples{


#1-d example


tps( rat.diet$t, rat.diet$trt) # lambda found by GCV
tps( rat.diet$t, rat.diet$trt, df=6) # lambda chosen so that spline has 6 
                                     # degrees of freedom


#2-d example
tps(ozone$x, ozone$y) -> fit # fits a surface to ozone measurements.
plot(fit) # plots fit and residuals.


#4-d example
tps(BD[,1:4],BD$lnya,scale.type="range") -> fit # fits a surface to
# DNA strand displacement amplification as a function of various 
# buffer components.
surface(fit)  
# plots fitted surface and contours




#2-d example using a reduced set of basis functions
r1 <- range(flame$x[,1])
r2 <-range( flame$x[,2])
g.list <- list(seq(r1[1], r1[2],6), seq(r2[1], r2[2], 6))
make.surface.grid(g.list) -> knots  # these knots are a 6X6 grid over
# the ranges of the two flame variables
tps(flame$x, flame$y, knots=knots, m=4) -> out 


# here is an example using a link matrix
# 
x<- seq( 0,1,.02)
f<- 8*(x**2)*(1-x)
M <- matrix( .02, ncol=50, nrow=50)
M[ col(M)< row( M)]<- 0 
set.seed( 123)
y<- M\%*\%f + error*.1 # So y is approximately the integral of f pus error
ex.out<-tps( x,y,link.matrix=M)
#Note: predict will give the spline NOT the predicted values for y!








}
\keyword{splines,}
% Converted by Sd2Rd version 0.2-a3.




