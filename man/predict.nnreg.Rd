\name{predict.nnreg}
\title{
Evaluation of neural net surface.
}
\usage{
predict.nnreg(out, x, model=NA, derivative=0, type="full")
}
\arguments{
\item{out}{
Fitted nnreg object.
}
\item{x}{
Matrix of x values on which to evaluate the neural net surface.
}
\item{model}{
Model number to use in predicting. Default is the best model based on GCV(2).
}
\item{derivative}{
Derivative of function is returned if derivative=1. 
}
\item{type}{
Form of predictions. Default is the prediction for the independent
variable. If type="terms" the individual values for the hidden units
are calculated.


}}
\value{
Vector of predicted responses. If derivative=1 a vector of derivatives
or a matrix of partial derivatives is returned. If type="terms" a
list with components:  u a matrix with the projections of the
independent vectors plus the offset ( X\%*\% gamma + gamma_0) for each
hidden unit,   
yhat, a matrix where the columns vectors are the evaluation of each hidden
unit and constant, the value of the constant (intercept) in the model. 
}
\seealso{
nnreg, predict.surface
}
\examples{
nnreg(ozone$x,ozone$y,1,2) -> fit           # nnreg fit
cbind(seq(87,89,,10),seq(40,42,,10)) -> x   # new x matrix
predict(fit,x) -> out                       # evaluate fit at x
}
\keyword{neural}
% Converted by Sd2Rd version 0.2-a3.
