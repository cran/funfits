\name{predict.netfit}
\title{
Evaluation of a neural net fit based on a model component 
}
\usage{
predict.netfit(fit, x, derivative=0)
}
\arguments{
\item{fit}{
A model component of an nnreg object
}
\item{x}{
A matrix of x values on which to evaluate the neural net surface.
}
\item{derivative}{
Derivative of function are returned if derivative=1.
}}
\value{
A vector of predicted responses. If derivative=1 a vector of derivatives 
or a matrix of partial derivatives is returned.
}
\seealso{
nnreg, predict.nnreg
}
\examples{
nnreg(ozone$x,ozone$y,1,2) -> fit           # nnreg fit
cbind(seq(87,89,,10),seq(40,42,,10)) -> x   # new x matrix
predict(fit$model[[1]],x) -> out            # evaluate model 1 at x
}
\keyword{neural}
% Converted by Sd2Rd version 0.2-a3.
