\name{predict.se}
\title{
Standard errors of predictions
}
\usage{
predict.se(object, \dots)
}
\arguments{
\item{object}{
A fitted model object of a certain class
}}
\value{
A vector of standard errors for the predicted values.
}
\description{
This function is generic and will call the appropriate function to calculate
the standard errors for the object class.
}
\seealso{
predict, predict.se.tps, predict.se.krig
}
\examples{
tps(ozone$x,ozone$y) -> fit                 # tps fit
predict.se(fit)                             # std errors of predictions
cbind(seq(87,89,,10),seq(40,42,,10)) -> x   # new x matrix
predict.se(fit,x) -> out                    # std errors of predictions


tps(BD[,1:4],BD$lnya,scale.type="range") -> fit # fitting a DNA strand
# displacement amplification surface to various buffer compositions
predict.se(fit) -> out # std erros of predictions
}
\keyword{FUNFITS}
% Converted by Sd2Rd version 0.2-a3.
