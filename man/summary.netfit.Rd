\name{summary.netfit}
\title{
Summary of parameters of netfit object
}
\usage{
summary.netfit(fit, standardized=F)
}
\arguments{
\item{fit}{
A model component of an nnreg object
}
\item{standardized}{
Data is standardized to have mean zero and standard deviation one.
}}
\value{
A list of netfit paramters.
}
\seealso{
nnreg, summary.nnreg, predict.netfit
}
\examples{
nnreg(ozone$x,ozone$y,1,3) -> fit           # nnreg fit
summary(fit$model[[1]],fit$x)               # summary of model 1 
# parameters
}
\keyword{neural}
% Converted by Sd2Rd version 0.2-a3.
