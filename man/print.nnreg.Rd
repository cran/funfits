\name{print.nnreg}
\title{
Print nnreg fit results
}
\usage{
print.nnreg(fit)
}
\arguments{
\item{fit}{
Object from nnreg fitting.
}}
\value{
Results from nnreg. This includes input parameters and results of each of
the polished fits.
}
\seealso{
print, summary.nnreg, nnreg
}
\examples{
nnreg(ozone$x,ozone$y,1,3) -> fit # nnreg fit
print(fit) # print the summary
fit # this will too
}
\keyword{neural}
% Converted by Sd2Rd version 0.2-a3.
