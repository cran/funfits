\name{plot.nnreg}
\title{
Diagnoistic and summary plots of a neural net regression object
}
\usage{
plot.nnreg(out, model=out$best.model, main=NA, digits=4, graphics.reset=T, \dots)
}
\arguments{
\item{out}{
A nnreg object
}
\item{model}{
Model number to plot. Default is the best model based on GCV(2).
}
\item{main}{
Title of the plot. Default is the function call.
}
\item{digits}{
Number of significant digits for the RMSE label.
}
\item{graphics.reset}{
Reset to original graphics parameters after function plotting.
}
\item{...}{
Any plotting arguments.
}}
\description{
This function creates four summary plots of the tps object. The first
plot is a plot of the nnreg fit to the data. The second is predicted
values vs residuals. The third is number of parameters vs GCV with cost=1 
and cost=2. The fourth plot  is number of parameters vs root mean squared 
error.
}
\seealso{
nnreg, summary.nnreg
}
\examples{
nnreg(ozone$x,ozone$y,1,2) -> fit # fitting a surface to ozone 
# measurements, from 1 to 2 hidden units
plot(fit) # plots fit and residuals


nnreg(as.matrix(BD[,1:4]),BD[,5],1,5) -> fit # fitting DNA strand
# displacement amplification surface to various buffer compositions
plot(fit) # plots fit and residuals
}
\keyword{neural}
% Converted by Sd2Rd version 0.2-a3.














