\name{xline}
\title{
Draw a vertical line
}
\usage{
xline(x, \dots)
}
\arguments{
\item{x}{
Value on x axis specifying location of vertical line.
}
\item{...}{
Any ploting options for abline function.
}}
\examples{
plot(ozone$x[,2],ozone$y) # plot ozone longitude vs ozone value
xline(42) # plots vertical line at longitude value of 42
}
\keyword{graphics,}
% Converted by Sd2Rd version 0.2-a3.
