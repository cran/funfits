\name{yline}
\title{
Draw a horizontal line
}
\usage{
yline(y, \dots)
}
\arguments{
\item{y}{
Value on y axis specifying location of horizontal line.
}
\item{...}{
Any ploting options for abline function.
}}
\examples{
plot(ozone$x[,2],ozone$y) #plot ozone longitude vs ozone value
yline(40) #plots horizontal line at ozone value of 40
}
\keyword{graphics,}
% Converted by Sd2Rd version 0.2-a3.
