\name{lplot}
\title{
Labeled plot
}
\usage{
lplot(x, y, labels, tcex=0.7, \dots)


}
\arguments{
\item{x}{
vector of x values
}
\item{y}{
vector of y values
}
\item{labels}{
The labels to put at each of the points in the scatterplot.
If only one value or symbol is given, this value
is used at all the points.
}
\item{tcex}{
The size of the labels relative to the text in the rest of the plot.
}
\item{...}{
Any other arguments to pass to the plot function.


}}
\description{
This function is a convenient combination of the plot and text
functions. If either x or y is categorical then the points are
equally spaced along the relevant axis and categories are labeled
along the axis. This feature is useful for making interaction
plots from an ANOVA or making dot plots of several data sets
side by side.


}
\seealso{
plot, text,
}
\examples{
Plot the latitude verses the minimum average January temperature for
50 US cites and label the plot with the city names.  Also include some
axes labels and make the city labels smaller than the default size so
that they are more readable:


lplot(climate$lat,climate$jan,climate$city,xlab="Latitude",
ylab="Jan. Temp.",tcex=.4)


Plot force measurements from the actuator experiment by the supply
line pressure. Label the points by the type of actuator:


lplot(actuator.exp$press,actuator.exp$force,actuator.exp$act,
xlab="pressure",ylab="force")






}
\keyword{graphics}
% Converted by Sd2Rd version 0.2-a3.
