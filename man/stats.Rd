\name{stats}
\title{
Calculate summary statistics
}
\usage{
stats(x, by)
}
\arguments{
\item{x}{
The data structure to compute the statistics. This can either be a
vector, matrix (data sets are the columns), or a list (data sets are
the components).
}
\item{by}{
If x is a vector, an optional vector
specifying the categories to divide x into separate data sets.


}}
\value{
A matrix where rows index the summary statistics and the columns
index the separate data sets.


}
\description{
Stats breaks x up into separate data sets and then calls describe
to calculate the statistics. The default set of statistics are the number of
(nonmissing) observations, mean, standard deviation, minimum, lower quartile,
median, upper quartile, maximum, and number of missing observations.  If any
data set is nonnumeric, missing values are returned for the statistics.
The by argument is a useful way to calculate statistics on parts of a
data set according to different cases.


}
\seealso{
describe


}
\examples{
Statistics for a normal random sample of 200 observations:


rnorm(200) -> zork
stats(zork) -> look
look


Summary statistics for the components of the climate data:


stats(climate) -> look
look


Just list the means of the climate variables (second row of look):


look[2,]


Statistics for the actuator thrust by pressure level (there are 10 different
pressure levels).


stats(actuator.calibrate$force,by=actuator.calibrate$pressure)


}
\keyword{statistics,}
% Converted by Sd2Rd version 0.2-a3.
