\name{make.surface.grid}
\title{
Makes a grid surface that can be plotted.
}
\usage{
make.surface.grid(grid.list, X, nx=30, ny=30, info.list=F, FUN=median)
}
\arguments{
\item{grid.list}{
A list with as many components as variables describing the surface.
All components should have a single value except the two that give the
grid points for evaluation. If the matrix or data frame has column names, 
these must appear in the grid list.
}
\item{X}{
Matrix of independent variables.
}
\item{nx}{
Number of grid points evenly spaced over the range of x-values.
}
\item{ny}{
Number of grid points evenly spaced over the range of y-values.
}
\item{info.list}{
Surface information, component of attributes.
}
\item{FUN}{
}}
\value{

}
\section{Side Effects}{
}
\description{
}
\seealso{
grid.list, predict.surface, predict.surface.se
}
\examples{
}
\keyword{FUNFITS}
% Converted by Sd2Rd version 0.2-a3.
