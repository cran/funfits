\name{set.panel}
\title{
Specify a panel of plots
}
\usage{
set.panel(m=1, n=1, relax=F)
}
\arguments{
\item{m}{
Number of rows in the panel of plots
}
\item{n}{
Number of columns in the panel.
}
\item{relax}{
If true and the par command is already set for multiple plots, 
then the set.panel command is ignored. The default is 
relax set to false. 
}}
\description{
After set.panel is called, the graphics screen is reset to
put  plots  according to a m x n table. Plotting starts in
the upper left hand corner and proceeds row by row.  After
m x n plots have been drawn, the next plot will erase the
window and start in the 1,1 position again.  This function
is just a repackaging for specifying the mfrow argument to
par.


Setting up a panel of plots is a quick way to  change  the
aspect  ratio  of  the graph (ratio of height to width) or
the size. For example, plotting 2 plots to a page produces
a  useful  size  graph  for including in a report. You can
print out the graphs at any stage without having  to  fill
up the entire window with plots.
}
\section{Side Effects}{
The function will echo your choice of m and n to the terminal. 
}
\seealso{
par
}
\examples{
set.panel(5,2) #divide screen to hold 10 plots where there are 5 rows 
	       #and 2 columns
set.panel() #reset screen to one plot per screen
}
\keyword{FUNFITS}
% Converted by Sd2Rd version 0.2-a3.
