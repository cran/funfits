\name{bplot}
\title{
boxplot
}
\usage{
bplot(x, \dots, xpos=NA, width, label, by, srt=0, add=F, space=0.25, sort.names=T, xlab="", ylab="")
}
\arguments{
\item{x}{
Vector, matrix,  list or data frame
}
\item{...}{
Additional data sets to be plotted.
}
\item{xpos}{
The boxes will be plotted vertically and xpos gives the x locations for 
their centers. If omitted the boxes are equally spaced at integer 
values. 
}
\item{width}{
Width of boxplots if omitted then the width is a reasonable fraction of the 
distance between boxes and is set by the space argument.
}
\item{label}{
Labels under each boxplot. If missing the columns names or components of 
x are used. 
}
\item{by}{
If x is a vector a vector of categories that will be used to split the 
x values into separate groups. Of course box plots will be made of the 
groups. 
}
\item{srt}{
Rotate the labels. Default is to put them horizontal. 
}
\item{add}{
Do not create a new plots just add the boxplots to a current plot. Note 
that the xpos argument may be useful in this case. 
}
\item{space}{
Space between boxplots. 
}
\item{sort.names}{
Plot the boxplots sorted in alphabetic order by their labels. 
}
\item{xlab}{
Label for the x-axis
}
\item{ylab}{
Label for the y-axis
}}
\description{
This function was created as a complement to the usual S function for 
boxplots. The current function makes it possible to put the boxplots at 
unequal x-positions. This is useful for visually grouping a large set of 
boxplots into several groups. Also placement of the boxplots with 
respect to the x-axis can add information to the plot.  Another aspect 
is the emphasis on data structures for groups of data. One useful 
feature is the by option to break up the x vector into distinct  groups.  
If 5 or less observations are in a group the points themselves are 
plotted instead of a box.  


Finally it should be noted that this function constructs the boxplots 
based on percentiles of the data rather than using the hinge and fence 
arcana from EDA. This is probably heretical to some but it makes it 
much easier to  teach undergraduate students what these are!  


The function is simple enough so that one _could_ hingeify the boxes. 
Also it would be simple to modify this function to get boxes 
horizontally instead of vertically.  


 
}
\seealso{
boxplot,lplot 
}
\examples{
# boxplots of the times for different events in a minitriathlon
bplot( minitri)
# some boxplots of uniform random sample with different power 
# transformations
temp<- runif(50)
bplot( temp, temp**1.5, temp**2, xpos=c(1,1.5,2.0), labels=c("uniform","1.5 
power", "square") )  
}
\keyword{FUNFITS}
% Converted by Sd2Rd version 0.2-a3.






