\name{write.data}
\title{
Writes data set out to a file
}
\usage{
write.data(x, file)
}
\arguments{
\item{x}{
The data set to be printed.
}
\item{file}{
Name of UNIX output file in quotes. If this omitted, the data is written to
a file with the name of the data set with a .output added to the
end.
}}
\section{Side Effects}{
A UNIX file will be created  with  a printed version of the data set.
}
\description{
This is a friendly version for creating output files using the sink
and print functions.
To avoid funny printing due to "page breaks," the page length is
reset within this function. Besides this detail, the output to
the file will be the same as what would appear on the screen
if the data set name was typed. Thus, if
a matrix is written out, there will be row and column numbers
included.  Lists will have the same component formatting that
appears when the list is printed to the screen.
}
\seealso{
read.data, write, print, sink
}
\examples{


Write the west coast climate data out to a file called temp.out:


write.data(climate.w,"temp.out")


Write out the triathlon data. The default UNIX file name will be
minitri.output:


write.data(minitri)


}
\keyword{IO,}
% Converted by Sd2Rd version 0.2-a3.
