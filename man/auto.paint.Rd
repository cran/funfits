\name{auto.paint}
\title{
Data frame of paint measurements from the hood of a car.
}
\section{SUMMARY}{
The \code{auto.paint} data frame has 80 rows and 2 columns.
Two paint measurements are the thickness and the distinctness of image  
taken from the hood of a car at the GM research lab.
}
\section{"DATA}{
This data frame contains the following columns:
}
\arguments{
\item{thick}{
Thickness of the paint.
}
\item{DOI}{
Distinctness of image of the paint.
}}
\section{SOURCE}{
James Wendelberger


General Motors research lab
}
\seealso{
nkreg
}
\examples{
nkreg( auto.paint$thick, auto.paint$DOI, 5.0)-> look
}
\keyword{FUNFITS}
% Converted by Sd2Rd version 0.2-a3.
