\name{BD}
\title{
Data frame of the effect of buffer compositions on DNA strand displacement
amplification.
}
\section{SUMMARY}{
The \code{BD} data frame has 89 rows and 5 columns. There are 89 runs with 
four buffer components (KCL, MgCl2, KP04, dnTP) systematically varied in
a space-filliing design. The response is the DNA amplification rate. 
}
\section{"DATA}{
This data frame contains the following columns:
}
\arguments{
\item{KCl}{
Buffer component.
}
\item{MgCl2}{
Buffer component.
}
\item{KPO4}{
Buffer component.
}
\item{dNTP}{
Buffer component, deoxyribonucleotides.
}
\item{lnya}{
Exponential amplification rate on a log scale, i.e. the actual amplification
rate.
}}
\section{SOURCE}{
Becton Dickinson Research Center
Research Triangle Park, NC
}
\seealso{
tps, krig, nnreg, optim
}
\examples{
tps(BD[,1:4],BD$lnya,scale.type="range") -> fit # fitting a DNA strand
# displacement amplification  surface to various buffer compositions
surface(fit)  # plots fitted surface and contours
}
\keyword{FUNFITS}
% Converted by Sd2Rd version 0.2-a3.
