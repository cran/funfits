\name{plot.window}
\title{
Creates a motif graphics window.
}
\usage{
plot.window(screen.name)


}
\arguments{
\item{screen.name}{
The name of the screen enclosed in quotes that will have the plot
window.


}}
\description{
This function will take screen names either from the UNIX environment
variable DISPLAY or from the name supplied when plot.window is called.
This name will override the name in the DISPLAY variable.  


After successfully opening the window, some instructions about
printing are given. The graphics window is set to be the same
ratio as an 8 1/2 by 11 inch sheet of paper, and the plots are in
portrait format. The window defaults are in the data
set .motif.options.


}
\section{NOTE}{
If you are running Splus on a machine different from the one you are
sitting at, the remote computer may need permission to write to the screen.
On your local machine you can set the priviledges by typing xhost +
in UNIX before running S. This will allow any other machine on your
network to write to your screen.
}
\seealso{
set.panel


}
\examples{
Open a plot window on your machine:


plot.window()


You are running S on a remote machine but would like a graphics window
on your local machine. The machine you are sitting at is called
RA_FISHER.stat. In UNIX on your local machine type in  xhost +  .
Now run Splus
in a terminal window that is connected to the remote machine. 


plot.window("RA_FISHER.stat:0")


Note the addition of the :0 to the machine name.


If this does not work, the name may not be in the right form. The
sure fire way is to use the IP number of your machine in place of
the alphabetic name. Typically an IP number is a set of 4 1-3 digit
numbers separated by periods. E.g 151.1.62.60 add :0 to this to
specify the screen:


plot.window("151.1.62.60:0")


}
\keyword{graphics,}
% Converted by Sd2Rd version 0.2-a3.
