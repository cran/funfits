\name{rossler.state}
\title{
A multivariate  time series of the form W_t=f(W_t-1) + e_t.


The Rossler system is a system of three first order nonlinear ODEs.
The system has been numerically integrated with a time step size of .1 
and noise that is distributed normal with mean =0 and sigma=.1 is added
at every time step. The map f is highly nonlinear and has complicated 
structure not easily fit by simple regression methods. 
The time series is 
created by 
sampling every 5 time steps. The rossler object is a matrix with columns
the x,y and z values of the system at the sampling times. 
}
\keyword{FUNFITS}
% Converted by Sd2Rd version 0.2-a3.
