\name{make.lags}
\title{
Lags vectors and covariates correctly so that an autoregressive model 
can be estimated by regression. 
}
\usage{
make.lags(x, lags, cov,nobs=3500)
}
\arguments{
\item{x}{
Vector or matrix representing a univariate or multivariate time series.
(rows are assumed to idex time)
}
\item{lags}{
Vector of time delays used in reconstruction.
}
\item{nobs}{
Maximum length of time series.
}
\item{cov}{
A vector or matrix of covariates that will be matched with the times for 
the independent varaible
}}
\value{


\item{x}{
Matrix of lagged values of the time series, independent variables.
The covaraites are the last columns of this matrix
}
\item{y}{
Vector of time series values, dependent variables.
}
\item{nvar}{
Number of variables or dimension of x matrix.
}
\item{lags}{
Time delays used in constructing the x matrix.
}
\item{start}{
Observation number of univariate time series used for the start of the 
y vector.
}
\item{end}{
Observation number of univariate time series used for the end of the 
y vector. 
}
\item{skip}{
Information about which columns of the returned X matrix are covariates.


}}
\description{
This function is used to create the appropriate data structure for
a nonlinear autoregressive process of the form X_t = F(X_t-1) + e_t.


}
\seealso{
nnreg, rossler
}
\examples{


make.lags(rossler.state[,1],c(1,2,3)) -> data  
# create  
# 3-d time delay vector model of the x variable of rossler system.
nnreg(data$x,data$y,5,5) -> fit # fit time series model using nnreg.


# fitting a state space model to the rossler state vector
# only one lag is neede in this case. 
make.lags(rossler.state, lags=c(1))-> data
nnreg( data$x, data$y[,1], 5,5)-> fit1
nnreg( data$x, data$y[,2], 5,5)-> fit2
nnreg( data$x, data$y[,3], 5,5)-> fit3


}
\keyword{FUNFITS}
% Converted by Sd2Rd version 0.2-a3.






