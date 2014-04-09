\name{DQ}
\alias{DQ}
\title{
  Dynamic quantile test
}
\description{
Backtest a series of Value at Risk forecasts using the Dynamic Quantile Test of Engle and Manganelli
Define a variabele HIT_t = Y_t < VAR_t - prob
Then it should not be possbile to predict HIT_t based on information known at time = t-1. \cr
This may be tested with a regression: \cr
HIT_t = b1 + b2*VAR_t + b3*VAR_t-1 + .... + bn*VAR_t-x + bn+1*HIT_t-1 + ..... + bn+z*Y_t ..... \cr
None of the included variables should be significant in explaining HIT_t
}
\usage{
DQ(Y, VAR, prob, intercept=T, nVAR=0, nHIT=c(1,2,3,4), nY=F)
}
\arguments{
\item{Y}{A series of observed values.}
\item{VAR}{A series of predicted Value at Risk estimates.}
\item{prob}{The probability level for which one do not expect the series to exceed the VaR estimate.}
\item{intercept}{Logical, should an intercept be included in the regression}
\item{nVAR}{Which lagged VaR estimates should be included in the regression. FALSE means that no lagged values are included. 
A single numeric, n, gives all lagged values between 1 and n. An array with numerics, c(x,y,z,...) includes lag x,y,z,..}
\item{nHIT}{Which lagged HIT variables should be included in the regression. FALSE means that no lagged values are included. 
A single numeric, n, gives all lagged values between 1 and n. An array with numerics, c(x,y,z,...) includes lag x,y,z,..}
\item{nY}{Which lagged values of the observed series should be included in the regression. FALSE means that no lagged values are included. 
A single numeric, n, gives all lagged values between 1 and n. An array with numerics, c(x,y,z,...) includes lag x,y,z,..}
}
\references{
Robert F. Engle and Simone Manganelli 
CAViaR: Conditional Autoregressive Value at Risk by Regression Quantiles
Journal of Business & Economic Statistics
Vol. 22, No. 4 (Oct., 2004), pp. 367-381
Published by: American Statistical Association
Article Stable URL: http://www.jstor.org/stable/1392044
}
\value{
Returns a list with the value of the coefficients of the regression, the test observator of Engle and Manganelli (2004) and the probability 
level of observering those coefficients given that their "true" conterparts all are zero.
}
\author{Steinar Veka}
\examples{
}