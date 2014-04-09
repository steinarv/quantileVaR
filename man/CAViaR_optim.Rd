\name{CAViaR_optim}
\alias{CAViaR_optim}
\title{
  Optimize a model
}
\description{
Finds optimal parameters for one of the following quantile regression models: \cr
	Symmetric Absolute Value: f_t(b) = b1 + b2 f_t-1(b) + b3 |y_t-1| \cr
	Asymmetric Slope: f_t(b) = b1 + b2 f_t-1(b) + b3 (y_t-1)[+] + b4 (y_t-1)[-] \cr
	Indirect GARCH(1,1): f_t(b) = (b1 + b2 f_t-1(b)^2 + b3 y_t-1^2)^0.5 \cr
	Adaptive: f_t(b1) = f_t-1(b1) + b1 {[1-exp(G[y_t-1 - f_t-1(b1)])]^-1 - prob} \cr
	
	f(b) is the quantile for probability level (prob) given b1,b2,... \cr
	(x)[+] = max(x,0) \cr
	(x)[-] = min(x,0) \cr
}
\usage{
CAViaR_optim <- function(y, Model=1, prob=0.05)
}
\arguments{
\item{y}{
A numeric array with values intended for estimating the model.}
\item{Model}{
Scalar between 1 and 5. 1 = Symmetric Absolute Value, 2 = Asymmetric Slope \cr
3 = Indirect GARCH(1,1), 4 = Adaptive, 5 = Symmetric Absolute Value with b2 set to zero.}
\item{prob}{
Scalar between 0 and 1 gives the quantile.}
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
Returns a list with 4 elements: \cr
$bestRQ: The score used to optimize the model parameter. \cr
$bestVals: The ten best parameter vectors and their corresponding RQ score. \cr
$bestPar: The optimal parameter vector. \cr
$VAR: In sample Value at Risk (or quantile values) given the optimal parameter vector. \cr

}
\author{Steinar Veka}
\examples{
dax <- c(EuStockMarkets[,2]) \cr
r_dax <- diff(log(dax)) \cr
fit = CAViaR_optim(r_dax, Model=1, prob=0.05) \cr
mean(r_dax<fit$VAR) \cr
plot(r_dax, type="l", ylim=range(r_dax)*1.2) \cr
lines(fit$VAR, col=3) \cr
}