% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/inequalityMeasureBoot.R
\name{ineq.weighted.boot}
\alias{ineq.weighted.boot}
\title{Weighted inequality measures with bootstrap}
\usage{
ineq.weighted.boot(
  X,
  W = rep(1, length(X)),
  B = 100,
  AF.norm = TRUE,
  Atkinson.e = 1,
  Jenkins.alfa = 0.8,
  Entropy.power = 0.5,
  zeroes = "include",
  Kolm.p = 1,
  Kolm.scale = "Standardization",
  Leti.norm = T,
  AN_Y.a = 1,
  AN_Y.b = 1,
  Apouey.a = 2/(1 - length(W[!is.na(W) & !is.na(X)])),
  Apouey.b = length(W[!is.na(W) & !is.na(X)])/(length(W[!is.na(W) & !is.na(X)]) - 1),
  BL.withsqrt = FALSE,
  keepSamples = FALSE,
  keepMeasures = FALSE,
  conf.alpha = 0.05,
  calib.boot = FALSE,
  Xs = rep(1, length(X)),
  total = sum(W),
  calib.method = "truncated",
  bounds = c(low = 0, upp = 10)
)
}
\arguments{
\item{X}{is a data vector}

\item{W}{is a vector of weights}

\item{B}{is a number of bootstrap samples.}

\item{AF.norm}{(logical). If TRUE (default) then index is divided by its maximum possible value}

\item{Atkinson.e}{is a parameter for Atkinson coefficient}

\item{Jenkins.alfa}{is a parameter for Jenkins coefficient}

\item{Entropy.power}{is a generalized entropy index parameter}

\item{zeroes}{defines what to do with zeroes in the data vector. Possible options are "remove" and "include". See Entropy function for details.}

\item{Kolm.p}{is a parameter for Kolm index}

\item{Kolm.scale}{method of data standardization before computing}

\item{Leti.norm}{(logical). If TRUE (default) then Leti index is divided by a maximum possible value}

\item{AN_Y.a}{is a positive parameter for Abul Naga and Yalcin inequality measure}

\item{AN_Y.b}{is a parameter for Abul Naga and Yalcin inequality measure}

\item{Apouey.a}{is a parameter for Apouey inequality measure}

\item{Apouey.b}{is a parameter for Apouey inequality measure}

\item{BL.withsqrt}{if TRUE function returns index given by BL2, elsewhere by BL (default). See more in details of BL function.}

\item{keepSamples}{if TRUE, it returns bootstrap samples of data (Xb) and weights (Wb)}

\item{keepMeasures}{if TRUE, it returns values of all inequality measures for each bootstrap sample}

\item{conf.alpha}{significance level for confidence interval}

\item{calib.boot}{if FALSE, then naive bootstrap is performed, calibrated bootstrap elsewhere}

\item{Xs}{matrix of calibration variables. By default it is a vector of 1's, applied if calib.boot is TRUE}

\item{total}{vector of population totals. By default it is a sum of weights, applied if calib.boot is TRUE}

\item{calib.method}{weights' calibration method for function calib (sampling)}

\item{bounds}{vector of bounds for the g-weights used in the truncated and logit methods; 'low' is the smallest value and 'upp' is the largest value}
}
\value{
This functions returns a data frame from ineq.weighted extended with bootstrap results: expected value, bias (in \%), standard deviation, coefficient of variation, lower and upper bound of confidence interval.
If keepSamples=TRUE or keepMeasures==TRUE then the output becomes a list. If keepSamples=TRUE, the functions returns  Xb and Wb, which are the samples of vector data and the samples of weights, respectively.
If keepMeasures==TRUE, the functions returns Mb, which is a set of inequality measures from bootstrapping.
}
\description{
For weighted mean and weighted total of X (or median of X) as well as for each relevant inequality measure, returns outputs from ineq.weighted and bootstrap outcomes: expected value, bias (in \%), standard deviation, coefficient of variation, lower and upper bound of confidence interval.
}
\details{
By default, naive bootstrap is performed, that is no weights calibration is conducted.
You can choose calibrated bootstrap to calibrate weights with respect to provided variables (Xs) and totals (total).
Confidence interval is simply derived with quantile of order \eqn{\alpha} and \eqn{1-\alpha} where \eqn{\alpha} is a significance level for confidence interval.
}
\examples{
# Inequality measures with additional statistics for numeric variable
X=1:10
W=1:10
ineq.weighted.boot(X,W,B=10)

# Inequality measures with additional statistics for ordered factor variable
X=factor(c('H','H','M','M','L','L'),levels = c('L','M','H'),ordered = TRUE)
W=c(2,2,3,3,8,8)
ineq.weighted.boot(X,W,B=10)

}
