% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/inequalityMeasure.R
\name{BL}
\alias{BL}
\title{Blair and Lacy index}
\usage{
BL(X, W = rep(1, length(X)), withsqrt = FALSE)
}
\arguments{
\item{X}{is a data vector (numeric or ordered factor)}

\item{W}{is a vector of weights}

\item{withsqrt}{if TRUE function returns index given by BL2, elsewhere by BL (default). See more in details.}
}
\value{
The value of Blair and Lacy coefficient.
}
\description{
Computes Blair and Lacy inequality measure of a given variable taking into account weights.
}
\details{
Let \eqn{m} be the median category, \eqn{n} be the number of categories and \eqn{P_i} be the cumulative distribution of \eqn{i}-th category.
The indices of Blair and Lacy (2000) are the following:
\deqn{BL =  1-\frac{\sum_{i=1}^{n-1}(P_{i}-0.5)^2}{\frac{n-1}{4}}}
\deqn{BL2 =  1-\left(\frac{\sum_{i=1}^{n-1}(P_{i}-0.5)^2}{\frac{n-1}{4}}\right)^{\frac{1}{2}}}
}
\examples{
# Compare weighted and unweighted result
X=1:10
W=1:10
BL(X)
BL(X,W)

data(Well_being)
# Blair and Lacy index for health assessment with sample weights
X=Well_being$V1
W=Well_being$Weight
BL(X,W)


}
\references{
Blair J, Lacy M G. (2000): Statistics of ordinal variation, Sociological Methods and Research 28(251);251-280.
}
