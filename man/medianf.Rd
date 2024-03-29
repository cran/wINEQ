% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/inequalityMeasure.R
\name{medianf}
\alias{medianf}
\title{Median of ordered factor or numeric}
\usage{
medianf(X, W = rep(1, length(X)))
}
\arguments{
\item{X}{is a data vector (numeric or ordered factor)}

\item{W}{is a vector of weights}
}
\value{
The median category (number or label) of ordered factor.
}
\description{
Computes median of ordered factor or numeric variable taking into account weights.
}
\details{
Calculates median based on cumulative distribution. Tailored for ordered factors.
}
\examples{
# Compare weighted and unweighted result
X=factor(c('H','H','M','M','L','L'),levels = c('L','M','H'),ordered = TRUE)
W=c(2,2,3,3,8,8)
medianf(X)
medianf(X,W)



}
