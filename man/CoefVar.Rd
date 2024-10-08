% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/inequalityMeasure.R
\name{CoefVar}
\alias{CoefVar}
\title{Coefficient of Variation}
\usage{
CoefVar(X, W = rep(1, length(X)), square = FALSE)
}
\arguments{
\item{X}{is a data vector}

\item{W}{is a vector of weights}

\item{square}{logical, argument of the function CoefVar, for details see below}
}
\value{
The value of CoefVar coefficient.
}
\description{
Computes Coefficient of Variation inequality measure of a given variable taking into account weights.
}
\details{
Coefficient of variation is given by:
\deqn{CV= \frac{\sigma}{\mu}\times 100}
where \eqn{\sigma} is a standard deviation and \eqn{\mu} is arithmetic mean.
}
\examples{
# Compare weighted and unweighted result
X=1:10
W=1:10
CoefVar(X)
CoefVar(X,W)

data(Tourism)
#Coefficient of variation for Total expenditure with sample weights
X=Tourism$Total_expenditure
W=Tourism$Sample_weight
CoefVar(X,W)


}
\references{
Sheret M.: (1984) Social Indicators Research, An International and Interdisciplinary Journal for Quality-of-Life Measurement, Vol. 15, No. 3, Oct. ISSN 03038300

Coulter P. B.: (1989) Measuring Inequality ISBN 0-8133-7726-9
}
