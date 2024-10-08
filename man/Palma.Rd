% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/inequalityMeasure.R
\name{Palma}
\alias{Palma}
\title{Palma index}
\usage{
Palma(X, W = rep(1, length(X)))
}
\arguments{
\item{X}{is a data vector (numeric or ordered factor)}

\item{W}{is a vector of weights}
}
\value{
The value of Palma coefficient.
}
\description{
Palma proportion - originally the ratio of the total income of the 10\% richest people to the 40\% poorest people.
}
\details{
Palma index is calculated by the following formula:
\deqn{Palma =\frac{H}{L}}
where \eqn{H} is share of 10\% of the highest values,
\eqn{L} is share of 40\% of the lowest values.
}
\examples{
# Compare weighted and unweighted result
X=1:10
W=1:10
Palma(X)
Palma(X,W)

data(Tourism)
#Palma index for Total expenditure with sample weights
X=Tourism$Total_expenditure
W=Tourism$Sample_weight
Palma(X,W)



}
\references{
Cobham A., Sumner A.: (2013) Putting the Gini Back in the Bottle? 'The Palma' as a Policy-Relevant Measure of Inequality

Palma J. G.: (2011) Homogeneous middles vs. heterogeneous tails, and the end of the ‘Inverted-U’: the share of the rich is what it’s all about
}
