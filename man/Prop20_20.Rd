% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/inequalityMeasure.R
\name{Prop20_20}
\alias{Prop20_20}
\title{Proportion 20:20}
\usage{
Prop20_20(X, W = rep(1, length(X)))
}
\arguments{
\item{X}{is a data vector (numeric or ordered factor)}

\item{W}{is a vector of weights}
}
\value{
The value of 20:20 ratio coefficient.
}
\description{
20:20 ratio - originally the ratio of the total income of the 20\% richest people to the 20\% poorest people.
}
\details{
20:20 ratio is calculated as follows:
\deqn{Prop =\frac{H}{L}}
where \eqn{H} is share of 20\% of the highest values,
\eqn{L} is share of 20\% of the lowest values.
}
\examples{
# Compare weighted and unweighted result
X=1:10
W=1:10
Prop20_20(X)
Prop20_20(X,W)

data(Tourism)
#Prop20_20 proportion for Total expenditure with sample weights
X=Tourism$Total_expenditure
W=Tourism$Sample_weight
Prop20_20(X,W)


}
\references{
Panel Data Econometrics: Theoretical Contributions And Empirical Applications edited by Badi Hani Baltag

Notes on Statistical Sources and Methods - The Equality Trust.
}
