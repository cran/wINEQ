% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/inequalityMeasure.R
\name{RicciSchutz}
\alias{RicciSchutz}
\title{Ricci and Schutz index}
\usage{
RicciSchutz(X, W = rep(1, length(X)))
}
\arguments{
\item{X}{is a data vector}

\item{W}{is a vector of weights}
}
\value{
The value of Ricci and Schutz coefficient.
}
\description{
Computes Ricci and Schutz inequality measure of a given variable taking into account weights.
}
\details{
In the case of an empirical distribution with n elements where \eqn{y_{i}} denotes the wealth of household \eqn{i} and \eqn{\overline{y}} the sample average, the Ricci and Schutz coefficient can be expressed as:
\deqn{RS =  \frac{1}{2n} \sum_{i=1}^{n} \frac{\mid y_{i} - \overline{y} \mid}{\overline{y}}}
}
\examples{
# Compare weighted and unweighted result
X=1:10
W=1:10
RicciSchutz(X)
RicciSchutz(X,W)

data(Tourism)
#Ricci and Schutz index for Total expenditure with sample weights
X=Tourism$Total_expenditure
W=Tourism$Sample_weight
RicciSchutz(X,W)


}
\references{
Coulter P. B.: (1989) Measuring Inequality ISBN 0-8133-7726-9

Eliazar I. I., Sokolov I. M.: (2010) Measuring statistical heterogeneity: The Pietra index

Costa R. N., Pérez-Duarte S.: (2019) Not all inequality measures were created equal, Statistics Paper Series, No 31
}
