\name{ORIONZ.G-package}
\alias{ORIONZ.G-package}
\title{
EAP scoring in exploratory FA solutions with correlated residuals
}
\description{
\code{ORIONZ.G} implements a procedure for obtaining Bayes EAP score estimates based on extended Exploratory Factor Analysis for both linear and nonlinear solutions with a residual structure.
}
\details{
For more information about the methods used, please go to main page \code{\link{ORIONz}}.
}
\value{
\item{\code{\link{ORIONz}}}{EAP scoring in exploratory FA solutions with correlated residuals}
}
\references{
Ferrando, P.J., Navarro-Gonzalez, D., & Lorenzo-Seva, U. (2024). ORION-z: a simple and comprehensive proposal for EAP score estimation in Exploratory FA solutions with correlated residuals (forthcoming).

}
\author{
David Navarro-Gonzalez

Urbano Lorenzo-Seva

Pere Joan Ferrando
}
\examples{
## Example data set, containing the fitted model of a simulated set of 1000 responders and
## 10 items, with a 2 factor confirmatory structure, and containing 2 doublets
## (items 1-2 and 5-6) :

ORIONz(X = dataconti, LAM = LAM, PHI = PHI, Ree = Ree, PSI = PSI)
}
\keyword{ package }
