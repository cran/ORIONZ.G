\name{ORIONz}
\alias{ORIONz}
\title{
EAP scoring in exploratory FA solutions with correlated residuals
}
\description{
\code{ORIONz} implements procedures for obtaining Bayes EAP: a) individual score estimates; b) conditional standard errors; c) conditional and marginal reliability estimates derived from linear and non-linear extended Exploratoy Factor Analysis solutions that include a correlated-residual structure.
}
\usage{
ORIONz(X, LAM, PHI, THRES, Ree, PSI, model = "linear", disp = TRUE)
}
\arguments{
  \item{X}{
Raw data matrix (dimension: respondents x items).
}
  \item{LAM}{
A matrix containing the item loading estimates in the factor-analytic (FA) solution (dimension: items x factors).
}
  \item{PHI}{
Inter-factor correlation matrix (dimension: factors x factors).
}
  \item{THRES}{
A matrix containing the item thresholds (dimension: thresholds x items). Only required when a non-linear graded FA solution is requested.
}
  \item{Ree}{
Residual inter-item correlation matrix (dimension: items x items).
}
  \item{PSI}{
Diagonal matrix contanining the residual item standard deviations (dimension: items x items).
}
  \item{model}{
The factor-analytic model used for item calibration, where: 'linear' refers to the standard FA model and 'graded' refers to the non-linear ordered-categorical FA model.
}
  \item{disp}{
Determines if the output will be displayed in the console, TRUE by default. If it is TRUE, the output is printed in console and if it is FALSE, the output is returned silently to the output variable.
}
}
\details{
ORIONz implements procedures for obtaining individual score estimates and accompanying measures of accuracy that are based on extended exploratory factor-analytic oblique linear and nonlinear solutions with a residual correlation structure. In the linear case, ORIONz provides closed-form results that use all the structural information available from the calibration stage (including the inter-factor correlations). In the nonlinear case, ORIONz proposes a simple scaling bias correction that, allows more realistic estimates of the accuracy of the factor score estimates to be obtained. In both types of solutions, the outcome includes: individual score estimates, conditional standard errors, and measures of conditional and marginal reliability. An article describing the foundations and technical details of the ORIONz procedures is in preparation, and hopefully, all this information will be available for the interested practitioners very soon.
}
\value{
\item{th_corrected}{Unbiased Bayes EAP individual factor scores estimates obtained by taking into account the residual correlations.}
\item{reli_corrected}{Only if a non-linear graded solutions is required. Unbiased conditional reliability estimates obtained by taking into account the residual correlations. Because residual correlations are generally positive, these corrected estimates are usually deflated compared to the uncorrected ones.}
\item{marginal_reli_corrected}{Unbiased marginal overall reliability estimates obtained by taking into account the residual correlations. Because residual correlations are generally positive, these corrected estimates are usually deflated compared to the uncorrected ones.}
\item{se_corrected}{Only if a non-linear graded solution is required. Unbiased standard errors of measurement corresponding to each individual score estimate. Strictly speaking, these standard errors are Bayes Posterior Standard Deviations.}
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
