\name{sem.inla}
\alias{sem.inla}
\alias{slm.inla}
\alias{sdm.inla}
\title{
Fit spatial econometrics models with INLA
}
\description{
These functions fit some spatial econometrics models for a given
value of \code{rho} (the spatial autocorrelation parameter).
\code{sem.inla} fits a spatial error model, \code{slm} fits a spatial lag model
and \code{sdm.inla} fits a spatial Durbin model.
}
\usage{

sem.inla(formula, d, W, rho, improve = TRUE, impacts = FALSE, fhyper = NULL, probit = FALSE, ...)
slm.inla(formula, d, W, rho, mmatrix = NULL, improve = TRUE, impacts = FALSE, fhyper = NULL, probit = FALSE, ...)
sdm.inla(formula, d, W, rho, mmatrix = NULL, intercept = TRUE, impacts = FALSE, improve = TRUE, fhyper = NULL, probit = FALSE, ...)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{
Formula with the response variable, the fixed effects and, possibly, other
non-linear effects.
}
  \item{d}{
Data.frame with the data.
}
  \item{W}{
Adjacency matrix.
}
  \item{rho}{
Value of the spatial autocorrelation paramter.
}
  \item{improve}{
Logical. Whether improve model fitting (this may require more computing time).
}
  \item{impacts}{
Logical. Whether impacts are computed.
}
  \item{fhyper}{
Options to be passed to the definition of the hyper-parameters in the
spatial effects.
}
  \item{probit}{
Logical. Whether a probit model is used. Note this is only used
when computing the impacts and that argument \code{family} must
be set accordingly.
}
  \item{\dots}{
Other arguments passed to function \code{inla}.
}
}
\details{
These functions fit a spatial econometrics model with a fixed value
of the spatial autocorrelation parameter \code{rho}. 

In addition, the marginal -log-likelihood is corrected to account for
the variance-covariance matrix of the error term or random effects.

}
\value{
An \code{inla} object.
}
\references{
Bivand et al. (2013).
}
\author{
Virgilio Gómez-Rubio <virgilio.gomez@uclm.es>
}
%%\note{
%%  ~~further notes~~
%%}

\seealso{
\code{\link{leroux.inla}}
}
%%\examples{
%%}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}
