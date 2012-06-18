\name{mcdet_prep}
\alias{mcdet_prep}
\alias{mcdet_in}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Monte Carlo Jacobian approximation for spatial regression}
\description{
  Find the Jacobian term in spatial regression maximum likelihood by Monte Carlo approximation.
}
\usage{
mcdet_prep(W, p = 16, m = 30)
mcdet_in(alpha, clx)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{W}{a \code{CSparseMatrix} object representing a row-standardised spatial weights matrix}
  \item{p}{number of sample values to be averaged}
  \item{m}{number of lower order traces of \code{W}}
  \item{alpha}{single spatial autocorrelation coefficient}
  \item{clx}{object returned by \code{mcdet_prep}}
}

\value{
    \code{mcdet_prep} returns the pre-calculated components for \code{mcdet_in}, which returns the vector needed to compute the Jacobian for a given value of the spatial coefficient.
}
\references{LeSage J and RK Pace (2009) Introduction to Spatial Eceonometrics. CRC Press, Boca Raton, pp. 96-105.}
\author{Roger Bivand \email{Roger.Bivand@nhh.no}}

\examples{
data(house)
listw <- nb2listw(LO_nb)
rho <- seq(-0.9, 0.99, 0.01)
form <- formula(log(price) ~ age + I(age^2) + I(age^3) + log(lotsize) +
   rooms + log(TLA) + beds + syear)
system.time(res <- ml_env_setup(form, house, listw, verbose=TRUE))
system.time(Cholesky_J <- sapply(rho, function(lambda) J_fn(res, lambda)))
system.time(clx <- mcdet_prep(W=get("dWd", res), p = 16, m = 30))
system.time(MC_J <- sapply(rho, function(lambda) -mean(mcdet_in(lambda, clx))))
plot(rho, Cholesky_J, type="l", col="grey50", ylab="Jacobian")
lines(rho, MC_J, lty=2)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{spatial}

