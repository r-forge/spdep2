\name{get_trT}
\alias{get_trT}
\alias{lndetcheb_in}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Chebyshev Jacobian approximation for spatial regression}
\description{
  Find the Jacobian term in spatial regression maximum likelihood by Chebyshev approximation.
}
\usage{
get_trT(W, q = 5)
lndetcheb_in(trT, alpha)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{W}{a \code{CSparseMatrix} object representing a row-standardised spatial weights matrix}
  \item{q}{Chebyshev order}
  \item{trT}{object returned by \code{get_trT}}
  \item{alpha}{single spatial autocorrelation coefficient}
}

\value{
  \code{get_trT} returns the pre-calculated components for \code{lndetcheb_in}, which returns the Jacobian for a given value of the spatial coefficient.
}
\references{LeSage J and RK Pace (2009) Introduction to Spatial Eceonometrics. CRC Press, Boca Raton, pp. 105-108.}
\author{Roger Bivand \email{Roger.Bivand@nhh.no}}

\examples{
data(house)
listw <- nb2listw(LO_nb)
rho <- seq(-0.9, 0.99, 0.01)
form <- formula(log(price) ~ age + I(age^2) + I(age^3) + log(lotsize) +
   rooms + log(TLA) + beds + syear)
system.time(res <- ml_env_setup(form, house, listw, verbose=TRUE))
system.time(Cholesky_J <- sapply(rho, function(lambda) J_fn(res, lambda)))
system.time(trT <- get_trT(W=get("dWd", res), q=5))
system.time(Chebyshev_J <- sapply(rho, function(lambda) lndetcheb_in(trT, lambda)))
plot(rho, Cholesky_J, type="l", col="grey50", ylab="Jacobian")
lines(rho, Chebyshev_J, lty=2)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{spatial}

