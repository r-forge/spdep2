\name{ml_env_setup}
\alias{ml_env_setup}
\alias{do_LL}
\alias{sse_fn}
\alias{J_fn}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Error SAR computation testbed functions}
\description{
  These functions are to permit the testing of different components of maximum likelihood estimation for spatial regression models using \pkg{Matrix} class objects and methods at the R and C interface levels.
}
\usage{
ml_env_setup(formula, data, listw, verbose = TRUE)
do_LL(val, env, interp=c(FALSE, FALSE))
sse_fn(env, lambda)
J_fn(env, lambda)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{a symbolic description of the model to be fit. The details 
of model specification are given for \code{lm()}}
  \item{data}{an optional data frame (or data frame like object) containing the variables in the model. By default the variables are taken from the environment which the function is called}
  \item{listw}{a \code{listw} object created for example by \code{nb2listw}}
  \item{verbose}{default=FALSE; if TRUE, reports function values during optimization.}
  \item{lambda,val}{single spatial autocorrelation coefficient}
  \item{env}{an environment object as returned by \code{ml_env_setup}}
  \item{interp}{default \code{c(FALSE, FALSE)}, if the first value is TRUE, \code{sse_fn} is used instead of compiled code, if the second, \code{J_fn} is used instead of compiled code}
}
\details{
\code{ml_env_setup} populates an environment with the objects needed for compting the sum of squared errors and Jacobian terms in an error SAR model.
\code{dl_LL} is an R function passed to \code{optimize} using the environment from \code{ml_env_setup}, and permitting interpreted or compiled functions to be used for the two terms; the two compiled functions use QR methods and updating Cholesky methods respectively.
\code{sse_fn} calculates the SSE term using interpreted code and QR methods.
\code{J_fn} calculates the Jacobian term using interpreted code and updating Cholesky methods.
}
\value{
  \code{ml_env_setup} returns an environment object containing the objects needed for maximum likelihood optimization; \code{dl_LL} retuns the value of the log-likelihood function at \code{val}, \code{sse_fn} and \code{J_fn} are interpreted functions returning the sum of squared errors and Jacobian terms respectively.
}

\author{Roger Bivand \email{Roger.Bivand@nhh.no}}

\examples{
data(house)
listw <- nb2listw(LO_nb)
form <- formula(log(price) ~ age + I(age^2) + I(age^3) + log(lotsize) +
   rooms + log(TLA) + beds + syear)
system.time(res <- ml_env_setup(form, house, listw, verbose=TRUE))
system.time(optimize(do_LL, interval=c(-1, 1), env=res,
 interp=c(TRUE, TRUE), maximum=TRUE, tol=.Machine$double.eps^0.5))
system.time(optimize(do_LL, interval=c(-1, 1), env=res,
 interp=c(FALSE, TRUE), maximum=TRUE, tol=.Machine$double.eps^0.5))
system.time(optimize(do_LL, interval=c(-1, 1), env=res,
 interp=c(TRUE, FALSE), maximum=TRUE, tol=.Machine$double.eps^0.5))
system.time(optimize(do_LL, interval=c(-1, 1), env=res,
 interp=c(FALSE, FALSE), maximum=TRUE, tol=.Machine$double.eps^0.5))
system.time(.Call("do_LL", res, as.double(c(-1, 1)), .Machine$double.eps^0.5,
 PACKAGE="spdep2"))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{spatial}

