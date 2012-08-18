\name{sart_g}
\alias{sart_g}

\title{Bayesian estimates of the spatial autoregressive tobit model}
\description{
          y = rho*W*y + XB + e, e = N(0,sige*Omega), Omega = inv[(I_n-rho*W)'*(I_n -rho*W)]
          y is a binary 0,1 nx1 vector
          B = N(c,T), 
          1/sige = Gamma(nu,d0), 
          rho = beta(a1,a2) prior 
}
\usage{
result =  sart_g(y,W,ndraw,nomit,prior)
}
\arguments{
\item{y}{nobs x 1 independent variable vector (mean = 0)}
\item{W}{nobs x nobs 1st-order contiguity matrix (standardized)}
\item{ndraw}{# of draws}
\item{nomit}{# of initial draws omitted for burn-in}
\item{prior}{a structure variable for prior information input
            prior$beta   = prior means for beta,   c above (default 0)
            priov$bcov   = prior beta covariance , T above (default 1e+12)
            prior$mhflag = 1 for M-H sampling of rho (default = 0)
            prior$nsteps = # of Gibbs steps for Geweke z procedure
            prior$eig    = 0 for computing eigenvalues of W-matrix
                          (defaults to 1, uses rmin = -1, rmax = 1)
            prior$a1     = parameter for beta(a1,a2) prior on rho see: 'help beta_prior'
            prior$a2     = (default = 1.0, a uniform prior on rmin,rmax) 
            prior$rmin   = (optional) min rho used in sampling (default = -1)
            prior$rmax   = (optional) max rho used in sampling (default = 1)  
            prior$lflag  = 0 for full lndet computation (default = 1, fastest)
                         = 1 for MC approx (fast for large problems)
                         = 2 for Spline approx (medium speed)
            prior$order  = order to use with prior.lflag = 1 option (default = 50)
            prior$iter   = iters to use with prior.lflag = 1 option (default = 30) 
            prior$lndet  = a matrix returned by sar, sar_g, sarp_g, etc.
                           containing log-determinant information to save time
}
}
\value{
a structure:
          results$meth   = 'sart_g'
          results$bdraw  = bhat draws (ndraw-nomit x nvar)
          results$pdraw  = rho  draws (ndraw-nomit x 1)
          results$total    = a 3-d matrix (ndraw,nvars-1,ntrs) total x-impacts
          results$direct   = a 3-d matrix (ndraw,nvars-1,ntrs) direct x-impacts
          results$indirect = a 3-d matrix (ndraw,nvars-1,ntrs) indirect x-impacts
                             ntrs defaults to 101 trace terms
          results$nsteps   = nsteps from input
          results$nobs   = # of observations
          results$nvar   = # of variables in x-matrix
          results$ndraw  = # of draws
          results$nomit  = # of initial draws omitted
          results$y      = y-vector from input (nobs x 1)
          results$yhat   = mean of predicted based on posterior parameter means
                           i.e., yhat = inv(I_n - rho*W)*x*beta (nobs x 1)
          results$yprob  = stdn_cdf(yhat)
          results$zip    = # of zero y-values
          results$rmax   = 1/max eigenvalue of W (or rmax if input)
          results$rmin   = 1/min eigenvalue of W (or rmin if input)          
          results$tflag  = 'plevel' (default) for printing p-levels
                         = 'tstat' for printing bogus t-statistics 
          results$lflag  = lflag from input
          results$iter   = prior.iter option from input
          results$order  = prior.order option from input
          results$wlimit  = matrix of [rho lower95,logdet approx, upper95] 
                           intervals for the case of lflag = 1
          results$lndet = a matrix containing log-determinant information
                          (for use in later function calls to save time)
          results$priorb= a flag for diffuse or informative prior on b
          results$acc    = an ndraw x 1 vector of acceptance rates for M-H sampling
          results$cflag  = 1 if there is an intercept, 0 otherwise
}
\source{
Spatial Econometrics Toolbox, Matlab
}
\reference{
LeSage and Pace (2009) 
}

\seealso {
}
\examples{
see test_sart
}
\keyword{Spatial,Econometrics,Bayesian}
