\name{sarp_g}
\alias{sarp_g}

\title{Bayesian estimates of the spatial autoregressive probit model}
\description{
          y = rho*W*y + XB + e, e = N(0,I_n)
          y is a binary 0,1 nx1 vector
          B = N(c,T), 
          1/sige = Gamma(nu,d0), 
          rho = Uniform(rmin,rmax), or rho = beta(a1,a2); 
}
\usage{
result =  sarp_g(y,W,ndraw,nomit,prior)
}
\arguments{
\item{y}{nobs x 1 independent variable vector (mean = 0)}
\item{W}{nobs x nobs 1st-order contiguity matrix (standardized)}
\item{ndraw}{# of draws}
\item{nomit}{# of initial draws omitted for burn-in}
\item{prior}{a structure variable for prior information input
            prior$nsteps = # of samples used by truncated normal Gibbs sampler
            prior$beta  = prior means for beta,   c above (default 0)
            priov$bcov  = prior beta covariance , T above (default 1e+12)
            prior$a1    = parameter for beta(a1,a2) prior on rho see: 'help beta_prior'
            prior$a2    = (default = 1.0, a uniform prior on rmin,rmax) 
            prior$eig   = 0 for default rmin = -1,rmax = +1, 1 for eigenvalue calculation of these
            prior$rmin  = (optional) min rho used in sampling (default = -1)
            prior$rmax  = (optional) max rho used in sampling (default = 1)  
            prior$lflag = 0 for full lndet computation (default = 1, fastest)
                        = 1 for MC approx (fast for large problems)
                        = 2 for Spline approx (medium speed)
            prior$order = order to use with prior.lflag = 1 option (default = 50)
            prior$iter  = iters to use with prior.lflag = 1 option (default = 30) 
            prior$lndet = a matrix returned by sar, sar_g, sarp_g, etc.
                          containing log-determinant information to save time
}
}
\value{
a structure:
          results$meth     = 'sarp_g'
          results$beta     = posterior mean of bhat based on draws
          results$rho      = posterior mean of rho based on draws
          results$sige     = posterior mean of sige based on draws
          results$sigma    = posterior mean of sige based on (e'*e)/(n-k)
          results$bdraw    = bhat draws (ndraw-nomit x nvar)
          results$pdraw    = rho  draws (ndraw-nomit x 1)
          results$sdraw    = sige draws (ndraw-nomit x 1)
          results$total    = a matrix (ndraw,nvars-1) total x-impacts
          results$direct   = a matrix (ndraw,nvars-1) direct x-impacts
          results$indirect = a matrix (ndraw,nvars-1) indirect x-impacts
          results$total_obs= a matrix (ndraw,nvars-1) observation-level total x-impacts
          results$vmean  = mean of vi draws (nobs x 1) 
          results$rdraw  = r draws (ndraw-nomit x 1) (if m,k input)
          results$bmean  = b prior means, prior.beta from input
          results$bstd   = b prior std deviations sqrt(diag(prior.bcov))
          results$novi   = 1 for prior.novi = 1, 0 for prior.rval input
          results$nobs   = # of observations
          results$nvar   = # of variables in x-matrix
          results$ndraw  = # of draws
          results$nomit  = # of initial draws omitted
          results$nsteps = # of samples used by Gibbs sampler for TMVN
          results$y      = y-vector from input (nobs x 1)
          results$zip    = # of zero y-values
          results$yhat   = mean of posterior predicted (nobs x 1)
          results$resid  = residuals, based on posterior means
          results$rsqr   = r-squared based on posterior means
          results$rbar   = adjusted r-squared
          results$a1     = a1 parameter for beta prior on rho from input, or default value
          results$a2     = a2 parameter for beta prior on rho from input, or default value
          results$rmax   = 1/max eigenvalue of W (or rmax if input)
          results$rmin   = 1/min eigenvalue of W (or rmin if input)          
          results$tflag  = 'plevel' (default) for printing p-levels
                         = 'tstat' for printing bogus t-statistics 
          results$lflag  = lflag from input
          results$cflag  = 1 for intercept term, 0 for no intercept term
          results$iter   = prior.iter option from input
          results$order  = prior.order option from input
          results$limit  = matrix of [rho lower95,logdet approx, upper95] 
                           intervals for the case of lflag = 1
           results$lndet = a matrix containing log-determinant information
                           (for use in later function calls to save time)
           results$mlike = log marginal likelihood (a vector ranging over
                           rho values that can be integrated for model comparison)
}
\source{
Spatial Econometrics Toolbox, Matlab
}
\references{
LeSage and Pace (2009) Chapter 10 on Bayesian estimation 
             of spatial probit regression models.
}

\examples{
#see test_sarp
}
\keyword{Spatial,Econometrics,Bayesian}
