\name{sem_g}
\alias{sem_g}

\title{Bayesian estimates of the spatial error model}
\description{
          y = XB + u, u = rho*W + e
          e = N(0,sige*V), V = diag(v1,v2,...vn) 
          r/vi = ID chi(r)/r, r = Gamma(m,k)
          B = N(c,T), 
          1/sige = Gamma(nu,d0), 
          rho = Uniform(rmin,rmax), or rho = beta(a1,a2); 
}
\usage{
result =  sem_g(y,W,ndraw,nomit,prior)
}
\arguments{
\item{y}{nobs x 1 independent variable vector (mean = 0)}
\item{W}{nobs x nobs 1st-order contiguity matrix (standardized)}
\item{ndraw}{# of draws}
\item{nomit}{# of initial draws omitted for burn-in}
\item{prior}{a structure variable for prior information input
            prior$beta  = prior means for beta,   c above (default 0)
            priov$bcov  = prior beta covariance , T above (default 1e+12)
            prior$novi  = 1 turns off sampling for vi, producing homoscedastic model            
            prior$rval  = r prior hyperparameter, default=4
            prior$m     = informative Gamma(m,k) prior on r
            prior$k     = (default: not used)
            prior$nu    = informative Gamma(nu,d0) prior on sige
            prior$d0    = default: nu=0,d0=0 (diffuse prior)
            prior$a1    = parameter for beta(a1,a2) prior on rho see: 'help beta_prior'
            prior$a2    = (default = 1.0, a uniform prior on rmin,rmax) 
            prior$rmin  = (optional) min rho used in sampling (default = -1)
            prior$rmax  = (optional) max rho used in sampling (default = +1)  
            prior$eigs  = 0 to compute rmin/rmax using eigenvalues, (1 = don't compute default)
            prior$lflag = 0 for full lndet computation (default = 1, fastest)
                        = 1 for MC approx (fast for large problems)
                        = 2 for Spline approx (medium speed)
            prior$dflag = 1 for Metropolis-Hastings sampling for rho (default)
                        = 0 for griddy gibbs with univariate numerical integration
            prior$order = order to use with prior.lflag = 1 option (default = 50)
            prior$iter  = iters to use with prior.lflag = 1 option (default = 30)   
            prior$lndet = a matrix returned by sar, sar_g, sarp_g, etc.
                          containing log-determinant information to save time
            prior$mlog  = 0 for no log-marginal likelihood, 
                        = 1 for log-marginal likelihood, default = 1
}
}
\value{
a structure:
          results$meth   = 'sem_g'
          results$beta   = posterior mean of bhat
          results$rho    = posterior mean of rho
          results$sige   = posterior mean of sige
          results$bdraw  = bhat draws (ndraw-nomit x nvar)
          results$pdraw  = rho  draws (ndraw-nomit x 1)
          results$sdraw  = sige draws (ndraw-nomit x 1)
          results$vmean  = mean of vi draws (nobs x 1) 
          results$rdraw  = r draws (ndraw-nomit x 1) (if m,k input)
          results$bmean  = b prior means, prior.beta from input
          results$bstd   = b prior std deviations sqrt(diag(prior.bcov))
          results$r      = value of hyperparameter r (if input)
          results$nobs   = # of observations
          results$nvar   = # of variables in x-matrix
          results$ndraw  = # of draws
          results$nomit  = # of initial draws omitted
          results$y      = y-vector from input (nobs x 1)
          results$yhat   = mean of posterior predicted (nobs x 1)
          results$resid  = residuals, based on posterior means
          results$rsqr   = r-squared based on posterior means
          results$rbar   = adjusted r-squared
          results$nu     = nu prior parameter
          results$d0     = d0 prior parameter
          results$a1     = a1 parameter for beta prior on rho from input, or default value
          results$a2     = a2 parameter for beta prior on rho from input, or default value
          results$time1  = time for eigenvalue calculation
          results$time2  = time for log determinant calcluation
          results$time3  = time for sampling
          results$time   = total time taken  
          results$rmax   = 1/max eigenvalue of W (or rmax if input)
          results$rmin   = 1/min eigenvalue of W (or rmin if input)          
          results$tflag  = 'plevel' (default) for printing p-levels
                         = 'tstat' for printing bogus t-statistics 
          results$lflag  = lflag from input
          results$iter   = prior.iter option from input
          results$order  = prior.order option from input
          results$limit  = matrix of [rho lower95,logdet approx, upper95] 
                           intervals for the case of lflag = 1
          results$lndet = a matrix containing log-determinant information
                          (for use in later function calls to save time)
          results$mlike = log marginal likelihood for model comparisons,
                          (a vector ranging over rho-values from rmin to rmax that can be
                          integrated for model comparison)
          results$acc   = acceptance rate for M-H sampling (ndraw x 1) vector
}
\source{
Spatial Econometrics Toolbox, Matlab
}
\reference{
James P. LeSage, `Bayesian Estimation of Spatial Autoregressive
             Models',  International Regional Science Review, 1997 
             Volume 20, number 1\&2, pp. 113-129.
}

\seealso {
}
\examples{
see test_sem
}
\keyword{Spatial,Econometrics,Bayesian}
