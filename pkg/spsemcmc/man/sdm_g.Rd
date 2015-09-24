\name{sdm_g}
\alias{sdm_g}

\title{Bayesian estimates of the spatial Durbin model}
\description{
          y = rho*W*y + a + X*b1 + W*X*b2 + e, e = N(0,sige*V), V = diag(v1,v2,...vn) 
          r/vi = ID chi(r)/r, r = Gamma(m,k)
          a,b1,b2 = N(prior.beta,prior.bcov)
          1/sige = Gamma(nu,d0), 
          rho = Uniform(rmin,rmax), or rho = beta(a1,a2); 
          a1 = a2 = 1 produces uniform prior
}
\usage{
result =  sdm_g(y,W,ndraw,nomit,prior)
}
\arguments{
\item{y}{nobs x 1 independent variable vector (mean = 0)}
\item{W}{nobs x nobs 1st-order contiguity matrix (standardized)}
\item{ndraw}{# of draws}
\item{nomit}{# of initial draws omitted for burn-in}
\item{prior}{a structure variable for prior information input
            prior$beta  = prior means for beta,   (default 0)
            prior$bcov  = prior beta covariance , (default eye(nvars)*1e+12)
            prior$rval  = r prior hyperparameter, default = 4
            prior$novi  = 1 turns off sampling for vi, producing homoscedastic model V = eye(nobs)           
            prior$m     = informative Gamma(m,k) prior on r
            prior$k     = (default: not used)
            prior$nu    = informative Gamma(nu,d0) prior on sige
            prior$d0    = default: nu=0,d0=0 (diffuse prior)
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
            prior$logm  = 0 for no log marginal calculation, = 1 for log marginal (default = 1)
}
}
\value{
a structure:
        results$meth   = 'sdm_g'
        results$beta   = posterior mean of bhat based on draws
        results$rho    = posterior mean of rho based on draws
        results$sige   = posterior mean of sige based on draws
        results$beta_std = std deviation of beta draws
        results$sige_std = std deviation of sige draws
        results$rho_std  = std deviation of rho draws
        results$total  = a 3-d matrix (ndraw-nomit,nvars-1,101) total x-impacts
                           (ndraw-nomit,nvars,101) if no constant term in the model
        results$direct = a 3-d matrix (ndraw-nomit,nvars-1,ntrs) direct x-impacts
                           (ndraw-nomit,nvars,101) if no constant term in the model
        results$indirect = a 3-d matrix (ndraw-nomit,nvars-1,ntrs) indirect x-impacts
                           (ndraw-nomit,nvars,101) if no constant term in the model
        results$sigma  = posterior mean of sige based on (e'*e)/(n-k)
        results$bdraw  = bhat draws (ndraw-nomit x nvar)
        results$pdraw  = rho  draws (ndraw-nomit x 1)
        results$sdraw  = sige draws (ndraw-nomit x 1)
        results$vmean  = posterior mean of vi draws (nobs x 1) 
        results$rdraw  = r draws (ndraw-nomit x 1) (if m,k input)
        results$bmean  = b prior means, prior.beta from input
        results$bstd   = b prior std deviations sqrt(diag(prior.bcov))
        results$r      = value of hyperparameter r (if input)
        results$novi   = 1 for prior.novi = 1, 0 for prior.rval input
        results$nobs   = # of observations
        results$nvar   = # of variables in [iota x W*x]
        results$ndraw  = # of draws
        results$nomit  = # of initial draws omitted
        results$cflag  = 0 for no intercept, = 1 for intercept term
        results$p      = # of variables in x-matrix (excluding constant
                           term if used)
        results$y      = y-vector from input (nobs x 1)
        results$yhat   = mean of posterior predicted (nobs x 1) inv(I-rho*W)*x*bhat
        results$resid  = residuals, based on posterior means, y - yhat
        results$rsqr   = r-squared based on posterior means (using conventional formula)
        results$rbar   = adjusted r-squared
        results$nu     = nu prior parameter for gamma(nu,d0) prior on sige noise variance
        results$d0     = d0 prior parameter for gamma(nu,d0) prior on sige noise variance
        results$a1     = a1 parameter for beta prior on rho from input, or default value
        results$a2     = a2 parameter for beta prior on rho from input, or default value 1.01  
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
        results$mlike = log marginal likelihood (a vector ranging over
                        rho values that can be integrated for model comparison)
}
\source{
Spatial Econometrics Toolbox, Matlab
}
\references{
LeSage and Pace (2009) Introduction to Spatial Econometrics
 Chapter 5 on Bayesian spatial regression models.
}

\examples{
#see test_sdm
}
\keyword{Spatial,Econometrics,Bayesian}
