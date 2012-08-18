\name{far_g}
\alias{far_g}

\title{Bayesian estimates for the 1st-order Spatial autoregressive model}
\description{
		y = rho*W*y + e,    e = N(0,sige*V), 
        V = diag(v1,v2,...vn), r/vi = ID chi(r)/r, r = Gamma(m,k)
        sige = gamma(nu,d0)    
        rho = Uniform(rmin,rmax), or rho = beta(a1,a2); 
}
\usage{
result =  far_g(y,W,ndraw,nomit,prior)
}
\arguments{
\item{y}{nobs x 1 independent variable vector (mean = 0)}
\item{W}{nobs x nobs 1st-order contiguity matrix (standardized)}
\item{ndraw}{# of draws}
\item{nomit}{# of initial draws omitted for burn-in}
\item{prior}{a structure variable for prior information input
prior$novi  = 1 turns off sampling for vi, producing homoscedastic model. default = 0            
        prior$nu,   = informative Gamma(nu,d0) prior on sige
        prior$d0      default: nu=0,d0=0 (diffuse prior)
        prior$a1    = parameter for beta(a1,a2) prior on rho (default = 1.01)
        prior$a2    = (default = 1.01) see: 'help beta_prior'
        prior$rval, = r prior hyperparameter, default=4
        prior$m,    = informative Gamma(m,k) prior on r
        prior$k,    = informative Gamma(m,k) prior on r
        prior$eig   = 0 for default rmin = -1,rmax = +1, 1 for eigenvalue calculation of these
        prior$rmin, = (optional) min value of rho to use in sampling
        prior$rmax, = (optional) max value of rho to use in sampling
        prior$lflag = 0 for full lndet computation (default = 1, fastest)
                    = 1 for MC approximation (fast for very large problems)
                    = 2 for Spline approximation (medium speed)
        prior$dflag = 0 for numerical integration, 1 for Metropolis-Hastings (default = 0)
        prior$order = order to use with info.lflag = 1 option (default = 50)
        prior$iter  = iters to use with info.lflag = 1 option (default = 30)   
        prior$lndet = a matrix returned by sar, far_g, sar_g, etc.
                      containing log-determinant information to save time
}
}
\value{
a structure:
          results$meth   = 'far_g'
          results$pdraw  = rho draws (ndraw-nomit x 1)
          results$sdraw  = sige draws (ndraw-nomit x 1)
          results$vmean  = mean of vi draws (nobs x 1)
          results$rdraw  = r-value draws (ndraw-nomit x 1)
          results$nu     = prior nu-value for sige (if prior input)
          results$d0     = prior d0-value for sige (if prior input)
          results$a1     = a1 parameter for beta prior on rho from input, or default value
          results$a2     = a2 parameter for beta prior on rho from input, or default value
          results$r      = value of hyperparameter r (if input)
          results$m      = m prior parameter (if input)
          results$k      = k prior parameter (if input)    
          results$nobs   = # of observations
          results$ndraw  = # of draws
          results$nomit  = # of initial draws omitted
          results$y      = actual observations
          results$yhat   = mean of posterior for y-predicted (nobs x 1)
          results$time   = total time taken
          results$time1  = time for log determinant calcluation
          results$time2  = time for eigenvalue calculation   
          results$time3  = time taken for sampling                 
          results$rmax   = 1/max eigenvalue of W (or rmax if input)
          results$rmin   = 1/min eigenvalue of W (or rmin if input) 
          results$tflag  = 'plevel' (default) for printing p-levels
                        = 'tstat' for printing bogus t-statistics      
          results$lflag  = lflag from input
          results$dflag  = dflag value from input (or default value used)
          results$iter   = info.iter option from input
          results$order  = info.order option from input
          results$limit  = matrix of [rho lower95,logdet approx, upper95] intervals
                           (for the case of lflag = 1)      
         results$lndet = a matrix containing log-determinant information
                          (for use in later function calls to save time)
          results$acc   = acceptance rate for M-H sampling
          results$mlike = log marginal likelihood (a vector ranging over
                          rho values that can be integrated for model comparison)
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
see test_far
}
\keyword{Spatial,Econometrics,Bayesian}
