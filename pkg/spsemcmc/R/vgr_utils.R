#
#Conditional distribution of rho given sigma^2(sige)
#
#rho: Value of spatial autocorrelation
#y: Vector of observed values
#sige: sigma^2
#W: Adjacency matrix
#detval:nx2 matrix with values of rho and log-det-Jacobian 
#vi: diag(v_1, ..., v_n)=V
#a1, a2: Parameters from beta prior for rho
#


#Virgilio: Name changed from c_far to c_rho as this is used by other models
#c_far <- function(rho,y,sige,W,detval,vi,a1,a2){
c_rho <- function(rho,y,sige,W,detval,vi,a1,a2){

	gsize = detval[2,1] - detval[1,1]
	# Note these are actually log detvalues
	i1 = which(detval[,1] <= rho + gsize)
	i2 = which(detval[,1] <= rho - gsize)
	i1 = max(i1)
	i2 = max(i2)
	index = round((i1+i2)/2)
	#if (length(index)==0)  index = 1
	if (!is.finite(index)) index = 1 #Fixed this

	detm = detval[index,2]

	n = length(y)
	e = (diag(rep(1,n)) - rho*W)%*%y 
	ev = e*sqrt(vi)
	epe = (t(ev)%*%ev)/(2*sige)
	#bprior = dbeta(detval[,1],a1,a2);
	bprior = dbeta(detval[index,1],a1,a2);#VIRGILIO: Changed this
	epe = log(epe) + log(bprior);

	cout =   detm - (n/2)*epe;
	cout
}


#
#Draws rho-values using Greedy Gibbs and inversion
#
#detval:
#epe0:
#eped:
#epe0d:
#n:
#k:
#rho: spatial autocorrelation
#a1, a2 parameters for beta prior on rho


draw_rho <-function(detval,epe0,eped,epe0d,n,k,rho,a1,a2){

	nmk = (n-k)/2
	nrho = length(detval[,1])
	iota = matrix(rep(1,nrho),nrho,1);

	z = epe0[1,1]*iota - 2*detval[,1]*epe0d[1,1] + detval[,1]*detval[,1]*eped[1,1];
	den = detval[,2] - nmk*log(z)
	bprior = dbeta(detval[,1],a1,a2)#VIRGILIO: Changed this
	den = den + log(bprior)

	n = length(den)
	y = detval[,1]
	adj = max(den)
	den = den - adj
	x = exp(den)

	## trapezoid rule
	#isum = sum((y[2:n,1] + y[1:n-1,1])*(x[2:n,1] - x[1:n-1,1])/2)
	isum = sum((y[2:n] + y[1:(n-1)])*(x[2:n] - x[1:(n-1)])/2)#VIRGILIO:FIXED
	z = abs(x/isum)
	den = cumsum(z)

	rnd = runif(1)*sum(z)
	ind = which(den <= rnd)
	idraw = max(ind)
	if (idraw > 0 & idraw < nrho) rho = detval[idraw,1]#FIXME: This sometimes fail...

	return(rho)
}	


#Compute eigenvalues (range) for the weight matrix
#
#eflag: 1=compute; 0=use the ones provided by the user
#W: weight matrix
#rmin: lower bound of eigenvalues
#rmax: upper bound of eigenvalues
#n: Number of observations (Not used right now...)

#As this is used in several functions, I have changed it to set_eigs()
#far_eigs <- function(eflag,W,rmin,rmax,n){
set_eigs <- function(eflag,W,rmin,rmax,n){
	if (eflag == 1){
		lambda = max(eigen(W)$values)  ##sparse to be used ??
		rmin = 1/lambda
		rmax = 1
#		time2 = 0 ##??VIRGILIO: Not used
	}
	#return(list(rmin=rmin,rmax=rmax,time2=time2))
	return(list(rmin=rmin,rmax=rmax))
}


#Compute log-determinant using interpolation
#
#wsw: Weight matrix (?)
#rmin: Minimium value of rho 
#rmax: Maximum value of rho

lndetint <- function(wsw,rmin=0,rmax=1){
	c=wsw;

	n=dim(c)[1];
	s1=diag(n);
#	z=s1-.1*c;
#	p=colamd(z); ##colmad in R ?
	#%this is the symmetric minimum degree ordering
	#Virgilio: No permutation is used so far. FIXME!!


	iter=100;
	alphavec=((1:iter)-1)/iter;
	##%selecting points to use for the interpolation
	alphaselect=c( 10, 20, 40, 50, 60,  70,  80, 85, 90, 95, 96, 97, 98, 99, 100);
	itersub=length(alphaselect);

	detsub=matrix(rep(0,itersub),itersub,1);
	for (i in 1:itersub) {
		alpha=alphavec[alphaselect[i]];
		z=s1-alpha*c;
#		out=lu(z[,p]); ##lu()
		out=expand(lu(z)) ##lu()
		l=out$L
		u=out$U
		#%LU decomposition
		detsub[i]=sum(log(abs(diag(u))));
	}

	#%stores grid and log-determinant for later use
	#%interpolating for finer grid of alpha
	#VIRGILIO: Added interpolation using splines
	interp<-spline(x=alphavec[alphaselect], y=detsub, xout=alphavec)
	out <- list(rho=interp$x, lndet=interp$y)
	#out$lndet = interp1(c(0,t(alphavec(alphaselect)),c(0,detsub),alphavec,'spline'); ## interp1 in R ? and t() on that?
#	out$rho = t(alphavec);
	return(out)
}


#
#Compute log-det using anohter method (sparse matrices?)
#
#W: adjacency matrix
#lmin: lower bound for rho
#lmax: upper bound for rho

lndetfull <- function(W,lmin=0,lmax=1){

	rvec = seq(lmin,lmax,by=0.01);
	##spparms('tight'); ## Abhirup: in R?
	n  = dim(W)[1];

#VIRGILIO: Not used
#	z = diag(n) - 0.1*W; #ommited speye
#	p = colamd(z); ## colamd in R?


	niter = length(rvec); 
	dettmp = matrix(rep(0,2*niter),niter,2);
	for (i in 1:niter) {
		rho = rvec[i];
		z = diag(n) - rho*W; #speye(n) - rho*sparse(W);

		out=expand(lu(z)) ##lu()
		l=out$L
		u=out$U

		#out = lu(z[,p]); ##lu()
		#l=out$l
		#u=out$r
		dettmp[i,1] = rho;
		dettmp[i,2] = sum(log(abs(diag(u))));
	}

	out <- list()
	out$lndet = dettmp[,2];
	out$rho = dettmp[,1];
return(out)
}


#
#Compute log-determinat 
#
#order: Number of moments
#iter: Number of realizations
#wsw: Weight matrix
#rmin: Lower bound for rho
#rmax: Upper bound for rho
#

#FIXME: There is something wrong in here. The test does not match
#the exact values of hte log-det
lndetmc <- function(order,iter,wsw,rmin=0,rmax=1){

	n=nrow(wsw)

	# Exact moments from 1 to oexact
	td=matrix(c(0,sum(sum(wsw^2))/2), ncol=1)
	oexact=nrow(td)

	o=order
	### Stochastic moments

	mavmomi=matrix(0,o,iter)
	for (j in 1:iter){ 
		u=matrix(rnorm(n), ncol=1)
		v=u
		utu=t(u)%*%u
		for (i in 1:o){
			v=wsw%*%(v)
			mavmomi[i,j]=n*((t(u)%*%v)/(i*utu))
		}
	}

	mavmomi[1:oexact,]=td[,rep(1,iter)] #FIXED.

	###averages across iterations
	#avmomi=t(mean(t(mavmomi)))
	avmomi=matrix(apply(mavmomi,1,mean), ncol=1)#Virgilio: Re-Fixed this



	###%alpha matrix

	alpha=seq(rmin,rmax,0.01)
	#valpha=vander(alpha); VIRGILIO: COmpute Vandermonde matrix
	#valpha = vandermonde.matrix(alpha, length(alpha))#From package 'matrixcalc'
	valpha=sapply((length(alpha)-1):0, function(X){alpha^X})

	#valphaf=fliplr(valpha); VIRGILIO: Flip matrix columns left to right
	valphaf = valpha[, ncol(valpha):1]

	alomat=-valphaf[,c(2:(o+1))]  #Abhirup: error: incorrect no of dimentions?

	##%Estimated ln|I-aD| using mixture of exact, stochastic moments
	##%exact from 1 to oexact, stochastic from (oexact+1) to o

	lndetmat=as.vector(alomat%*%avmomi)


	##%standard error computations
	srvs=t(alomat%*%mavmomi)
	sderr=sqrt((apply(srvs*srvs, 2, mean)-apply(srvs, 2, mean)^2)/iter)

	##%lower bound computation
	eps<-.Machine$double.eps#VIRGILIO: Is this right. YES, it is. FIXED!
	fbound=(n*alpha^(o+1))/((o+1)*(1-alpha+eps))#VIRGILIO: Fixed this

	##%confidendence limits, with lower biased downward (more conservative)
	low95=(lndetmat-1.96*sderr-fbound)
	high95=(lndetmat+1.96*sderr)

	##%AR parameter, lower confidence limit, estimated log-det, upper confidence limit
	##% confide=[alpha'  low95 lndetmat high95];
	out <- list()
	out$rho = t(alpha)
	out$lndet = lndetmat
	out$up95 = high95
	out$lo95 = low95
	return (out)
}



#
#Compute log-determinat of |I-rho*W| using different methods
#
#ldetflag: 0=lndetfull; 1=lndetmc; 2=lndetint
#W: adjacency matrix
#rmin: Lower bound for rho
#rmax: Upper bound for rho
#detval: User-supplied 'detval' matrix (currently not used)
#order: 'order' for MC method
#iter: 'iter' for MC method

#VIRGILIO: Name changed to set_lndet as this is used by different methods
#far_lndet <- function(ldetflag=0,W,rmin,rmax,detval,order,iter){
set_lndet <- function(ldetflag=0,W,rmin=0,rmax=1,detval=NULL,order=NULL,iter=NULL){

	tt=seq(rmin,rmax,by=.001)

	if (ldetflag == 0){ ##% no approximation
		#t0 = clock;    
		#t0=0
		out = lndetfull(W,rmin,rmax); ###lndetful to be used
		#time1 = etime(clock,t0);
		#time1=0
#		tt=seq(rmin,rmax,by=.001); ##% interpolate a finer grid
		#outi = interp1(out$rho,out$lndet,t(tt),'spline'); ##Abhirup: interp1 in R??

	}    

	if(ldetflag == 1){ ## % use Pace and Barry, 1999 MC approximation
		##t0 = clock;    
		out = lndetmc(order,iter,W,rmin,rmax)##lndetmc to be used
		#time1 = 0
		#results$limit = c(out$rho, out$lo95, out$lndet, out$up95);	
		#tt=seq(rmin,rmax,by=0.001) ## % interpolate a finer grid
		#outi = interp1(out$rho,out$lndet,t(tt)) ## spline fitting using out.rho,out.lndet,t(tt) ?? VIRGILIO: We'll check this option later
		#detval = c(t(tt), outi)
	}

	if (ldetflag == 2){ ##% use Pace and Barry, 1998 spline interpolation

		#t0 = clock;
		out = lndetint(W,rmin,rmax);
		#time1 = 0
		#tt=seq(rmin,rmax,by=.001); ##% interpolate a finer grid
		#outi = interp1(out$rho,out$lndet,t(tt),'spline'); ##interp1
		#detval = c(t(tt),outi);
	}

	interp<-spline(x=out$rho, y=out$lndet, xout=t(tt))
	detval = cbind(interp$x, interp$y)

#if(ldetflag == -1){ #% the user fed down a detval matrix
#    time1 = 0;
#        #% check to see if this is right
#        if( detval == 0){
#            print('far_g: wrgon lndet input argument');
#        }
#        n1 = dim(detval)[1];
#	n2 = dim(detval)[2];
#        if (n2 ~= 2){
#            error('far_g: wrong sized lndet input argument');
#        elseif n1 == 1
#            error('far_g: wrong sized lndet input argument');
#        end;
#	}
#}     



	#return (list(detval,time1))    
	return (detval)

}


#
#Compute marginal for rho.
#
#detval: 'detval' matrix with rho's and log-det's
#e0:y-x*b0
#ed:Wy-x*bd
#epe0:e0'*e0
#eped:ed'*ed
#epe0d:ed'*e0
#nobs: Number of observations
#nvar: Number of variances (?)
#logdetx: log(det(X%*%X')) only used when the model includes covariates. Set to
#        0 for the FAR model
#a1,a2: parameters of beta prior for rho

#VIRGILIO: Name changed to rho_marginal
#far_marginal<- function(detval,e0,ed,epe0,eped,epe0d,nobs,nvar,logdetx=0,a1,a2){
rho_marginal<- function(detval,e0,ed,epe0,eped,epe0d,nobs,nvar,logdetx=0,a1,a2){
	n = length(detval)
	nmk = (nobs-nvar)/2
	bprior = dbeta(detval[,1],a1,a2)#VIRGILIO: Changed this
	C = log(bprior) + lgamma(nmk) - nmk*log(2*pi)-0.5*logdetx
	iota = matrix(rep(1,n),n,1)
	z = epe0[1,1]*iota - 2*detval[,1]*epe0d[1,1] + detval[,1]*detval[,1]*eped[,1]
	den =  detval[,2] - nmk*log(z)
	##den = real(den)
	out = C + den;

	return(out)

}

