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


draw_rho <-function(detval,epe0,eped,epe0d,n,k,rho,a1=1.01,a2=1.01){

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
	if (idraw > 0 & idraw < nrho) 
		rho = detval[idraw,1]#FIXME: This sometimes fail...

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

#sar_marginal2 <- function(detval,e0,ed,epe0,eped,epe0d,nobs,nvar,a1,a2,c,TI,xs,ys,sige,W){
rho_marginal2 <- function(detval,e0,ed,epe0,eped,epe0d,nobs,nvar, a1,a2,c_beta,TI,xs,ys,sige){#VIRGILIO: W is not needed
n = nrow(detval)
nmk = (nobs-nvar)/2
bprior = dbeta(detval[,1],a1,a2)
C = log(bprior) + lgamma(nmk) - nmk*log(2*pi)
iota = matrix(rep(1,n),n,1)
z = as.numeric(epe0)*iota - 2*detval[,1]*epe0d + detval[,1]*detval[,1]*eped
Q1 = matrix(rep(0,n),n,1);
Q2 = matrix(rep(0,n),n,1);
xpxi = solve(t(xs)%*%xs);
sTI = sige*TI;
xpxis = solve(t(xs)%*%xs + sTI);
logdetx = log(det(xpxis));
C = C - 0.5*logdetx;
          for (i in 1:n){
           rho = detval[i,1];
           D = diag(nobs) - rho*W;	#speye
           bhat = xpxi%*%(t(xs)%*%D%*%ys);
           beta = xpxis%*%(t(xs)%*%D%*%ys + sTI%*%c_beta); 
           Q1[i,1] = t(c_beta - beta)%*%sTI%*%(c_beta - beta);
           Q2[i,1] = t(bhat - beta)%*%(t(xs)%*%xs)%*%(bhat - beta);
          }

den = C + detval[,2] - nmk*log(z + Q1 + Q2);
return(den)
}
#Abhirup: added matmul here

matmul <-function(x,y){
rx = nrow(x)
cx = ncol(x)
ry = nrow(y)
cy = ncol(y)

if((cx == cy) & (rx == ry)){
  out = x*y
}

if ((cx == cy) & (rx == 1)){
	out = y*x[rep(1,ry),]#matrix(rep(x,ry),ry,cy, byrow=TRUE)
}
if ((cx == cy) & (ry == 1)){
	out = x*y[rep(1,rx),]#matrix(rep(y,rx),rx,cx, byrow=TRUE)
}
if ((rx == ry) & (cx == 1)){
	out = y*x[,rep(1,cy)]#matrix(rep(x,cy),ry,cy)
}
if ((rx == ry) & (cy == 1)){
	out = x*y[,rep(1,cx)]#matrix(rep(x,cx),rx,cx)
}
return(out)
}

##code for norm_rnd taken from norm_rnd.m
norm_rnd <-function(sige){
h = chol(sige);
nrow = dim(sige)[1];
x = matrix(rnorm(nrow),nrow,1);
y = t(h)%*%x;
return(y)
}

#code for chis_rnd is written according to its purpose

chis_rnd <-function(nn,v){
if(length(nn)==1) {nrow=nn[1]
out=matrix(rchisq(nrow,df=v),nrow,1)
}
if(length(nn)!=1) {
nrow=nn[1]
ncol=nn[2]
out=matrix(rchisq(nrow*ncol,df=v),nrow,ncol)
}
return(out)
}

sem_marginal<-function(detval,y,x,Wy,Wx,nobs,nvar,a1,a2)
{
	nmk=(nobs-nvar)/2
	nrho=length(detval[,1])
	iota=matrix(rep(1,nrho),nrho,1)
	rvec=detval[,1]
	epe=matrix(rep(0,nrho),nrho,1)
	
	rgrid=as.matrix(seq((detval[1,1]+0.001),(detval[length(detval[,1]),1]-0.001),0.1))
	#rgrid=t(rgrid)
	epetmp=matrix(rep(0,length(rgrid)),length(rgrid),1)
	detxtmp=matrix(rep(0,length(rgrid)),length(rgrid),1)
	
	for(i in 1:length(rgrid))
	{
		xs=x-rgrid[i,1]*Wx
		ys=y-rgrid[i,1]*Wy
		bs=solve((t(xs)%*%xs),t(xs)%*%ys)
		e=ys-xs%*%bs
		epetmp[i,1]=t(e)%*%e
		detxtmp[i,1]=det(t(xs)%*%xs)
	}
	
	tt=rvec
	epe=(spline(x=rgrid,y=epetmp,xout=rvec))$y
	detx=(spline(x=rgrid,y=detxtmp,xout=rvec))$y
	bprior=dbeta(detval[,1],a1,a2)
	C=log(bprior)+log(gamma(nmk))-nmk*log(2*pi)
	den=detval[,2]-0.5*log(detx)-nmk*log(epe)
	out=den+C
	
	return(out)
	
}

sem_marginal2<-function(detval,y,x,Wy,Wx,nobs,nvar,a1,a2,c,TI,sige)
{
	nmk=(nobs-nvar)/2
	nrho=length(detval[,1])
	iota=matrix(rep(1,nrho),nrho,1)
	rvec=detval[,1]
	epe=matrix(rep(0,nrho),nrho,1)
	rgrid=seq((detval[1,1]+0.001),(detval[dim(detval[1]),1]-0.001),0.1)
	rgrid=t(rgrid)
	epetmp=matrix(rep(0,length(rgrid)),length(rgrid),1)
	detxtmp=matrix(rep(0,length(rgrid)),length(rgrid),1)
	Q1=matrix(rep(0,length(rgrid)),length(rgrid),1)
	Q2=matrix(rep(0,length(rgrid)),length(rgrid),1)
	sTI=sige*TI
	for(i in 1:length(rgrid))
	{
		xs=x-rgrid[i,1]*Wx
		ys=y-rgrid[i,1]*Wy
		bs=solve((t(xs)%*%xs),t(xs)%*%ys)
		beta=solve((t(xs)%*%xs+sTI),(t(xs)%*%ys+sTI%*%c))
		e=ys-xs%*%bs
		epetmp[i,1]=t(e)%*%e
		detxtmp[i,1]=det(t(xs)%*%xs)
		Q1[i,1]=t(c-beta)%*%sTI%*%(c-beta)
		Q2[i,1]=t(bs-beta)%*%(t(xs)%*%xs)%*%(c-beta)
	}
	
	tt=rvec
	epe=(spline(x=rgrid,y=epetmp,xout=rvec))$y
	detx=(spline(x=rgrid,y=detxtmp,xout=rvec))$y
	Q1=(spline(x=rgrid,y=Q1,xout=rvec))$y
	Q2=(spline(x=rgrid,y=Q2,xout=rvec))$y
	bprior=dbeta(detval[,1],a1,a2)
	C=log(bprior)+log(gamma(nmk))-nmk*log(2*pi)
	den=detval[,2]-0.5*log(detx)-nmk*log(epe+Q1+Q1)
	out=den+C
	
	return(out)
}

c_rho_sem <- function(rho,y,x,b,sige,W,detval,vi,a1,a2){

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
#	n=dim(index)[1]
#	k=dim(index)[2]
#	nmk=(n-k)/2
#	z=diag(n)-rho*W
#	xs=z%*%x
#	ys=z%*%y
	
#	detx=0.5*log(det(t(xs)%*%xs))
#	n=length(y)
#	e=ys-xs%*%b
#	ev=e*sqrt(vi)
#	epe=nmk%*%log(t(ev)%*%ev)
#	cout=detm-detx-epe

n=length(y)
z=diag(n)-rho*W
e=z%*%y-z%*%x%*%b
ev=e*sqrt(vi)
epe=(t(ev)%*%ev)/(2*sige)
cout=detm-epe
	
	return(cout)

}

draw_rho_sem<-function(detval,y,x,Wy,Wx,V,n,k,rmin,rmax,rho)	
{
	nmk=(n-k)/2
	nrho=length(detval[,1])
	rgrid=seq(rmin+0.01,rmax-0.01,0.01)
	ng=length(rgrid)
	iota=matrix(rep(1,nrho),nrho,1)
	rvec=detval[,1]
	epet=matrix(rep(0,ng),ng,1)
	detxt=matrix(rep(0,ng),ng,1)
	for(i in 1:ng)
	{
		xs=x-rgrid[i]*Wx
		xs=matmul(xs,sqrt(V))
		ys=y-rgrid[i]*Wy
		ys=ys*sqrt(V)
		bs=solve((t(xs)%*%xs),(t(xs)%*%ys))
		e=ys-xs%*%bs
		epet[i,1]=t(e)%*%e
		detxt[i,1]=det(t(xs)%*%xs)
	}
	#Fit on log-scale to avoid having negative values
	epe=exp((spline(x=rgrid,y=log(epet),xout=detval[,1]))$y)
	detx=(spline(x=rgrid,y=detxt,xout=detval[,1]))$y

#Code used for testing
#	if(sum(is.nan(log(epe)))>0)
#	{
#		print(rho)
#		print(min(epe))
#		print(summary(epe))
#		print(sort(epe)[1:5])
#		print(sum(epe<0))
#		print(summary(epet))
#	}

	
	den=detval[,2]-0.5*log(detx)-nmk*log(epe)
	adj=max(den)
	den=den-adj
	den=exp(den)
	
	n=length(den)
	y=detval[,1]
	x=den
	
	isum=sum((y[2:n]+y[1:(n-1)])*(x[2:n]-x[1:(n-1)])/2)
	z=abs(x/isum)
	den=cumsum(z)
	#den=apply(z,2,cumsum)
	
	rnd=runif(1)*sum(z)
	ind=which(den<=rnd)
	idraw=max(ind)
	if((idraw>0) & (idraw<nrho))
	{
		rho=detval[idraw,1]
	}
	return(rho)
}
	
	
	
c_lambda_sac <- function(rho,lambda,y,x,b,sige,W1,W2,detval,P,a1,a2)
{
	gsize = detval[2,1] - detval[1,1]
	# Note these are actually log detvalues
	i1 = which(detval[,1] <= lambda + gsize)
	i2 = which(detval[,1] <= lambda - gsize)
	i1 = max(i1)
	i2 = max(i2)
	index = round((i1+i2)/2)
	#if (length(index)==0)  index = 1
	if (!is.finite(index)) index = 1 #Fixed this

	detm = detval[index,2]

	n=dim(x)[1]
	k=dim(x)[2]
	nmk=(n-k)/2
	B=diag(n)-lambda*W2
	A=diag(n)-rho*W1
	Bx=B%*%x
	
	b=solve((t(Bx)%*%Bx),(t(Bx)%*%B%*%A%*%y))
	
	e=B%*%(A%*%y-x%*%b)
	ev=sqrt(P)*e
	epe=(t(ev)%*%ev)/(2*sige)
	cout=detm-epe
	
	return(cout)
	}

c_rho_sac <- function(rho,lambda,y,x,b,sige,W1,W2,detval,P,a1,a2)	
{
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

	n=dim(x)[1]
	k=dim(x)[2]
	nmk=(n-k)/2
	B=diag(n)-lambda*W2
	A=diag(n)-rho*W1
	Bx=B%*%x
	
	b=solve((t(Bx)%*%Bx),(t(Bx)%*%B%*%A%*%y))
	
	e=B%*%(A%*%y-x%*%b)
	ev=sqrt(P)*e
	epe=(t(ev)%*%ev)/(2*sige)
	cout=detm-epe
	
	return(cout)
	}

matdiv <- function(x,y){
rx=dim(x)[1]
cx=dim(x)[2]

ry=dim(y)[1]
cy=dim(y)[2]

if((cx == cy) & (rx == ry)){
  out = x/y
}
if ((cx == cy) & (rx == 1)){
	out = y/matrix(rep(x,ry),ry,cy)
}
if ((cx == cy) & (ry == 1)){
	out = x/matrix(rep(y,rx),rx,cx)
}
if ((rx == ry) & (cx == 1)){
	out = y/matrix(rep(x,cy),ry,cy)
}
if ((rx == ry) & (cy == 1)){
	out = x/matrix(rep(x,cx),rx,cx)
}
return(out)
}
	
tnorm_rnd <-function(n,amu,sigma,a,b,la,lb,d,kstep){
niter=10
z=matrix(0,n,1)
dinv=solve(d)
anu=d%*%amu
tau=d%*%sigma%*%t(d)
tauinv=solve(tau)
a1=a-anu
b1=b-anu
c=matrix(0,n,n)
h=matrix(0,n,1)
for(i in 1:n){
	aa=tauinv[i,i]
	h[i,1]=1/sqrt(aa)
	for(j in 1:n){
		c[i,j]=-tauinv[i,j]/aa
		}
	}
for(initer in 1:niter){
	for(i1 in 1:n){
		i=kstep[i1,1]
		aa=0
		}
	for(j in 1:n){
		if(i != j)	aa=aa+c[i,j]*z[j,1]	
		}
		
	if(la[i,1]==1)	t1=normrt_rnd(0,1,(b1[i,1]-aa)/h[i,1])
	if(lb[i,1]==1)	t1=normlt_rnd(0,1,(a1[i,1]-aa)/h[i,1])
	if(la[i,1]!=1&lb[i,1]!=1)	t1=normt_rnd(0,1,(a1[i,1]-aa)/h[i,1],(b1[i,1]-aa)/h[i,1])
	z[i,1]=aa+h[i,1]*t1
	}
xdraw=dinv%*%z
for(i in 1:n){
	xdraw[i,1]=xdraw[i,1]+amu[i,1]
	}
xdraw		
}	
	
normrt_rnd <-function(mu,sigma2,right){
nobs=length(mu)
left=-999*matrix(1,nobs,1)
result=normt_rnd(mu,sigma2,left,right)
return(result)
}

normlt_rnd <- function(mu,sigma2,left){
nobs=length(mu)
right=-999*matrix(1,nobs,1)
result=normt_rnd(mu,sigma2,left,right)
return(result)
}	
	
normt_rnd <-function(mu,sigma2,left,right){

std=sqrt(sigma2)
lowerProb=pnorm((left-mu)/std)
upperProb=pnorm((right-mu)/std)
u=runif(length(mu),lowerProb,upperProb) ##??
result=mu+qnorm(u)*std
return(result)
}	
	
###############spdiags################
	
## it works only for square matrices 
## it could work with sparse matrices but it spits a tedious warning 
## it is definitely inefficient compared to the original matlab code 

## choose below different matrices to test the function. 
# r = c(2,3,5,5); c = c(2,1,4,5) 
# A = sparseMatrix(r, c) 
# A = replicate(1000, rnorm(1000) ) 
# A = rbind(c(1,2,3),c(2,3,4),c(3,4,5)) 

spdiags = function(A){ 

     # Find all nonzero diagonals 
     i = rep(seq(1, nrow(A),1),nrow(A)); 
     j = sort(i); 
     d = sort(j-i); 

       # d = d(find(diff([-inf; d]))); ## from Matlab ... 
       # d = c(d[which(diff(d) == 1)], d[length(d)] ) ## this emulate above but needs to stick in last element 

     d = unique(d); ##this should work just fine and it is simpler 
     p = length(d); ##the nr. col of the new matrix 
     m = nrow(A); n = ncol(A); 

     B = matrix(0, nrow = min(c(m,n)), ncol = p); 

   for (k in 1:p){ 
      # print(k) 
      cl = vector(); 

       if (m >= n){ 
          i = max(1, 1+d[k]):min(n, m+d[k]); 
       } else { i = max(1, 1-d[k]):min(m,n-d[k]); } 

       system.time( 
       if (length(i) != 0){ 
          B[i,k] = A[ col(A) == row (A) - d[k]] 
       } ) 
} 

return (list( B = B, d = d) ) 

} 
	
c_sar <-function(rho,y,xb,sige,W,detval,c_beta,T){
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

	z=diag(n)-rho*W
	e=z%*%y-xb
	n=length(y)
	T=T*sige
	z=(diag(n)-rho*W)%*%e
	epe=((t(z)%*%z)/2*sige)+0.5*((rho-c_beta)^2)/T
	
	cout=detm-epe
	cout
}	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	









