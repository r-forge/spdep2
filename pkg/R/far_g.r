
##### support functions ####

beta_prior <-function(rvec,a1,a2){
B = beta(a1,a2)
num = (1+rvec)^(a1-1)
num = num*(1-rvec)^(a2-1)
den = 2^(a1+a2-1)
out = (1/B)*num/den
out[1] = 0.001
out[length(out)] = 0.001
out
}

c_far <- function(rho,y,sige,W,detval,vi,a1,a2){
gsize = detval[2,1] - detval[1,1]
# Note these are actually log detvalues
i1 = which(detval[,1] <= rho + gsize)
i2 = which(detval[,1] <= rho - gsize)
i1 = max(i1)
i2 = max(i2)
index = round((i1+i2)/2)
if (length(index)==0)  index = 1

detm = detval[index,2]

n = length(y)
e = (diag(rep(1,n)) - rho*W)%*%y 
ev = e*sqrt(vi)
epe = (t(ev)*ev)/(2*sige)
bprior = beta_prior(detval[,1],a1,a2);
epe = log(epe) + log(bprior);

cout =   detm - (n/2)*epe;
cout
}

draw_rho <-function(detval,epe0,eped,epe0d,n,k,rho,a1,a2){
nmk = (n-k)/2
nrho = length(detval[,1])
iota = matrix(rep(1,nrho),nrho,1);

z = epe0%*%iota - 2*detval[,1]%*%epe0d + detval[,1]*detval[,1]%*%eped;
den = detval[,2] - nmk*log(z)
bprior = beta_prior(detval[,1],a1,a2)
den = den + log(bprior)

n = length(den)
y = detval[,1]
adj = max(den)
den = den - adj
x = exp(den)

## trapezoid rule
isum = sum((y[2:n,1] + y[1:n-1,1])*(x[2:n,1] - x[1:n-1,1])/2)
z = abs(x/isum)
den = cumsum(z)

rnd = runif(1)*sum(z)
ind = which(den <= rnd)
idraw = max(ind)
if (idraw > 0 & idraw < nrho) rho = detval[idraw,1]
}	

far_eigs <- function(eflag,W,rmin,rmax,n){
if (eflag == 1){
lambda = max(eigen(W)$values)  ##sparse to be used ??
rmin = 1/lambda   
rmax = 1
time2 = 0 ##??
}
return(list(rmin=rmin,rmax=rmax,time2=time2))
}

lndetmc <- function(order,iter,wsw,rmin,rmax){
[n,n]=size(wsw);

# Exact moments from 1 to oexact
td=c(0,sum(sum(wsw^2))/2)  ## full in R??
oexact=length(td)

o=order
### Stochastic moments

mavmomi=matrix(rep(0,o*iter),o,iter)
for (j in 1:iter){ 
u=t(rnorm(n))
v=u
utu=t(u)%*%u
for (i in 1:o){
v=wsw%*%v
mavmomi[i,j]=n*((t(u)%*%v)/(i*utu))
}
}

mavmomi[1:oexact,]=td[,rep(1,iter)]

###averages across iterations
avmomi=t(mean(t(mavmomi)))



###%alpha matrix

alpha=seq(rmin,rmax,0.01)
valpha=vander(alpha); ### ??
valphaf=fliplr(valpha); ### ???
alomat=-valphaf[,c(2:(o+1))]

##%Estimated ln|I-aD| using mixture of exact, stochastic moments
##%exact from 1 to oexact, stochastic from (oexact+1) to o

lndetmat=alomat%*%avmomi


##%standard error computations
srvs=t(alomat%*%mavmomi)
sderr=t(sqrt((mean(srvs*srvs)-mean(srvs)^2)/iter))

##%lower bound computation
fbound=t((n*alpha^(o+1))/((o+1)%*%(1-alpha+eps)))

##%confidendence limits, with lower biased downward (more conservative)
low95=(lndetmat-1.96*sderr-fbound)
high95=(lndetmat+1.96*sderr)

##%AR parameter, lower confidence limit, estimated log-det, upper confidence limit
##% confide=[alpha'  low95 lndetmat high95];

out.rho = t(alpha)
out.lndet = lndetmat
out.up95 = high95
out.lo95 = low95
return (list(out.rho,out.lndet,out.up95,out.lo95))
}

far_lndet <- function(ldetflag=0,W,rmin,rmax,detval,order,iter){
if(ldetflag == 1){ ## % use Pace and Barry, 1999 MC approximation
##t0 = clock;    
out = lndetmc(order,iter,W,rmin,rmax)
time1 = 0
#results.limit = [out.rho out.lo95 out.lndet out.up95];
tt=seq(rmin,rmax,0.001) ## % interpolate a finer grid
outi = spline(out$out.rho,out$out.lndet,t(tt)) ## spline fitting using out.rho,out.lndet,t(tt) ??
detval = c(t(tt), outi)
return (list(detval,time1))    

}

far_marginal <- function(detval,e0,ed,epe0,eped,epe0d,nobs,nvar,a1,a2){
n = length(detval)
nmk = (nobs-nvar)/2
bprior = beta_prior(detval[,1],a1,a2)
C = log(bprior) + lgamma(nmk) - nmk*log(2*pi)
iota = matrix(rep(1,n),n,1)
z = epe0*iota - 2*detval[,1]*epe0d + detval[,1]*detval[,1]*eped
den =  detval[,2] - nmk*log(z)
##den = real(den)
out = C + den;
return(out)

}

#################### Start far_g ####
######
far_g <-function(y,W,ndraw,nomit,prior){
results.order = prior$order
results.iter = prior$iter

V = matrix(rep(1,n),n,1)
In = matrix(rep(1,n),n,1)
vi = In

#### Allocate storage for results ####
psave = matrix(rep(0,ndraw-nomit),ndraw-nomit,1)
ssave = matrix(rep(0,ndraw-nomit),ndraw-nomit,1)
vmean = matrix(rep(0,n),n,1)
yhat = matrix(rep(0,n),n,1)
acc_rate = matrix(rep(0,ndraw),ndraw,1)

#### storage for draw on rvalue ####
if (mm!=0) rsave = matrix(rep(0,ndraw-nomit),ndraw-nomit,1)

tmp = far_eigs(eflag,W,rmin,rmax,n)
rmin = tmp$rmin
rmax = tmp$rmax
time2 = tmp$time2

tmp2 = far_lndet(lndetflag, W, rmin, rmax, detval, order, iter)
detval = tmp2$detval
time1 = tmp2$time1

iter=1
time3 = 0

#### The sampler starts ####

Wy=W%*%y
acc = 0
cc = 0.1

#### For novi_flag=0 or Heteroscedastic Model
if(novi_flag==0){ ##novi_flag_if

#start sampling
while(iter<=ndraw){
	#update sige
	nu1=n+nu
	Wys = sqrt(V)*Wy
	e = ys-rho*Wys
	d1 = d0 +t(e)%*%e
	chi = rchisq(1,nu1)
	t2 = chi/d1
	sige = 1/t2
	
	###update vi
	e = y-rho*Wy
	chiv = rchisq(n,rval+1)
	vi = ((e*e/sige)+In*rval)/chiv
	V = In/vi
	ys = y*sqrt(V)
	
	##update rval
	if(mm!=0){
		rval = rgamma(1,mm,1/kk)
		}
	if(metflag==1){
		#metropolis step to get rho update
		rhox = c_far(rho,y,sige,W,detval,In,a1,a2)
		accept = 0
		rho2 = rho+cc*rnorm(1)
		while(accept==0){
			if((rho2 > rmin)&(rho2 < rmax)){ 
				accept=1
			}
			else rho2 = rho+cc*rnorm(1)
		}
		rhoy = c_far(rho2,y,sige,W,detval,In,a1,a2)
		ru = runif(1)
		if((rhoy-rhox)>exp(1)){
			p=1
		}
		else {
			ratio = exp(rhoy-rhox)
			p = min(ratio)
		}
		if(ru < p){
			rho = rho2
			acc = acc+1
		}
		rtmp[iter,1]=rho
		acc_rate[iter,1]=acc/iter

#### update cc based on sd of rho draws

		if(acc_rate[iter,1]<0.4) cc=cc/1.1
		if(acc_rate[iter,1]>0.6) cc=cc*1.1
	}
	if(metflag==0){
	### update rho using numerical integration ##
	e0 = ys
	ed = Wys
	epe0 = t(e0)%*%e0
	eped = t(ed)%*%ed
	epe0d = t(ed)%*%e0
	rho = draw_rho(detval,epe0,eped,epe0d,n,1,rho,a1,a2)
	}	
	if(iter > nomit){
	ssave[iter-nomit,1] = sige
	psave[iter-nomit,1] = rho
	vmean = vmean+vi
	if(mm!=0) rsave[iter-nomit,1] = rval
	}
iter =iter +1
}
##time??

##compute posterior means and log marginal likelihood for return arguments 

rho = mean(psave)	
sigm = mean(ssave)
vmean = vmean/(ndraw-nomit)
V = In/vmean

ys = sqrt(V)*y
Wys = sqrt(V)*Wy
e0 = ys
ed = Wys
epe0 = t(e0)%*%e0
eped = t(ed)%*%ed
epe0d = t(ed)%*%e0		
e = (e0 - rho*ed)
yhat = y-e
sige = (1/n)*t(e0-rho*ed)%*%(e0-rho*ed)
mlike = far_marginal(detval,e0,ed,eped,epe0d,n,1,a1,a2)

####results

results.y = y      
results.nobs = n
results.nvar = 1   
results.meth = 'far_g'
results.pdraw = psave
results.sdraw = ssave
results.vmean = vmean
results.yhat = yhat
results.resid = e
results.tflag = 'plevel'
results.lflag = ldetflag
results.dflag = metflag
results.nobs  = n
results.ndraw = ndraw
results.nomit = nomit
results.y = y
results.nvar = 1
results.mlike = mlike
results.sige = sige
results.rho = rho
results.lndet = detval
results.acc = acc_rate
results.novi = novi_flag

if(mm != 0){
results.rdraw = rsave
results.m     = mm
results.k     = kk
}
else {
results.r     = rval
results.rdraw = 0
}

results.time = etime(clock,timet)
results.time1 = time1
results.time2 = time2
results.time3 = time3
results.lndet = detval
results.rmax = rmax 
results.rmin = rmin

} ##end novi_flag_if

if(novi_flag==1){ ##strt novi_flag_if
nu1 = n + nu 
e = y - rho*Wy
d1 = d0 + t(e)%*%e
chi = rchisq(1,nu1)
t2 = chi/d1
sige = 1/t2

	if(metflag==1){
		#metropolis step to get rho update
		rhox = c_far(rho,y,sige,W,detval,In,a1,a2)
		accept = 0
		rho2 = rho+cc*rnorm(1)
		while(accept==0){
			if((rho2 > rmin)&(rho2 < rmax)){ 
				accept=1
			}
			else rho2 = rho+cc*rnorm(1)
		}
		rhoy = c_far(rho2,y,sige,W,detval,In,a1,a2)
		ru = runif(1)
		if((rhoy-rhox)>exp(1)){
			p=1
		}
		else {
			ratio = exp(rhoy-rhox)
			p = min(ratio)
		}
		if(ru < p){
			rho = rho2
			acc = acc+1
		}
		rtmp[iter,1]=rho
		acc_rate[iter,1]=acc/iter
		#### update cc based on sd of rho draws

		if(acc_rate[iter,1]<0.4) cc=cc/1.1
		if(acc_rate[iter,1]>0.6) cc=cc*1.1
	} ###end of metflag==1

	if(metflag==0){
	### update rho using numerical integration ##
	e0 = ys
	ed = Wys
	epe0 = t(e0)%*%e0
	eped = t(ed)%*%ed
	epe0d = t(ed)%*%e0
	rho = draw_rho(detval,epe0,eped,epe0d,n,1,rho,a1,a2)
	}	
	if(iter > nomit){
	ssave[iter-nomit,1] = sige
	psave[iter-nomit,1] = rho
	
	if(mm!=0) rsave[iter-nomit,1] = rval
	}
iter =iter +1
} ###end of metflag==0

##time??

##compute posterior means and log marginal likelihood for return arguments 

rho = mean(psave)	

e0 = y
ed = Wy
epe0 = t(e0)%*%e0
eped = t(ed)%*%ed
epe0d = t(ed)%*%e0		
e = (e0 - rho*ed)
yhat = y-e
sige = (1/(n-1))*t(e0-rho*ed)%*%(e0-rho*ed)
mlike = far_marginal(detval,e0,ed,eped,epe0d,n,1,a1,a2)

####results

results.y = y      
results.nobs = n
results.nvar = 1   
results.meth = 'far_g'
results.pdraw = psave
results.sdraw = ssave
results.vmean = vmean
results.yhat = yhat
results.resid = e
results.tflag = 'plevel'
results.lflag = ldetflag
results.dflag = metflag
results.nobs  = n
results.ndraw = ndraw
results.nomit = nomit
results.y = y
results.nvar = 1
results.mlike = mlike
results.sige = sige
results.rho = rho
results.lndet = detval
results.acc = acc_rate
results.novi = novi_flag

if(mm != 0){
results.rdraw = rsave
results.m     = mm
results.k     = kk
}
else {
results.r     = rval
results.rdraw = 0
}

results.time = etime(clock,timet)
results.time1 = time1
results.time2 = time2
results.time3 = time3
results.lndet = detval
results.rmax = rmax 
results.rmin = rmin

} ##end novi_flag_if
return (list(results.meth,results.pdraw,results.sdraw,results.vmean,results.rdraw,results.nu,results.d0,results.a1,results.a2,results.r,results.m,results.k,results.nobs,results.ndraw,results.nomit,results.y,results.yhat,results.time,results.time1,results.time2,results.time3,results.rmax,results.rmin,results.tflag,results.lflag,results.dflag,results.iter,results.order,results.limit,results.lndet,results.acc,results.mlike))
}###end of far_g

