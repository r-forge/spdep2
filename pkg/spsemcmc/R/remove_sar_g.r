sar_parse <- function(prior,k){
#% set defaults

eflag = 0;    # % default to not computing eigenvalues
ldetflag = 1;  #% default to 1999 Pace and Barry MC determinant approx
mflag = 1;     #% default to compute log marginal likelihood
order = 50;    #% there are parameters used by the MC det approx
iter = 30;     #% defaults based on Pace and Barry recommendation
rmin = -1;     #% use -1,1 rho interval as default
rmax = 1;
detval = 0;    #% just a flag
rho = 0.5;
sige = 1.0;
rval = 4;
mm = 0;
kk = 0;
nu = 0;
d0 = 0;
a1 = 1.01;
a2 = 1.01;
c = matrix(rep(0,k),k,1); #defuse prior for beta
T = diag(k)*1e+12; ## Abhirup: ??
prior_beta = 0;   #% flag for diffuse prior on beta
novi_flag = 0; #% do vi-estimates
inform_flag = 0;

# check with user input
if(length(prior$novi)!=0) novi_flag=prior$novi
if(length(prior$nu)!=0) nu=prior$nu
if(length(prior$d0)!=0) d0=prior$d0
if(length(prior$rval)!=0) rval=prior$rval
if(length(prior$a1)!=0) a1=prior$a1
if(length(prior$a2)!=0) a2=prior$a2
if(length(prior$m)!=0){ 
	mm=prior$m
	kk=prior$k
	rval=rgamma(1,mm,1/kk)
	}
if(length(prior$m)!=0){ 
	c=prior$beta
	inform_flag=1
	}
if(length(prior$rmin)!=0) {
	rmin=prior$rmin
	eflag=0
	}
if(length(prior$rmax)!=0) {
	rmin=prior$rmax
	eflag=0
	}
if(length(prior$lndet)!=0){
	detval = prior$lndet;
	ldetflag = -1;
	eflag = 0;
	rmin = detval; ##detval[1,1]
	nr = length(detval);
	rmax = detval[nr]; ###detval[nr,1]
	}
if(length(prior$lflag)!=0){
	tst = prior$lflag;
        if (tst == 0)        ldetflag = 0; 
        if (tst == 1)        ldetflag = 1; 
        if (tst == 2)        ldetflag = 2;
	}
if(length(prior$order)!=0) order=prior$order
if(length(prior$iter)!=0) iter=prior$iter
if(length(prior$eig)!=0) eflag=prior$eig
return (c(nu,d0,rval,mm,kk,rho,sige,rmin,rmax,detval,ldetflag,eflag,order,iter,novi_flag,c,T,inform_flag,a1,a2))
}

sar_eigs <- function(eflag,W,rmin,rmax,n){
if (eflag == 1){
lambda = max(eigen(W)$values)  ##sparse to be used 
rmin = 1/lambda   
rmax = 1
time2 = 0 
}
}

### functions for supporting lndet ###
## reused code from far_g ##

lndetint <- function(wsw,rmin=0,rmax=1){
c=wsw;

n=dim(c)[1];
s1=diag(n);
z=s1-.1*c;
p=colamd(z); ##colmad in R ?
#%this is the symmetric minimum degree ordering


iter=100;
alphavec=((1:iter)-1)/iter;
##%selecting points to use for the interpolation
alphaselect=c( 10, 20, 40, 50, 60,  70,  80, 85, 90, 95, 96, 97, 98, 99, 100);
itersub=length(alphaselect);

detsub=matrix(rep(0,itersub),itersub,1);
for (i in 1:itersub) {
alpha=alphavec[alphaselect[i]];
z=s1-alpha*c;
out=lu(z[,p]); ##lu()
l=out$l
r=out$r
#%LU decomposition
detsub[i]=sum(log(abs(diag(u))));
}

#%stores grid and log-determinant for later use
#%interpolating for finer grid of alpha
out <- list()
#out$lndet = interp1(c(0,t(alphavec(alphaselect)),c(0,detsub),alphavec,'spline'); ## interp1 in R ? and t() on that?
out$rho = t(alphavec);
return(out)
}

lndetfull <- function(W,lmin,lmax){
rvec = seq(lmin,lmax,by=0.01);
##spparms('tight'); ## Abhirup: in R?
n  = dim(W)[1];
z = diag(n) - 0.1*W; #ommited speye
p = colamd(z); ## colamd in R?
niter = length(rvec); 
dettmp = matrix(rep(0,2*niter),niter,2);
for (i in 1:niter) {
    rho = rvec[i];
    z = diag(n) - rho*W; #speye(n) - rho*sparse(W);
    out = lu(z[,p]); ##lu()
    l=out$l
    u=out$r
    dettmp[i,1] = rho;
    dettmp[i,2] = sum(log(abs(diag(u))));
}
out <- list()
out$lndet = dettmp[,2];
out$rho = dettmp[,1];
return(out)
}

lndetmc <- function(order,iter,wsw,rmin,rmax){
n=nrow(wsw);#VIRGILIO:Fixed this ,##ABHIRUP: Changed this.

# Exact moments from 1 to oexact
td=c(0,sum(sum(wsw^2))/2)  ## full in R??
oexact=length(td)

o=order
### Stochastic moments

mavmomi=matrix(rep(0,o*iter),o,iter)
for (j in 1:iter){ 
u=rnorm(n)
v=u
utu=t(u)%*%u
for (i in 1:o){
v=wsw%*%(v)
mavmomi[i,j]=n*((t(u)%*%v)/(i*utu))
}
}

#mavmomi[1:oexact,]=td[,rep(1,iter)] # ABHIRUP: is this allright?? dimentions doesn't match...

###averages across iterations
avmomi=t(mean(t(mavmomi)))



###%alpha matrix

alpha=seq(rmin,rmax,0.01)
#valpha=vander(alpha); VIRGILIO: COmpute Vandermonde matrix
valpha = vandermonde.matrix(alpha, length(alpha))#From package 'matrixcalc'

#valphaf=fliplr(valpha); VIRGILIO: Flip matrix columns left to right
valphaf = valpha[, ncol(valpha):1]

alomat=-valphaf[,c(2:(o+1))]  #Abhirup: error: incorrect no of dimentions?

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
out <- list()
out$rho = t(alpha)
out$lndet = lndetmat
out$up95 = high95
out$lo95 = low95
return (out)
}

######### lndet
sar_lndet <- function(ldetflag,W,rmin,rmax,detval,order,iter){
if (ldetflag == 0){ ##% no approximation
#t0 = clock;    
t0=0
out = lndetfull(W,rmin,rmax); ###lndetful to be used
#time1 = etime(clock,t0);
time1=0
tt=seq(rmin,rmax,by=.001); ##% interpolate a finer grid
outi = interp1(out$rho,out$lndet,t(tt),'spline'); ##Abhirup: interp1 in R??
detval = c(t(tt),outi)
}    

if(ldetflag == 1){ ## % use Pace and Barry, 1999 MC approximation
##t0 = clock;    
out = lndetmc(order,iter,W,rmin,rmax)##lndetmc to be used
time1 = 0
results$limit = c(out$rho, out$lo95, out$lndet, out$up95);	
tt=seq(rmin,rmax,by=0.001) ## % interpolate a finer grid
outi = interp1(out$rho,out$lndet,t(tt)) ## spline fitting using out.rho,out.lndet,t(tt) ?? VIRGILIO: We'll check this option later
detval = c(t(tt), outi)
}

if (ldetflag == 2){ ##% use Pace and Barry, 1998 spline interpolation

#t0 = clock;
out = lndetint(W,rmin,rmax);
time1 = 0
tt=seq(rmin,rmax,by=.001); ##% interpolate a finer grid
outi = interp1(out$rho,out$lndet,t(tt),'spline'); ##interp1
detval = c(t(tt),outi);
}
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
#}     
return (list(detval,time1))    
}

draw_rho <- function(detval,epe0,eped,epe0d,n,k,rho,a1,a2){
nmk = (n-k)/2
nrho = length(detval[,1])
iota = matrix(rep(1,nrho),nrho,1);

z = epe0%*%iota - 2*detval[,1]%*%epe0d + detval[,1]*detval[,1]%*%eped;
den = detval[,2] - nmk*log(z)
bprior = dbeta(detval[,1],a1,a2)#VIRGILIO: Changed this
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

sar_marginal <- function(detval,e0,ed,epe0,eped,epe0d,nobs,nvar,logdetx,a1,a2){
n = length(detval)
nmk = (nobs-nvar)/2
bprior = dbeta(detval[,1],a1,a2)#VIRGILIO: Changed this
C = log(bprior) + lgamma(nmk) - nmk*log(2*pi)
iota = matrix(rep(1,n),n,1)
z = epe0*iota - 2*detval[,1]*epe0d + detval[,1]*detval[,1]*eped
den =  detval[,2] - nmk*log(z)
##den = real(den)
out = C + den;
return(out)
}

sar_marginal2 <- function(detval,e0,ed,epe0,eped,epe0d,nobs,nvar,a1,a2,c,TI,xs,ys,sige,W){
n = length(detval)
nmk = (nobs-nvar)/2
bprior = dbeta(detval[,1],a1,a2)#VIRGILIO: Changed this
C = log(bprior) + lgamma(nmk) - nmk*log(2*pi)
iota = matrix(rep(1,n),n,1)
z = epe0*iota - 2*detval[,1]*epe0d + detval[,1]*detval[,1]*eped
Q1 = matrix(rep(0,n),n,1);
Q2 = matrix(rep(0,n),n,1);
xpxi = inv(t(xs)%*%xs);
sTI = sige%*%TI;
xpxis = inv(t(xs)%*%xs + sTI);
logdetx = log(det(xpxis));
C = C - 0.5*logdetx;
          for (i in 1:n){
           rho = detval[i,1];
           D = diag(nobs) - rho*W;	#speye
           bhat = xpxi%*%(t(xs)%*%D%*%ys);
           beta = xpxis%*%(t(xs)%*%D%*%ys + sTI%*%c); 
           Q1(i,1) = t(c - beta)%*%sTI%*%(c - beta);
           Q2(i,1) = t(bhat - beta)%*%(t(xs)%*%xs)%*%(bhat - beta);
          }

den = C + detval[,2] - nmk*log(z + Q1 + Q2);
return(den)
}

#### start sar_g #####
sar_g <- function(y,x,W,ndraw,nomit,prior){


#% error checking on inputs
n=dim(y)[1]
junk = dim(y)[2]
n1=dim(x)[1]
k=dim(x)[2]
n2=dim(W)[1]
n4=dim(W)[2]


time1 = 0
time2 = 0
time3 = 0
time4 = 0
results <-list()
results$nobs  = n;
results$nvar  = k;
results$y = y; 
  
#if nargin == 5
#    prior.lflag = 1;
#end;

#[nu,d0,rval,mm,kk,rho,sige,rmin,rmax,detval,ldetflag, ...
#eflag,order,iter,novi_flag,c,T,inform_flag,a1,a2] = sar_parse(prior,k);
    
temp=sar_parse(prior,k)
nu=temp[1]
d0=temp[2]
rval=temp[3]
mm=temp[4]
kk=temp[5]
rho=temp[6]
sige=temp[7]
rmin=temp[8]
rmax=temp[9]
detval=temp[10]
ldetflag=temp[11]
eflag=temp[12]
order=temp[13]
iter=temp[14]
novi_flag=temp[15]
c=temp[16]
T=temp[17]
inform_flag=temp[18]
a1=temp[19]
a2=temp[20]

# check if the user handled the intercept term okay
#    n = ncol(y)
#    if (sum(x[,1]) != n){
#    tst = sum(x); % we may have no intercept term
#    ind = find(tst == n); % we do have an intercept term
#     if length(ind) > 0
#     error('sar_g: intercept term must be in first column of the x-matrix');
#     elseif length(ind) == 0 % case of no intercept term
#     cflag = 0;
#     p = size(x,2);
#     end;
#    elseif sum(x(:,1)) == n % we have an intercept in the right place
#     cflag = 1;
#     p = size(x,2)-1;
#    end;
#     
#    results.cflag = cflag;
#    results.p = p;
    
#    if n1 ~= n2
#    error('sar_g: wrong size weight matrix W');
#    elseif n1 ~= n
#    error('sar_g: wrong size weight matrix W');
#    end;
#    [nchk junk] = size(y);
#    if nchk ~= n
#    error('sar_g: wrong size y vector input');
#    end;
#   

#results.order = order;
#results.iter = iter;

#timet = clock; % start the timer

out_temp = sar_eigs(eflag,W,rmin,rmax,n)
rmin=out_temp$rmin
rmax=out_temp$rmax
time1=out_temp$time1
#results$time1 = time1;

out_temp = sar_lndet(ldetflag,W,rmin,rmax,detval,order,iter)
detval=out_temp$detval
time2=out_temp$time2
#results$time2 = time2;

results$order = order;
results$iter = iter;


#% storage for draws
bsave = matrix(rep(0,(ndraw-nomit)*k),ndraw-nomit,k)
if (mm !=0) rsave = matrix(rep(0,ndraw-nomit),ndraw-nomit,1)
psave = matrix(rep(0,ndraw-nomit),ndraw-nomit,1)
ssave = matrix(rep(0,ndraw-nomit),ndraw-nomit,1)
vmean= matrix(rep(0,n),n,1)


#% ====== initializations
#% compute this stuff once to save time

TI = solve(T);
TIc = TI%*%c;

In = matrix(rep(0,n),n,1);
V = In;
vi = In;
Wy = W%*%y;

#switch novi_flag
    
#case{0} 
##% we do heteroscedastic model    
if (novi_flag==0){

#hwait = waitbar(0,'sar: MCMC sampling ...');

#t0 = clock;                  
iter = 1;
          while (iter <= ndraw){ #% start sampling;
                  
          #% update beta   
          xs = matmul(x,sqrt(V)); ## code for matmul in support functions
          ys = sqrt(V)*y;
          Wys = sqrt(V)*Wy;
          AI = solve((t(xs)%*%xs + sige%*%TI),diag(rep(1,k));         
          yss = ys - rho*Wys;          
          xpy = t(xs)%*%yss;
          b = t(xs)%*%yss + sige%*%TIc;
          b0 = solve((t(xs)%*%xs + sige%*%TI),b);
          bhat = norm_rnd(sige%*%AI) + b0; ##code for norm_rnd ?? check  
          xb = xs%*%bhat; 
                    
          #% update sige
          nu1 = n + 2*nu; 
          e = (yss - xb);
          d1 = 2*d0 + t(e)*e;
          chi = chis_rnd(1,nu1); ##code for chis_rnd
          sige = d1/chi; 

	  #% update vi
          ev = y - rho*Wy - x%*%bhat; 
          #chiv = chis_rnd(n,rval+1);  
          chiv = matrix(rchisq(n,rval+1),n,1); #% Statistics Toolbox function ##why not use this all the time ??
          vi = ((ev*ev/sige) + In%*%rval)/chiv;
          V = In/vi; 
                        
          #% update rval
          if (mm != 0)   rval = gamm_rnd(1,1,mm,kk);  
          
          
      #% we use griddy Gibbs to perform rho-draw
#          b0 = (t(xs)%*%xs + sige%*%TI )\(t(xs)%*%ys + sige%*%TIc);
#          bd = (t(xs)%*%xs + sige%*%TI)\(t(xs)%*%Wys + sige%*%TIc);
          b0 = solve( (t(xs)%*%xs + sige%*%TI ), (t(xs)%*%ys + sige%*%TIc) );
          bd = solve( (t(xs)%*%xs + sige%*%TI), (t(xs)%*%Wys + sige%*%TIc) );
          e0 = ys - xs%*%b0;
          ed = Wys - xs%*%bd;
          epe0 = t(e0)%*%e0;
          eped = t(ed)%*%ed;
          epe0d = t(ed)%*%e0;
          rho = draw_rho(detval,epe0,eped,epe0d,n,k,rho,a1,a2);  #function later

          
    if (iter > nomit){ ## if we are past burn-in, save the draws
        bsave[iter-nomit,1:k] = t(bhat);
        ssave[iter-nomit,1] = sige;
        psave[iter-nomit,1] = rho;
        vmean = vmean + vi; 
    	if (mm != 0)   rsave(iter-nomit,1) = rval
    	         
    }
                    
iter = iter + 1; 
##waitbar(iter/ndraw);         
}   ##% end of sampling loop
##close(hwait);

#time3 = etime(clock,t0);
#results$time3 = time3;

##case{1} % we do homoscedastic model 
    
#hwait = waitbar(0,'sar: MCMC sampling ...');

if(novi_flag==1){

#t0 = clock;                  
iter = 1;
xpx = t(x)%*%x;
xpy = t(x)%*%y;
Wy = W%*%y;
xpWy = t(x)%*%Wy;


          while (iter <= ndraw){ #% start sampling;
                  
          ##% update beta   
          AI = solve((xpx + sige%*%TI),diag(rep(1,k)));        
          ys = y - rho*Wy;          
          b = t(x)*ys + sige%*%TIc;
          b0 = solve((xpx + sige%*%TI),b);
          bhat = norm_rnd(sige%*%AI) + b0;  
          xb = x%*%bhat;
          
          ##% update sige
          nu1 = n + 2*nu; 
          ##%e = e0 - rho*ed;
          e = (ys - xb);
          d1 = 2*d0 + t(e)%*%e;
          chi = chis_rnd(1,nu1);
          sige = d1/chi;
          
          ###% update rho using griddy Gibbs
          AI = solve((xpx + sige%*%TI),diag(rep(1,k)));
          b0 = solve((xpx + sige%*%TI),(xpy + sige%*%TIc));
          bd = solve((xpx + sige%*%TI),(xpWy + sige%*%TIc));
          e0 = y - x%*%b0;
          ed = Wy - x%*%bd;
          epe0 = t(e0)%*%e0;
          eped = t(ed)%*%ed;
          epe0d = t(ed)%*%e0;
          rho = draw_rho(detval,epe0,eped,epe0d,n,k,rho,a1,a2);

    if (iter > nomit){ ###% if we are past burn-in, save the draws
    bsave[iter-nomit,1:k] = t(bhat);
    ssave[iter-nomit,1] = sige;
    psave[iter-nomit,1] = rho;
    vmean = vmean + vi; 
    }
                    
iter = iter + 1; 
##waitbar(iter/ndraw);         
}### % end of sampling loop
##close(hwait);

##time3 = etime(clock,t0);
##results.time3 = time3;


##otherwise
##error('sar_g: unrecognized prior.novi_flag value on input');
##% we should never get here

}### % end of homoscedastic vs. heteroscedastic vs. log-marginal options

##% calculate effects estimates
        
##t0 = clock; 

##% pre-calculate traces for the x-impacts calculations
uiter=50;
maxorderu=100;
nobs = n;
rv=matrix(rnorm(nobs*uiter),nobs,uiter);
tracew=matrix(rep(0,maxorderu),maxorderu,1);
wjjju=rv;
for (jjj=1:maxorderu){
    wjjju=W%*%wjjju;
    tracew[jjj]=mean(mean(rv*wjjju));
    
}

traces=tracew;
traces[1]=0;
traces[2]=sum(sum(t(W)*W))/nobs;
trs=matrix(c(1,traces));
ntrs=length(trs);
trbig=t(trs);
                 
        if (cflag == 1){
        bdraws = bsave[:,2:length(bsave));
        }
	if(cflag == 0){
        bdraws = bsave;
        }
        pdraws = psave;

        ree = 0:ntrs-1;

        rmat = matrix(rep(0,ntrs),1,ntrs);####three dimentional matrix in R??
        total = matrix(rep(0,(ndraw-nomit)*p*ntrs)),ndraw-nomit,p,ntrs);
        direct = zeros(ndraw-nomit,p,ntrs);
        indirect = zeros(ndraw-nomit,p,ntrs);
        
for i=1:ndraw-nomit;
    rmat = pdraws(i,1).^ree;
    for j=1:p;
            beta = [bdraws(i,j)];
            total(i,j,:) = beta(1,1)*rmat;
    direct(i,j,:) = (beta*trbig).*rmat;
    indirect(i,j,:) = total(i,j,:) - direct(i,j,:);
    end;

end;
#### ?? ###
#time4 = etime(clock,t0);
#results.time4 = time4;


#% compute posterior means and log marginal likelihood for return arguments
bmean = mean(bsave);
beta = t(bmean);
rho = mean(psave);
sige = mean(ssave);
vmean = vmean/(ndraw-nomit);
V = In/vmean;


results$sige = sige;
nobs=dim(x)$1
nvar=dim(x)$2
#[nobs,nvar] = size(x);
          xs = matmul(x,sqrt(V));
          ys = sqrt(V)*y;
          Wys = W%*%ys;
          AI = inv(t(xs)%*%xs + sige%*%TI);
          b0 = AI*(t(xs)%*%ys + sige%*%TIc);
          bd = AI*(t(xs)%*%Wys + sige%*%TIc);
          e0 = ys - xs%*%b0;
          ed = Wys - xs%*%bd;
          epe0 = t(e0)%*%e0;
          eped = t(ed)%*%ed;
          epe0d = t(ed)%*%e0;
 logdetx = log(det(t(xs)%*%xs + sige%*%TI));
  if (inform_flag == 0){
   mlike = sar_marginal(detval,e0,ed,epe0,eped,epe0d,nobs,nvar,logdetx,a1,a2);}
  if(inform_flag == 1){
   mlike = sar_marginal2(detval,e0,ed,epe0,eped,epe0d,nobs,nvar,logdetx,a1,a2,c,TI,xs,ys,sige,W);
  }
 yhat = solve((diag(rep(1,nobs)) - rho*W),(x%*%beta));
 e = y - yhat;

##% compute R-squared
epe = t(e)%*%e;
sige = epe/(n-k);
results$sigma = sige;
ym = y - mean(y);
rsqr1 = epe;
rsqr2 = t(ym)%*%ym;
results$rsqr = 1- rsqr1/rsqr2; ##% r-squared
rsqr1 = rsqr1/(nobs-nvar);
rsqr2 = rsqr2/(nobs-1.0);
results$rbar = 1 - (rsqr1/rsqr2); ###% rbar-squared

##time = etime(clock,timet);


results$meth  = 'sar_g';
results$total = total;
results$direct = direct;
results$indirect = indirect;
results$beta_std = std(bsave)';
results$sige_std = std(ssave);
results$rho_std = std(psave);
results$beta = beta;
results$rho = rho;
results$bdraw = bsave;
results$pdraw = psave;
results$sdraw = ssave;
results$mlike = mlike;
results$vmean = vmean;
results$yhat  = yhat;
results$resid = e;
results$bmean = c;
results$bstd  = sqrt(diag(T));
results$ndraw = ndraw;
results$nomit = nomit;
results$time  = time;
results$nu = nu;
results$d0 = d0;
results$a1 = a1;
results$a2 = a2;
results$tflag = 'plevel';
results$rmax = rmax; 
results$rmin = rmin;
results$lflag = ldetflag;
results$lndet = detval;
results$novi  = novi_flag;
results$priorb = inform_flag;

if (mm!= 0){
results$rdraw = rsave;
results$m     = mm;
results$k     = kk;}
if(mm==0){
results$r     = rval;
results$rdraw = 0;
}
}#### end of sar_g


