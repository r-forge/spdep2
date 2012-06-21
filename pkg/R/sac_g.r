sar_parse <- function(prior,k){
#% set defaults

eflag = 1;    # % default to not computing eigenvalues
ldetflag = 1;  #% default to 1999 Pace and Barry MC determinant approx
#mflag = 1;     #% default to compute log marginal likelihood
order = 50;    #% there are parameters used by the MC det approx
iter = 30;     #% defaults based on Pace and Barry recommendation
rmin = -1;     #% use -1,1 rho interval as default
rmax = 1;
detval1 = 0;    #% just a flag
detval2 = 0;
rho = 0.5;
sige = 1.0;
rval = 4;
mm = 0;
kk = 0;
nu = 0;
d0 = 0;
a1 = 1.0;
a2 = 1.0;
c = matrix(rep(0,k),k,1); #defuse prior for beta
T = diag(k)*1e+12; ## Abhirup: ??
#prior_beta = 0;   #% flag for diffuse prior on beta
novi_flag = 0; #% do vi-estimates
cc1 = 0.2;        ## % initial tuning parameter for M-H sampling
cc2 = 0.2;
inform_flag = 0;

# check with user input
if(length(prior$novi)!=0) novi_flag=prior$novi
if(length(prior$nu)!=0) nu=prior$nu
if(length(prior$d0)!=0) d0=prior$d0
if(length(prior$rval)!=0) rval=prior$rval
if(length(prior$eigs)!=0) eflag=prior$eigs
if(length(prior$a1)!=0) a1=prior$a1
if(length(prior$a2)!=0) a2=prior$a2
if(length(prior$beta)!=0) {
	c=prior$eigs
	inform_flag=1
	}
if(length(prior$m)!=0){ 
	mm=prior$m
	kk=prior$k
	rval=rgamma(1,mm,1/kk)
	}
if(length(prior$bcov)!=0){ 
	T=prior$bcov
	inform_flag=1
	}
if(length(prior$rmin)!=0) {
	rmin=prior$rmin
	
	}
if(length(prior$rmax)!=0) {
	rmin=prior$rmax
	
	}
if(length(prior$lmin)!=0) {
	lmin=prior$lmin
	
	}
if(length(prior$lmax)!=0) {
	lmax=prior$lmax
	
	}
if(length(prior$lndet)!=0){
	detval1 = prior$lndet1;
	detval2 = prior$lndet2;
	ldetflag = -1;
	eflag = 0;
	rmin = detval1; ##detval[1,1]
	nr = length(detval1);
	rmax = detval1[nr]; ###detval[nr,1]
	lmin = detval2
	nl = length(detval2)
	lmax = detval2[nl]
	}
#if(length(prior$lflag)!=0){
#	tst = prior$lflag;
#        if (tst == 0)        ldetflag = 0; 
#        if (tst == 1)        ldetflag = 1; 
#        if (tst == 2)        ldetflag = 2;
#	}
if(length(prior$order)!=0) order=prior$order
if(length(prior$iter)!=0) iter=prior$iter

return (c(nu,d0,rval,rho,lambda,sige,rmin,rmax,lmin,lmax,detval1,detval2,ldetflag,eflag,order,iter,novi_flag,c,T,cc1,cc2,a1,a2,inform_flag))
}

sac_eigs <- function(eflag,W,rmin,rmax,n){
if (eflag == 1){
lambda = max(eigen(W)$values)  ##sparse to be used 
rmin = 1/lambda   
rmax = 1
time2 = 0 
}
}

## reused code block from far_g ##
## lndet functions ##
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

far_lndet <- function(ldetflag=0,W,rmin,rmax,detval,order,iter){

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




## Begining of sac_g ##
sac_g <- function(y,x,W1,W2,ndraw,nomit,prior){
results <-list()
##% error checking on inputs
n=dim(y)[1]
junk = dim(y)[2]
#[n junk] = size(y);
results$y = y;
n1=dim(x)[1]
k=dim(x)[2]
n3=dim(W1)[1]
n4=dim(W1)[2]
n5=dim(W2)[1]
n6=dim(W2)[2]



#if n1 ~= n
#error('sac_g: x-matrix contains wrong # of observations');
#elseif n3 ~= n4
#error('sac_g: W1 matrix is not square');
#elseif n3~= n
#error('sac_g: W1 matrix is not the same size at y,x');
#elseif n5~= n6
#error('sac_g: W2 matrix is not the same size at y,x');
#elseif n5~= n
#error('sac_g: W2 matrix is not the same size at y,x');

#end;

#% check if the user handled the intercept term okay
#    n = length(y);
#    if sum(x(:,1)) ~= n
#    tst = sum(x); % we may have no intercept term
#    ind = find(tst == n); % we do have an intercept term
#     if length(ind) > 0
#     error('sar: intercept term must be in first column of the x-matrix');
#     elseif length(ind) == 0 % case of no intercept term
#     cflag = 0;
#     p = size(x,2);
#     end;
#    elseif sum(x(:,1)) == n % we have an intercept in the right place
#     cflag = 1;
#     p = size(x,2)-1;
#    end;
     
    results$cflag = cflag;
    results$p = p;
    

#if nargin == 6
#    prior.lflag = 1;
#end;

#[nu,d0,rval,rho,lambda,sige,rmin,rmax,lmin,lmax,detval1,detval2,ldetflag,eflag,order,iter,novi_flag,c,T,cc1,cc2,a1,a2,inform_flag] = sac_parse(prior,k);


results$order = order;
results$iter = iter;

#% error checking on prior information inputs
checkk = dim(c)$1;
#if checkk ~= k
#error('sac_g: prior means are wrong');
#elseif junk ~= 1
#error('sac_g: prior means are wrong');
#end;

#[checkk junk] = size(T);
#if checkk ~= k
#error('sac_g: prior bcov is wrong');
#elseif junk ~= k
#error('sac_g: prior bcov is wrong');
#end;

P = matrix(rep(1,n),n,1); 
In = matrix(rep(1,n),n,1); ##; % initial value for V   
ys = y;
          
vmean = matrix(rep(0,n),n,1); 
yhat = matrix(rep(0,n),n,1);


#[rmin,rmax,time1a] = sac_eigs(eflag,W1,rmin,rmax,n);
tmp_eigs=sac_eigs(eflag,W1,rmin,rmax,n);
rmin=tmp_eigs$rmin
rmax=tmp_eigs$rmax
#time1a=tmp_eigs$time

#[lmin,lmax,time1b] = sac_eigs(eflag,W2,lmin,lmax,n);
tmp_eigs=sac_eigs(eflag,W2,lmin,lmax,n);
lmin=tmp_eigs$rmin
lmax=tmp_eigs$rmax
#time1b=tmp_eigs$time
#results.time1 = time1a + time1b;

results$rmin = rmin;
results$rmax = rmax;
results$lmin = lmin;
results$lmax = lmax;

results$lflag = ldetflag;

tmp_lndet = sac_lndet(ldetflag,W1,rmin,rmax,detval1,order,iter);
detval1=tmp_lndet$detval
tmp_lndet = sac_lndet(ldetflag,W2,lmin,lmax,detval2,order,iter);
detval2=tmp_lndet$detval

#results.time2 = time2a+time2b;


##% storage for draws
          bsave = matrix(rep(0,(ndraw-nomit)*k),ndraw-nomit,k)
          psave = matrix(rep(0,(ndraw-nomit)),ndraw-nomit,1)
          lsave = matrix(rep(0,(ndraw-nomit)),ndraw-nomit,1)
          bsave = matrix(rep(0,(ndraw-nomit)),ndraw-nomit,1)
          vmean = matrix(rep(0,n),n,1)
          acc_rate1 = matrix(rep(0,ndraw),ndraw,1)
          acc_rate2 = matrix(rep(0,ndraw),ndraw,1)

##% ====== initializations
##% compute this stuff once to save time
TI = inv(T);
TIc = TI%*%c;
iter = 1;
In = diag(rep(1,n))
V = In;
Wy = W1%*%y;  ##use sparse here
Wx = W1%*%x;
vi = In;
P = In;
V = vi;
In = diag(rep(1,n));   ### use sparse !!

B = diag(rep(1,n)) - lambda*W2;
A = diag(rep(1,n)) - rho*W1;

nu1 = n + 2*nu;   ##% nu1 is the df of the posterior for sige; n is the # of obs; 


##switch (novi_flag) 
  
  
##case{0} % we do heteroscedastic model
if(novi_flag==0){    
##hwait = waitbar(0,'sac\_g: MCMC sampling ...');
#t0 = clock;                  
iter = 1;
acc1 = 0;
acc2 = 0;

   while (iter <= ndraw); % start sampling;
                  
#   % update beta 
   xtil = B%*%x;                            ##% xtil is used with vi, so not a waste creating here
   ytil = B%*%A%*%y;                          ##% ytil is used with vi, so not a waste creating here
   xstar = matmul(P,xtil);                ##% P here is the sqrt of inv of covar matrix, nx1
   ystar = P*ytil;                       ##% P is all ones when assuming homoscedasticity
   Hinv = inv(t(xstar)%*%xstar + sige%*%TI);    ##% Hinv is covariance of beta
   b0 = Hinv*(t(xstar)%*%ystar + sige%*%TIc);   ##% b0 is mean of beta
   bhat = norm_rnd(Hinv) + b0;            ##% bhat is simulated beta; norm_rnd is MVN, mean 0, covar Hinv  ##code for norm_rnd
   xb = x%*%bhat;

            
##   % update sige (here we take vi into account)
   Bx = (diag(rep(1,n)) - lambda*W2)%*%x; ##use sparse
   b = solve((t(Bx)%*%Bx),(t(Bx)%*%B%*%A%*%y));
   e = B%*%(A%*%y - x%*%b);
   ev = P*e;
   d1 = 2*d0 + t(ev)%*%ev;
   chi = chis_rnd(1,nu1); 	##code for chis_rnd
   sige = d1/chi;

##   % update vi (based on e, without taking vi into account)
          chiv = chis_rnd(n,rval+1);
          vi = ((e*e/sige) + In%*%rval)/chiv;
          P = In/vi; 
              

##   % update lambda using metropolis-hastings
##          % numerical integration is too slow here
          xb = x%*%bhat;
          rhox = c_lambda(rho,lambda,y,x,bhat,sige,W1,W2,detval2,P,a1,a2);   	##c_lambda
          accept = 0;
          lambda2 = lambda + cc1*rnorm(1);		##cc1 matrix??
          while (accept == 0){
           if ((lambda2 > lmin) & (lambda2 < lmax))    accept = 1;  
           else  lambda2 = lambda + cc1*rnorm(1);
            
          } 
           rhoy = c_lambda(rho,lambda2,y,x,bhat,sige,W1,W2,detval2,P,a1,a2);
          ru = unif_rnd(1,0,1);			##unif_rnd
          if ((rhoy - rhox) > exp(1))   pp = 1;
          else {         
          ratio = exp(rhoy-rhox);
          pp = min(ratio); 				##ratio??
          }
              if (ru < pp){
              lambda = lambda2;
              acc1 = acc1 + 1;
              }
      acc_rate1[iter,1] = acc1/iter;		##acc_rate1 ??
      ###% update cc based on std of rho draws
       if (acc_rate1(iter,1) < 0.4){
       cc1 = cc1/1.1;
       }
       if (acc_rate1(iter,1) > 0.6){
       cc1 = cc1*1.1;
       }
       
       B = (diag(rep(1,n)) - lambda*W2;

       
     ##% update rho using metropolis-hastings
       ##   % numerical integration is too slow here
          xb = x%*%bhat;
          rhox = c_rho(rho,lambda,y,x,bhat,sige,W1,W2,detval1,P,a1,a2);
          accept = 0;
          rho2 = rho + cc2*rnorm(1);
          while (accept == 0){
           if ((rho2 > lmin) & (rho2 < lmax))   accept = 1;  
           else  rho2 = rho + cc2*randn(1,1);
           rhoy = c_rho(rho2,lambda,y,x,bhat,sige,W1,W2,detval1,P,a1,a2);
          ru = unif_rnd(1,0,1);
          if ((rhoy - rhox) > exp(1))   pp = 1;
          else{  
          ratio = exp(rhoy-rhox);
          pp = min(1,ratio);
          }
              if (ru < pp){
              rho = rho2;
              acc2 = acc2 + 1;
              }
      acc_rate2[iter,1] = acc2/iter;			##acc_rate2
      ##% update cc based on std of rho draws
       if (acc_rate2(iter,1) < 0.4)       cc2 = cc2/1.1;
       
       if (acc_rate2(iter,1) > 0.6)       cc2 = cc2*1.1;
       

       A = (diag(rep(1,n)) - rho*W1;

                                                         
    if (iter > nomit){ ##% if we are past burn-in, save the draws
    bsave[iter-nomit,1:k] = t(bhat);
    ssave[iter-nomit,1] = sige;
    psave[iter-nomit,1] = rho;
    lsave[iter-nomit,1] = lambda;
    vmean = vmean + vi;        
    }
                    

iter = iter + 1; 
##waitbar(iter/ndraw);         
} ####% end of sampling loop
###close(hwait);

##time3 = etime(clock,t0);
#results$time3 = time3;

##case{1} % we do homoscedastic model 
    
if(novi_flag==1){
#hwait = waitbar(0,'sac\_g: MCMC sampling ...');
#t0 = clock;                  
iter = 1;
acc1 = 0;
acc2 = 0;
          while (iter <= ndraw){ 	##% start sampling;
           
   ##% update beta 
   xtil = B%*%x;                          ##% xtil is used with vi, so not a waste creating here
   ytil = B%*%A%*%y;                   ##       % ytil is used with vi, so not a waste creating here
   Hinv = inv(t(xtil)%*%xtil + sige%*%TI); ##     % Hinv is covariance of beta
   b0 = Hinv%*%(t(xtil)%*%ytil + sige%*%TIc);     ##% b0 is mean of beta
   bhat = norm_rnd(sige%*%Hinv) + b0;       ##% bhat is simulated beta; norm_rnd is MVN, mean 0, covar Hinv
   xb = x%*%bhat;

            
  ## % update sige
   Bx = (diag(rep(1,n)) - lambda%*%W2)%*%x;
   b = solve((t(Bx)%*%Bx),(t(Bx)%*%B%*%A%*%y));
   e = B%*%(A%*%y - x%*%b);
   d1 = 2*d0 + t(e)%*%e;
   chi = chis_rnd(1,nu1);
   sige = d1/chi;


  ## % update lambda using metropolis-hastings
         ## % numerical integration is too slow here
          xb = x%*%bhat;
          rhox = c_lambda(rho,lambda,y,x,bhat,sige,W1,W2,detval2,P,a1,a2);
          accept = 0;
          lambda2 = lambda + cc1%*%rnorm(1);
          while (accept == 0){
           if ((lambda2 > lmin) & (lambda2 < lmax))           accept = 1;  
           else           lambda2 = lambda + cc1*randn(1,1);
	   }
          rhoy = c_lambda(rho,lambda2,y,x,bhat,sige,W1,W2,detval2,P,a1,a2);
          ru = unif_rnd(1,0,1);
          if ((rhoy - rhox) > exp(1))          pp = 1;
          else          ratio = exp(rhoy-rhox);
          pp = min(1,ratio);
          
              if (ru < pp){
              lambda = lambda2;
              acc1 = acc1 + 1;
              }
      acc_rate1[iter,1] = acc1/iter;
##      % update cc based on std of rho draws
       if (acc_rate1(iter,1) < 0.4){
       cc1 = cc1/1.1;
       }
       if (acc_rate1(iter,1) > 0.6){
       cc1 = cc1*1.1;
       }
       
       B = diag(rep(1,n)) - lambda*W2;

       
##     % update rho using metropolis-hastings
##          % numerical integration is too slow here
          xb = x%*%bhat;
          rhox = c_rho(rho,lambda,y,x,bhat,sige,W1,W2,detval1,P,a1,a2);
          accept = 0;
          rho2 = rho + cc2*rnorm(1);
          while (accept == 0){
           if ((rho2 > lmin) & (rho2 < lmax)){
           accept = 1}  
           else{
           rho2 = rho + cc2*rnorm(1);
           }
          }
           rhoy = c_rho(rho2,lambda,y,x,bhat,sige,W1,W2,detval1,P,a1,a2);
          ru = unif_rnd(1,0,1);
          if ((rhoy - rhox) > exp(1)){
          pp = 1;
          else{          
          ratio = exp(rhoy-rhox);
          pp = min(1,ratio);
          }}
              if (ru < pp){
              rho = rho2;
              acc2 = acc2 + 1;
              }
      acc_rate2[iter,1] = acc2/iter;
##      % update cc based on std of rho draws
       if (acc_rate2(iter,1) < 0.4){
       cc2 = cc2/1.1;
       }
       if (acc_rate2(iter,1) > 0.6){
       cc2 = cc2*1.1;
       }

       A = (diag(rep(1,n)) - rho*W1;

                         

    if (iter > nomit){ 	## % if we are past burn-in, save the draws
    bsave[iter-nomit,1:k] = t(bhat);
    ssave[iter-nomit,1] = sige;
    psave[iter-nomit,1] = rho;
    lsave[iter-nomit,1] = lambda;
    vmean = vmean + vi;
    }
                    
iter = iter + 1; 
#waitbar(iter/ndraw);         
} ##% end of sampling loop
##close(hwait);

#time3 = etime(clock,t0);
#results.time3 = time3;
#% ===============================================================


#otherwise
#error('sac_g: unrecognized novi_flag value on input');
#% we should never get here

#end; % end of homoscedastic vs. heteroscedastic options

#% calculate effects estimates
        
#t0 = clock; 

#% pre-calculate traces for the x-impacts calculations
uiter=50;
maxorderu=100;
nobs = n;
rv=matrix(rnorm(nobs*uiter),nobs,uiter);
tracew=matrix(rep(0,maxorderu),maxorderu,1);
wjjju=rv;
for jjj=1:maxorderu{
    wjjju=W1%*%wjjju;
    tracew[jjj]=mean(mean(rv%*%wjjju));
    
}

traces=[tracew];
traces[1]=0;
traces[2]=sum(sum(t(W1)*W1))/nobs;#####mark###
trs=c(1,traces);
ntrs=length(trs);
trbig=t(trs);
                 
        if (cflag == 1){
        bdraws = bsave[,2:length(bsave)];}
        if (cflag == 0){
        bdraws = bsave;
        } 
        pdraws = psave;

        ree = 0:ntrs-1;

        rmat = matrix(rep(0,ntrs),1,ntrs);
        total = matrix(rep(0,(ndraw-nomit)*p*ntrs),ndraw-nomit,p,ntrs); ## Three dimentional matrix in R ??
        direct = matrix(rep(0,(ndraw-nomit)*p*ntrs),ndraw-nomit,p,ntrs);
        indirect = matrix(rep(0,(ndraw-nomit)*p*ntrs),ndraw-nomit,p,ntrs);
        
for (i in 1:ndraw-nomit){
    rmat = pdraws[i,1]^ree;
    for (j in 1:p){
            beta = [bdraws[i,j]];
            total[i,j,] = beta[1,1]*rmat;
    direct[i,j,] = (beta%*%trbig)*rmat;
    indirect[i,j,] = total[i,j,] - direct[i,j,];
    }

}

time4=0 	#time4 = etime(clock,t0);
results$time4 = time4;

#% compute posterior means 
vmean = vmean/(ndraw-nomit);
bmean = mean(bsave);
bmean = t(bmean);
rho = mean(psave);
lam = mean(lsave);
sige = mean(ssave);


B = diag(n) - lam*W2;
A = diag(n) - rho*W1;
Bx = (diag(n) - lam*W2)%*%x;
b = solve((t(Bx)%*%Bx),(t(Bx)%*%B%*%A%*%y));
e = B%*%(A%*%y - x%*%b);
ev = P*e;

results$resid = e;
results$yhat = y-e;

nobs = dim(x)[1];
nvar = dim(x)[2];

sigu = t(e)%*%e;
sige = sigu/(nobs-nvar);
ym = y - mean(y);
rsqr1 = sigu;
rsqr2 = t(ym)%*%ym;
rsqr = 1.0 - rsqr1/rsqr2; ##% conventional r-squared
rsqr1 = rsqr1/(nobs-nvar);
rsqr2 = rsqr2/(nobs-1.0);

time=0	#time = etime(clock,timet);

results$meth  = 'sac_g';
results$beta = bmean;
results$rho = rho;
results$lambda = lam;
results$sige = sige;
results$bdraw = bsave;
results$pdraw = psave;
results$ldraw = lsave;
results$sdraw = ssave;
results$vmean = vmean;
results$total = total;
results$direct = direct;
results$indirect = indirect;
results$yhat  = yhat;
results$bmean = c;
results$bstd  = sqrt(diag(T));
results$rsqr  = rsqr;
results$rbar = 1 - (rsqr1/rsqr2); % rbar-squared
results$sige = sige;
results$nobs  = n;
results$nvar  = nvar;
results$ndraw = ndraw;
results$nomit = nomit;
results$time  = time;
results$acc1 = acc_rate1;
results$acc2 = acc_rate2;
results$nu = nu;
results$d0 = d0;
results$a1 = a1;
results$a2 = a2;
results$tflag = 'plevel';
results$novi = novi_flag;
results$lndet1 = detval1;
results$lndet2 = detval2;
results$priorb = inform_flag;

return(results)
} 
## end of sac_g ##



