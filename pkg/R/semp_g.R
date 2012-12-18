semp_g <- function(y,x,W,ndraw,nomit,prior){

n=nrow(y)
junk=ncol(y)
yin=y
results <-list()
results$y=y
n1=nrow(x)
k=ncol(x)
n3=nrow(W)
n4=ncol(W)



if(n1!=n) 
{
	stop('semp_g: x-matrix contains wrong # of observations');
}
else
{
	if(n3!=n4)
	{
		stop('semp_g: W matrix is not square');
	}
	else
	{
		if(n3!=n)
			stop('semp_g: W matrix is not the same size at y,x');
	}
}


pprior=prior_parse(prior,k)
attach(pprior)

results$order=order
results$iter=iter

# error checking on prior information inputs

checkk=nrow(c_beta)
junk = ncol(c_beta)

if(checkk!=k)
{
	stop('semp_g: prior means are wrong');
}
else
{
	if(junk!=1)
		stop('semp_g: prior means are wrong');
}

checkk=nrow(Tbeta)
junk = ncol(Tbeta)

if(checkk!=k)
{
	stop('semp_g: prior bcov are wrong');
}
else
{
	if(junk!=k)
		stop('semp_g: prior bcov are wrong');
}


out_temp = set_eigs(eflag,W,rmin,rmax,n)
rmin=out_temp$rmin
rmax=out_temp$rmax

results$rmin=rmin
results$rmax=rmax
results$lflag=ldetflag

detval = set_lndet(ldetflag, W, rmin, rmax, detval, order, iter)

bsave=matrix(0,ndraw-nomit,k)
psave=matrix(0,ndraw-nomit,1)
ssave=matrix(0,ndraw-nomit,1)
ymean=matrix(0,n,1)
yhat=matrix(0,n,1)
acc_rate=matrix(0,ndraw,1)

TI=solve(Tbeta)
TIc=TI%*%c_beta
iter=1
In=matrix(1,n,1)
Wy=W%*%y
Wx=W%*%x

zipv=which(yin==0)
zipo=which(yin==1)
nzip=length(zipv)
W2diag=	diag(t(W)%*%W)			##spdiags(t(W)%*%W,0) ##spdiags()

iter=1
acc=0

sige=1#This is set by default and not updated because of identif. problems

while(iter <= ndraw){

	#Update beta
	xs=x-rho*Wx
	AI=solve(t(xs)%*%xs+sige*TI)
	ys=y-rho*Wy
	b=t(xs)%*%ys+sige*TIc
	b0=AI%*%b
	bhat=mvrnorm(1, b0, sige*AI)
	
	#Update sige
	#nu1=n+2*nu
	#e=ys-xs%*%bhat
	#d1=2*d0+t(e)%*%e
	#chi=rchisq(1, nu1)
	#sige=as.numeric(d1/chi)
	#sige<-1
	
	#Update z-values
	mu=x%*%bhat
	ymu=y-mu
	dsig=matrix(1,n,1)-rho*rho*W2diag
	yvar=matrix(1,n,1)/dsig
	A=(1/sige)*(diag(n)-rho*W)%*%ymu
	B=t(diag(n)-rho*W)%*%A
	Cy=ymu-yvar*B
	ym=mu+Cy
	
	ind=which(yin==0)
	#y[ind,1]=normrt_rnd(ym[ind,1],yvar[ind,1],0)
	y[ind,1]=rtnorm(1, mean=ym[ind,1], sd=sqrt(yvar[ind,1]),-Inf, 0)
	ind=which(yin==1)
	#y[ind,1]=normlt_rnd(ym[ind,1],yvar[ind,1],0)
	y[ind,1]=rtnorm(1, mean=ym[ind,1], sd=sqrt(yvar[ind,1]),0, Inf)
	
	Wy=W%*%y
	
	xb=x%*%bhat
	rhox=c_rho_sem(rho,y,x,bhat,sige,W,detval,matrix(1,n,1),a1,a2)
	accept=0
	rho2=rho+cc*rnorm(1)
	while(accept==0){
		if((rho2>rmin)&(rho2<rmax))	
			accept=1
		else	
			rho2=rho+cc*rnorm(1)
		}
		
	rhoy=c_rho_sem(rho2,y,x,bhat,sige,W,detval,matrix(1,n,1),a1,a2)
	ru=runif(1)
#	if(exp(rhoy-rhox) > 1)	
#		p=1
#	else {
		ratio=exp(rhoy-rhox)
		p=min(1,ratio)
#		}
	if(ru<p){
		rho=rho2
		acc=acc+1
		}	
	acc_rate[iter,1]=acc/iter
	if(acc_rate[iter,1]<0.4)
		cc=cc/1.1
	if(acc_rate[iter,1]>0.6)
		cc=cc*1.1
	
	if(iter>nomit){
		bsave[iter-nomit,1:k]=t(bhat)
		ssave[iter-nomit,1]=sige
		psave[iter-nomit,1]=rho
		ymean=ymean+y
		}
	iter=iter+1	
	}
ymean=ymean/(ndraw-nomit)
Wy=W%*%ymean
bmean=apply(bsave,2,mean)
#bmean=t(bmean)
rho=apply(psave,2,mean)

nobs=nrow(x)
nvar=ncol(x)

if(inform_flag==0)	
	mlike=sem_marginal(detval,ymean,x,Wy,Wx,nobs,nvar,a1,a2)
else
	mlike=sem_marginal2(detval,ymean,x,Wy,Wx,nobs,nvar,a1,a2,c_beta,TI,sige)

yhat=x%*%bmean
yprob=pnorm(yhat)
y=results$y
n=length(y)
e=ymean-x%*%bmean
sigu=t(e)%*%e
sige=sigu/(nobs-nvar)
ym=ymean-mean(ymean)
rsqr1=sigu
rsqr2=t(ym)%*%ym
rsqr=1.0-rsqr1/rsqr2
rsqr1=rsqr1/(nobs-nvar)
rsqr2=rsqr2/(nobs-1.0)

results$meth  = 'semp_g';
results$beta = bmean;
results$rho = rho;
results$sige = sige;
results$bdraw = bsave;
results$pdraw = psave;
results$sdraw = ssave;
results$yhat  = yhat;
results$yprob = yprob;
results$ymean = ymean;
results$bmean = c;
results$bstd  = sqrt(diag(T));
results$rsqr  = rsqr;
results$rbar = 1 - (rsqr1/rsqr2); 
results$sige = sige;
results$nobs  = n;
results$zip = nzip;
results$nvar  = nvar;
results$ndraw = ndraw;
results$nomit = nomit;
results$acc = acc_rate;
results$dflag = metflag;
results$nu = nu;
results$d0 = d0;
results$a1 = a1;
results$a2 = a2;
results$mlike = mlike;
results$tflag = 'plevel';
results$novi = novi_flag;
results$lndet = detval;
results$priorb = inform_flag;

detach(pprior)
return(results)
}

