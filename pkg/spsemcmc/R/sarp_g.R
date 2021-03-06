#
sarp_g <- function(y,x,W,ndraw,nomit,prior){

n=nrow(y)
junk = ncol(y)
yin=y
n1=nrow(x)
k=ncol(x)
n2=nrow(W)
n4=ncol(W)

results <-list()
nobsa=n

results$nobs=n
results$nvar=k
results$y=y
results$zip=n-apply(y,2,sum)

pprior=prior_parse(prior,k)
attach(pprior)

n=length(y)
if(sum(x[,1])!=n){
	tst=apply(x, 2, sum)
	ind=which(tst==n)
	if(length(ind)>0){
		print("sarp_g: intercept term must be in first column of the x-matrix")
	}
	if(length(ind)==0){
		cflag=0
		p=ncol(x)
		}
	}
else#if(sum(x[,1]==n))
{
	#else {
	cflag=1
	p=ncol(x)-1
}	

	results$cflag=cflag
	results$p=p

if(n2!=n4)
{
    stop('sarp_g: wrong size weight matrix W');
}
else
{
	if(n1!=n)
	    stop('sarp_g: wrong size weight matrix W');
}

nchk=nrow(y)
junk=ncol(y)
if(nchk!=n)
    stop('sarp_g: wrong size y vector input');

results$order=order
results$iter=iter

out_temp = set_eigs(eflag,W,rmin,rmax,n)
rmin=out_temp$rmin
rmax=out_temp$rmax
detval = set_lndet(ldetflag, W, rmin, rmax, detval, order, iter)


#Pre-calculate traces for the x-impacts calculations
iiter=50
o=100

diag_ests=matrix(0,n,o)

for(iii in 1:iiter){
	u=matrix(rnorm(n),n,1)
	umat=u[,rep(1,o)]
	wumat=matrix(0,n,o)
	wu=u
	wumat[,1]=wu
	for(ii in 2:o){
		wu=W%*%wu
		wumat[,ii]=wu[,1]
		}
	diag_estimates_iii=umat*wumat
	diag_ests=diag_ests+diag_estimates_iii	
	}
estimated_diags=diag_ests/iiter

bsave = matrix(0,ndraw-nomit,k)
psave = matrix(0,ndraw-nomit,1)	
ymean= matrix(0,n,1)
acc_rate= matrix(0,ndraw,1)

total=matrix(0,ndraw-nomit,p)
total_obs=matrix(0,n,p)
direct=matrix(0,ndraw-nomit,p)
indirect=matrix(0,ndraw-nomit,p)

avg_total=matrix(0,p,1)
avg_direct=matrix(0,p,1)
avg_indirect=matrix(0,p,1)

TI = solve(Tbeta);
TIc = TI%*%c_beta;

In = matrix(1,n,1);

Wy=W%*%y
Wadd=W+t(W)
WtW=t(W)%*%W
sige=1

iter=1
xpx=t(x)%*%x
xpy=t(x)%*%y
Wy=W%*%y
xpWy=t(x)%*%Wy

ind1=which(yin==0)
nobs0=length(ind1)
ind2=which(yin==1)
nobs1=length(ind2)

while(iter <= ndraw){
	AI=solve(xpx+sige*TI)
	ys=y-rho*Wy
	b = t(x)%*%ys + sige*TIc;
	b0 = AI%*%b
	bhat = matrix(mvrnorm(1, b0, sige*AI), ncol=1);  
	xb = x%*%bhat; 
	
	b0 = solve((t(x)%*%x),(t(x)%*%y));
	bd = solve((t(x)%*%x),(t(x)%*%Wy));
	e0 = y - x%*%b0;
	ed = Wy - x%*%bd;
	epe0 = t(e0)%*%e0;
	eped = t(ed)%*%ed;
	epe0d = t(ed)%*%e0;
	rho = draw_rho(detval,epe0,eped,epe0d,n,k,rho,a1,a2)
	
	hh= diag(n)-rho*W
	mu=solve(hh,xb)
	tauinv=diag(n)-rho*Wadd+rho*rho*WtW
	
	aa=as.matrix(diag(tauinv)) 
	h=matrix(1,n,1)/sqrt(aa)
	cc=matdiv(-tauinv,aa) #matdiv
	ctilde=cc-diag(diag(cc))
	
	if(iter==1){
		z=matrix(0,n,1)
		}
	
	for(initer in 1:nsample){
			for(i in 1:n){
			aa=ctilde[i,]%*%z
			muuse=(-mu[i,1]-aa)/h[i,1]
			if(yin[i,1]==0)
				t1=rtnorm(1, mean=0, sd=1, -Inf,muuse)#t1=normrt_rnd(0,1,muuse)
			else
				t1=rtnorm(1, 0, 1, muuse, Inf)#t1=normlt_rnd(0,1,muuse)

			z[i,1]=aa+h[i,1]*t1	
			}
		}
	y=mu+z
	Wy=W%*%y
	
	if(iter>nomit){
		bsave[iter-nomit,1:k]=t(bhat)
		psave[iter-nomit,1]=rho
		ymean=ymean+y
		
		rhovec=(rho^c(0:(o-1)))
		if(cflag==1)	
			beff=matrix(bhat[2:nrow(bhat),], ncol=1) #?
		if(cflag==0)
			beff=bhat

		s=solve(hh,diag(n))
		pdfz=matrix(dnorm(mu[,1]), ncol=1) 	#used dnorm in place of stdn_pdf
		for(kk in 1:p){
			avg_direct[kk,1]=t(pdfz)%*%(estimated_diags%*%rhovec*beff[kk,1]/n)
			dd=diag(pdfz[,1])	##spdiags ??
			avg_total[kk,1]=mean(apply(dd%*%s*beff[kk,1],2,sum))
			total_obs[,kk]=total_obs[,kk]+apply(dd%*%s*beff[kk,1],2,sum)
			avg_indirect[kk,1]=avg_total[kk,1]-avg_direct[kk,1]
			}
		total[iter-nomit,]=t(avg_total)
		direct[iter-nomit,]=t(avg_direct)
		indirect[iter-nomit,]=t(avg_indirect)	
		}
	
	iter=iter+1
	}
total_obs=total_obs/(ndraw-nomit)	

#Compute posterior means
bbeta = apply(bsave,2,mean)
rho=mean(psave)
ymean=ymean/(ndraw-nomit)
results$sige=1
sige=1

nobs=nrow(x)
nvar=ncol(x)

Wy=W%*%ymean
AI=solve(t(x)%*%x+sige*TI)
b0=AI%*%(t(x)%*%ymean+sige*TIc)
bd=AI%*%(t(x)%*%Wy+sige*TIc)
e0=ymean-x%*%b0
ed=Wy-x%*%bd
epe0=t(e0)%*%e0
eped=t(ed)%*%ed
epe0d=t(ed)%*%e0

logdetx=log(det(t(x)%*%x+sige*TI))

if(inform_flag==0)	
	mlike=rho_marginal(detval,e0,ed,epe0,eped,epe0d,nobs,nvar,logdetx,a1,a2)
if(inform_flag==1)	
	mlike=rho_marginal2(detval,e0,ed,epe0,eped,epe0d,nobs,nvar,logdetx,a1,a2,c_beta,TI,x,y,sige)#,W)

results$meth  = 'sarp_g';
results$ymean = ymean;
results$total = total;
results$direct = direct;
results$indirect = indirect;
results$total_obs = total_obs;
results$beta = bbeta;
results$rho = rho;
results$bdraw = bsave;
results$pdraw = psave;
results$bmean = cc;
results$bstd  = sqrt(diag(Tbeta));
results$ndraw = ndraw;
results$nomit = nomit;
results$nsteps = nsample;
results$a1 = a1;
results$a2 = a2;
results$tflag = 'plevel';
results$rmax = rmax; 
results$rmin = rmin;
results$lflag = ldetflag;
results$lndet = detval;
results$priorb = inform_flag;
results$mlike = mlike;


detach(pprior)
return(results)
}

























