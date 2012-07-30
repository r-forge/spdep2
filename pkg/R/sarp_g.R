#
#sar_g <- function(y,x,W,ndraw,nomit,prior){

n=dim(y)[1]
junk = dim(y)[2]
yin=y
n1=dim(x)[1]
k=dim(x)[2]
n2=dim(W)[1]
n4=dim(W)[2]

results <-list()
nobsa=n

results$nobs=n
results$nvar=k
results$y=y
results$zip=n-apply(y,2,sum)
temp=prior_parse(prior,k)
attach(temp)

n=length(y)
if(sum(x[,1])!=n){
	tst=apply(x, 2, sum)
	ind=which(tst==n)
	if(length(ind)>0){
		print("sar_g: intercept term must be in first column of the x-matrix")
	}
	if(length(ind)==0){
		cflag=0
		p=ncol(x)
		}
	}
	#if(sum(x[,1]==n)){
	else {
	cflag=1
	p=ncol(x)-1
	}	

	results$cflag=cflag
	results$p=p

out_temp = set_eigs(eflag,W,rmin,rmax,n)
rmin=out_temp$rmin
rmax=out_temp$rmax
detval = set_lndet(ldetflag, W, rmin, rmax, detval, order, iter)

iiter=50
o=100

diag_ests=matrix(rep(0,n),n,0)

for(iii in 1:iiter){
	u=matrix(rnorm(n),n,1)
	umat=u[,rep(1,o)]
	wumat=matrix(rep(0,n*o),n,o)
	wu=u
	wumat[,1]=wu
	for(ii in 2:o){
		wu=W%*%wu
		wumat[,ii]=wu
		}
	diag_estimates_iii=umat*wumat
	diag_ests=diag_estimates_iii	
	}
estimated_diags=diag_ests/iiter

bsave = matrix(rep(0,(ndraw-nomit)*k),ndraw-nomit,k)
psave = matrix(rep(0,(ndraw-nomit)),ndraw-nomit,1)	
ymean= matrix(rep(0,n),n,1)
acc_rate= matrix(rep(0,ndraw),ndraw,1)

total=matrix(rep(0,(ndraw-nomit)*p),ndraw-nomit,p)
total_obs=matrix(rep(0,n*p),n,p)
direct=matrix(rep(0,(ndraw-nomit)*p),ndraw-nomit,p)
indirect=matrix(rep(0,(ndraw-nomit)*p),ndraw-nomit,p)

avg_total=matrix(rep(0,p),p,1)
avg_direct=matrix(rep(0,p),p,1)]
avg_indirect=matrix(rep(0,p),p,1)

TI = solve(T);
TIc = TI%*%c_beta;

In = matrix(rep(1,n),n,1);

Wy=W%*%y
Wadd=W+t(W)
WtW=t(W)%*%W
sige=1

iter=1
xpx=t(x)%*%x
xpy=t(x)%*%y
Wy=W%*%y
xpWy=t(x)%*%Wy

indl=which(yin==0)
nobs0=length(ind1)
ind2=which(yin==1)
nobs1=length(ind2)

while(iter <= ndraw){
	AI=solve(xpx+as.numeric(sige)*TI)
	ys=y-rho*Wy
	b = t(xs)%*%yss + as.numeric(sige)*TIc;
	b0 = solve((t(xs)%*%xs + as.numeric(sige)*TI),b);
	bhat = norm_rnd(as.numeric(sige)*AI) + b0;  
	xb = xs%*%bhat; 
	
	b0 = solve((t(xs)%*%xs + as.numeric(sige)*TI ),(t(xs)%*%ys + as.numeric(sige)*TIc));
	bd = solve((t(xs)%*%xs + as.numeric(sige)*TI),(t(xs)%*%Wys + as.numeric(sige)*TIc));
	e0 = ys - xs%*%b0;
	ed = Wys - xs%*%bd;
	epe0 = t(e0)%*%e0;
	eped = t(ed)%*%ed;
	epe0d = t(ed)%*%e0;
	rho = draw_rho(detval,epe0,eped,epe0d,n,k,rho,a1,a2)
	
	hh= diag(n)-rho%*%W
	mu=solve(hh,xb)
	tauinv=diag(n)-rho*Wadd+rho*rho*WtW
	
	aa=diag(tauinv) ##??
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
			if(yin[i,1]==0)	t1=normrt_rnd(0,1,muuse)
			if(yin[i,1]==1)	t1=normrt_rnd(0,1,muuse)
			z[i,1]=aa+h[i,1]*t1	
			}
		}
	y=mu+z
	Wy=W%*%y
	
	if(iter>nomit){
		bsave[iter-nomit,1:k]=t(bhat)
		psave[iter-nomit,1]=rho
		ymean=ymean+y
		
		rhovec=t(rho^c(0:(o-1)))
		if(cflag==1)	beff=bhat[2:length(bhat),1] #?
		if(cflag==0)	beff=bhat
		s=solve(hh,diag(n))
		pdfz=stdn_pdf(mu)
		for(kk=1:p){
			avg_direct[kk,1]=t(pdfz)%*%estimated_diags%*%rhovec*beff[kk,1]/n
			#dd=spdiags()	##spdiags ??
			avg_total[kk,1]=mean(apply(dd%*%s*beff[kk,1],2,sum))
			total_obs[,kk]=total_obs[,kk]+apply(dd%*%s*beff[kk,1],2,sum)
			avg_indirect[kk,1]=avg_total[kk,1]-avg_direct[kk,1]
			}
		total[iter-nomit,]=t(avg_total)
		direct[iter-nomit,]t(avg_direct)
		indirect[iter-nomit,]=t(avg_indirect)	
		}
	
	iter=iter+1
	}
beta = apply(bsave,2,mean)
rho=mean(psave)
ymean=ymean/(ndraw-nomit)
results$sige=1
sige=1

nobs=dim(x)[1]
nvar=dim(x)[2]
Wy=W%*%ymean
AI=solve(t(x)%*%x+as.numeric(sige)*TI)
b0=AI%*%(t(x)%*%ymean+as.numeric(sige)*TIc)
bd=AI%*%(t(x)%*%Wy+as.numeric(sige)*TIc)
e0=ymean-x%*%b0
ed=Wy-x%*%bd
epe0=t(e0)%*%e0
eped=t(ed)%*%ed
epe0d=t(ed)%*%e0

logdetx=log(det(t(x)%*%x)+as.numeric(sige)*TI)

if(inform_flag==0)	mlike=sar_marginal(detval,e0,ed,epe0,eped,epe0d,nobs,nvar,logdetx,a1,a2)
if(inform_flag==1)	mlike=sar_marginal2(detval,e0,ed,epe0,eped,epe0d,nobs,nvar,logdetx,a1,a2,c,TI,x,y,sige,W)

results$meth  = 'sarp_g';
results$ymean = ymean;
results$total = total;
results$direct = direct;
results$indirect = indirect;
results$total_obs = total_obs;
results$beta = beta;
results$rho = rho;
results$bdraw = bsave;
results$pdraw = psave;
results$bmean = cc;
results$bstd  = sqrt(diag(T));
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



























