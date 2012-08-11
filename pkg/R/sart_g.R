#sart_g.R

sart_g <-function(y,x,W,ndraw,nomit,prior){

n=dim(y)[1]
junk=dim(y)[2]
n1=dim(x)[1]
k=dim(x)[2]
n3=dim(W)[1]
n4=dim(W)[2]
yin=y

results <-list()

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
	if(sum(x[,1]==n)){
	#else {
	cflag=1
	p=ncol(x)-1
	}	
	
results$cflag=cflag
results$p=p

temp=prior_parse(prior,k)
attach(temp)

ind=which(y<=vflag)
nobsc=length(ind)

checkk=dim(c_beta)[1]

results$y=y
results$nobs=n
results$nvar=k
results$order=order
results$iter=iter

out_temp = set_eigs(eflag,W,rmin,rmax,n)
rmin=out_temp$rmin
rmax=out_temp$rmax
detval = set_lndet(ldetflag, W, rmin, rmax, detval, order, iter)

bsave = matrix(0,ndraw-nomit,k)
psave = matrix(0,ndraw-nomit,1)	
ssave = matrix(0,ndraw-nomit,1)
ymean = matrix(0,n,1)
acc_rate= matrix(0,ndraw,1)

TI=solve(T)
TIc=TI%*%c_beta
acc=0

Wy=W%*%y
iter=1
xpx=t(x)%*%x

ind1=which(yin==0)
nobs0=length(ind1)
ind2=which(yin >0)
nobs1=length(ind2)



while(iter <= ndraw){
	AI=solve(xpx+as.numeric(sige)*TI)
	ys=y-rho*Wy
	b=t(x)%*%ys+as.numeric(sige)*TIc
	b0=AI%*%b
	bhat=norm_rnd(as.numeric(sige)*AI)+b0
	xb=x%*%bhat
	
	nu1=n+2*nu
	e=ys-xb
	d1=2*d0+t(e)%*%e
	chi=chis_rnd(1,nu1)
	sige=d1/chi
	
	if(metflag==1){
		rhox=c_sar(rho,y,xb,sige,W,detval)
		accept=0
		rho2=rho+cc*rnorm(1)
		while(accept==0){
			if((rho2>rmin)&(rho2<rmax))		accept=1
			else 	rho2=rho+cc*rnorm(1)
			}
		rhoy=c_sar(rho2,y,xb,sige,W,detval)
		ru=runif(1)
		if((rhoy-rhox)>exp(1))	pp=1
		else { 	ratio=exp(rhoy-rhox)
			pp=min(1,ratio)
			}
		if(ru < pp){
			rho=rho2
			acc=acc+1
			}
		acc_rate[iter,1]=acc/iter
		
		if(acc_rate[iter,1]<0.4)	cc=cc/1.1
		if(acc_rate[iter,1]>0.6)	cc=cc*1.1	
		}
		
		if(metflag==0){
		b0=solve((t(x)%*%x),(t(x)%*%y))
		bd=solve((t(x)%*%x),(t(x)%*%Wy))
		e0=y-x%*%b0
		ed=Wy-x%*%bd
		epe0=t(e0)%*%e0
		eped=t(ed)%*%ed
		epe0d=t(ed)%*%e0
		rho=draw_rho(detval,epe0,eped,epe0d,n,k,rho)
		}
		
	if(iter==1)	z=matrix(0,n,1)	
	z=matrix(0,n,1)
	
	h=diag(n)-rho*W
	mu=solve(h,xb)
	tauinv=(t(h)%*%h)/as.numeric(sige)
	aa=diag(tauinv)
	h=matrix(1,n,1)/sqrt(aa)
	c_tmp=matdiv(-tauinv,as.matrix(aa))
	ctilde=c_tmp-diag(diag(c_tmp))
	
	for(initer in 1:nsample){
		for(i in 1:n){
			if(yin[i,1]==0){
				aa=ctilde[i,]%*%z
				muuse=(-mu[i,1]-aa)/h[i,1]
				t1=normrt_rnd(0,1,muuse)
				z[i,1]=aa+h[i,1]*t1
				}
			}
		}
	y[ind1,1]=mu[ind1,1]+z[ind1,1]
	Wy=W%*%y
	
	if(iter>nomit){
		bsave[iter-nomit,1:k]=t(bhat)
		psave[iter-nomit,1]=rho
		ssave[iter-nomit,1]=sige
		ymean=ymean+y
		}
	iter=iter+1
	}
uiter=50
maxorderu=100
nobs=n
rv=matrix(rnorm(nobs*uiter),nobs,uiter)
tracew=matrix(0,maxorderu,1)
wjjju=rv
for(jjj in 1:maxorderu){
	wjjju=W%*%wjjju
	tracew[jjj]=mean(apply(rv*wjjju,2,mean))
	}
traces=tracew
traces[1]=0
traces[2]=sum(apply((t(W)*W),2,sum))/nobs
trs=rbind(1,traces)
ntrs=length(trs)
trbig=t(trs)

if(cflag==1){
	bdraws=bsave[,2:dim(bsave)[2]]
	}	
if(cflag==0)	bdraws=bsave

pdraws=psave
ree=seq(0,(ntrs-1),1)
rmat=matrix(0,1,ntrs)
total=array(0,c(ndraw-nomit,p,ntrs))
direct=array(0,c(ndraw-nomit,p,ntrs))
indirect=array(0,c(ndraw-nomit,p,ntrs))

for(i in 1:ndraw-nomit){
	rmat=pdraws[i,1]^ree
	for(j in 1:p){
		beta=bdraws[i,j]
		total[i,j,]=beta*rmat
		direct[i,j,]=(beta%*%trbig)*rmat
		indirect[i,j,]=total[i,j,]-direct[i,j,]
		}
	}	

ymean=ymean/(ndraw-nomit)
bsave=bsave[(nomit+1):dim(bsave)[1],]
psave=t(psave[(nomit+1):length(psave),1])
ssave=t(ssave[(nomit+1):length(ssave),1])

bmean=apply(bsave,2,mean)#?
beta=t(bmean)
rho=apply(psave,2,mean)
sige=apply(ssave,2,mean)

#yhat=solve((diag(n)-rho*W),(x%*%(beta))) ##dimention mismatch ?

results$meth  = 'sart_g';
results$beta = beta;
results$rho = rho;
results$bdraw = bsave;
results$pdraw = psave;
results$sdraw = ssave;
#results$yhat  = yhat;
results$ymean = ymean;
results$nsteps = nsample;
results$total = total;
results$direct = direct;
results$indirect = indirect;
results$bmean = c;
results$bstd  = sqrt(diag(T));
results$nobs  = n;
results$nvar  = k;
results$ndraw = ndraw;
results$nomit = nomit;
results$tflag = 'plevel';
results$acc = acc_rate;
results$order = order;
results$rmax = rmax; 
results$rmin = rmin;
results$lflag = ldetflag;
results$lndet = detval;
#results$priorb = prior_beta;
results$limit = vflag;
results$trunc = cflag;
results$nobsc = nobsc;

return(results)

}




























