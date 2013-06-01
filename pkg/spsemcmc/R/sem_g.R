#sem_g.R

sem_g <-function(y,x,W,ndraw,nomit,prior){

results <- list()

n = nrow(y)
#junk=ncol(y)
results$y=y

n1=nrow(x)
k=ncol(x)
n3=nrow(W)
n4=ncol(W)

if(n1!=n)
	stop("sem_g: x-matrix constains wrong number of observations")
if(n3!=n4)
	stop("sem_g: W matrix is not square")
if(n3!=n)
	stop("sem_g: W matrix is not the same size at y,x")

pprior=prior_parse(prior,k=k)
attach(pprior)

results$order=order
results$iter=iter

checkk<-nrow(c_beta)
junk<-ncol(c_beta)

if((checkk!=k) | (junk!=1))
	stop("sem_g: prior means are wrong")

#V=matrix(1, n,1)
#In=matrix(1, n,1)
#vi=In

bsave=matrix(0, ndraw-nomit, k) ## changed from 1 to k
ssave=matrix(0, ndraw-nomit, 1)
psave=matrix(0, ndraw-nomit, 1)
vmean=matrix(1, n, 1)
yhat=matrix(1, n, 1)

acc_rate=matrix(0, ndraw,1)

if(mm!=0)	
	rsave=matrix(0, ndraw-nomit,1)

temp = set_eigs(eflag,W,rmin,rmax,n);
rmin=temp$rmin
rmax=temp$rmax

detval = set_lndet(ldetflag,W,rmin,rmax,detval,order,iter);

## Initialisation

TI = solve(Tbeta);
TIc = TI%*%c_beta; ##name changed from c to c_beta for less ambiguity
iter = 1;
In = matrix(1, n,1);  ##name changed in to In
V = In;
ys=y*sqrt(V)
Wy=W%*%y
Wx=W%*%x
vi=In
V=vi

#commment: homosc. and heteroc. cases are handled in a single loop

iter = 1;
acc=0;


while (iter <= ndraw){ #% start sampling;
			  
	#% update beta   
	xs = matmul(sqrt(V),x);
	ys = sqrt(V)*y;		
 	Wxs = W%*%xs;
        Wys = W%*%ys;
        xss = xs - rho*Wxs;
        AI = solve(t(xss)%*%xss+sige*TI)
        yss = ys - rho*Wys;
	b = t(xss)%*%yss + sige*TIc;
	b0 = AI%*%b;
	bhat = mvrnorm(1, b0, sige*AI); 
	
	#update sige:
	
	nu1 = n + 2*nu; 
	e = yss-xss%*%bhat;
	ed = e - rho*W%*%e;
	d1 = 2*d0 + t(ed)%*%ed;
	chi = rchisq(1,nu1);
	sige = as.numeric(d1/chi);
	
	#update vi when novi_flag==0
	if(novi_flag==0){
		#ev = y - rho*Wy - x%*%bhat; 
		ev = ys - xs%*%bhat; 
		#chiv = chis_rnd(n,rval+1);  
		chiv = matrix(rchisq(n,rval+1),n,1); 
		vi = ((ev*ev/sige) + In*rval)/chiv; 
		V = In/vi; 
					
		#% update rval
		if (mm != 0)   
			rval = rgamma(1,shape=mm,rate=kk);  
	  }
	if(metflag==0)
	{
		  rho = draw_rho_sem(detval,y,x,Wy,Wx,V,n,k,rmin,rmax,rho);
	}
	else 
	{
		xb = x%*%bhat;
		rhox = c_rho_sem(rho,y,x,bhat,sige,W,detval,V,a1,a2);
		accept = 0;
		rho2 = rho + cc*rnorm(1);
		while(accept==0)
		{
			if((rho2>rmin)&(rho2<rmax))	
				accept=1
			else 
				rho2=rho+cc*rnorm(1)
		}
		rhoy = c_rho_sem(rho2,y,x,bhat,sige,W,detval,V,a1,a2);
		ru = runif(1);
#		if(exp(rhoy-rhox)>1)	
#			p=1
#		else
#		{
			ratio=exp(rhoy-rhox)
			p=min(1,ratio)
#		}
		if(ru<p)
		{
			rho=rho2
			acc=acc+1
		}
		acc_rate[iter,1] = acc/iter
		
		#update cc based on sd of rho draws
		
		if(acc_rate[iter,1]<0.4)	cc=cc/1.1
		if(acc_rate[iter,1]>0.6)	cc=cc*1.1	
	}	
	
	if(iter>nomit)
	{
		bsave[iter-nomit,1:k] = t(bhat);
		ssave[iter-nomit,1] = sige;
		psave[iter-nomit,1] = rho;
		vmean = vmean + vi;
		if(mm!=0)	
			rsave[iter-nomit,1]=rval
	}
	iter=iter+1
}

vmean = vmean/(ndraw-nomit);
bmean = apply(bsave,2,mean);
#bmean = t(bmean);
rho = apply(psave,2,mean);
V = In/vmean;
ys = y*sqrt(V);
xs = matmul(x,sqrt(V));
Wys = W%*%ys;
Wxs = W%*%xs;	
nobs=nrow(xs)
nvar=nrow(xs)	

if (mlog == 1)
	{
    if (inform_flag == 0)
		{
			mlike = sem_marginal(detval,ys,xs,Wys,Wxs,nobs,nvar,a1,a2);
		}
    else
    {
    mlike = sem_marginal2(detval,ys,xs,Wys,Wxs,nobs,nvar,a1,a2,c,TI,sige);
    }
}

n=nrow(x)
nvar=ncol(x)
yhat=x%*%bmean
y=results$y
n=length(y)
e = y-yhat;
eD = e - rho*W%*%e;
epe = t(eD)%*%eD;

sigu = epe;
sige = sigu/(n-nvar);
ym = y - mean(y);
rsqr1 = sigu;
rsqr2 = t(ym)%*%ym;
rsqr = 1.0 - rsqr1/rsqr2; ##% conventional r-squared
rsqr1 = rsqr1/(n-nvar);
rsqr2 = rsqr2/(n-1.0);

results$meth  = 'sem_g';
results$beta = bmean;
results$rho = rho;
results$sige = sige;
results$bdraw = bsave;
results$pdraw = psave;
results$sdraw = ssave;
results$vmean = vmean;
results$yhat  = yhat;
results$bmean = c_beta;
results$bstd  = sqrt(diag(Tbeta));
results$rsqr  = rsqr;
results$rbar = 1 - (rsqr1/rsqr2); #% rbar-squared
results$sige = sige;
results$nobs  = n;
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
if(mm!=0)
	{
	results$rdraw=rsave
	results$m=mm
	results$k=kk
	}
else
	{
	results$r=rval
	results$rdraw=0
	}

detach(pprior)
return(results)
}	

