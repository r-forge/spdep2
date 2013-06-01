#sac_g

sac_g <- function(y,x,W1,W2,ndraw,nomit,prior){
results <-list()
n=dim(y)[1]
junk = dim(y)[2]
results$y = y;
n1=dim(x)[1]
k=dim(x)[2]
n3=dim(W1)[1]
n4=dim(W1)[2]
n5=dim(W2)[1]
n6=dim(W2)[2]

n=length(y)
if(sum(x[,1])!=n)
	{
	tst=apply(x, 2, sum)
	ind=which(tst==n)
	if(length(ind)>0)
		{
		print("sar_g: intercept term must be in first column of the x-matrix")
		}
	if(length(ind)==0)
		{
		cflag=0
		p=ncol(x)
		}
	}
if(sum(x[,1]==n))
	{
	cflag=1
	p=ncol(x)-1
	}	


results$cflag = cflag;
results$p = p;

temp=prior_parse(prior,k)
attach(temp)

results$order=order
results$iter=iter

P = matrix(rep(1,n),n,1); 
In = matrix(rep(1,n),n,1); ##; % initial value for V   
ys = y;
          
vmean = matrix(rep(0,n),n,1); 
yhat = matrix(rep(0,n),n,1);

out_temp = set_eigs(eflag,W,rmin,rmax,n)
rmin=out_temp$rmin
rmax=out_temp$rmax
out_temp = set_eigs(eflag,W,lmin,lmax,n)
lmin=out_temp$rmin
lmax=out_temp$rmax

results$rmin = rmin;
results$rmax = rmax;
results$lmin = lmin;
results$lmax = lmax;

results$lflag = ldetflag;

detval1 = set_lndet(ldetflag, W1, rmin, rmax, detval1, order, iter)
detval2 = set_lndet(ldetflag, W2, rmin, rmax, detval2, order, iter)

bsave = matrix(rep(0,(ndraw-nomit)*k),ndraw-nomit,k)
psave = matrix(rep(0,(ndraw-nomit)),ndraw-nomit,1)
lsave = matrix(rep(0,(ndraw-nomit)),ndraw-nomit,1)
ssave = matrix(rep(0,(ndraw-nomit)),ndraw-nomit,1)
vmean = matrix(rep(0,n),n,1)
acc_rate1 = matrix(rep(0,ndraw),ndraw,1)
acc_rate2 = matrix(rep(0,ndraw),ndraw,1)

TI = solve(T);
TIc = TI%*%c_beta;
iter = 1;
In = matrix(rep(1,n),n,1)
V = In;
Wy = W1%*%y;  
Wx = W1%*%x;
vi = In;
P = In;
V = vi;
#In = diag(n);   

B = diag(n) - lambda*W2;
A = diag(n) - rho*W1;

nu1 = n + 2*nu; 

##sampling

iter = 1;
acc1 = 0;
acc2 = 0;

while (iter <= ndraw)
	{

	#   % update beta 
	xtil = B%*%x;                            ##% xtil is used with vi, so not a waste creating here
	ytil = B%*%A%*%y;                          ##% ytil is used with vi, so not a waste creating here
	xstar = matmul(P,xtil);                ##% P here is the sqrt of inv of covar matrix, nx1
	ystar = P*ytil;                       ##% P is all ones when assuming homoscedasticity
	Hinv = solve(t(xstar)%*%xstar + as.numeric(sige)*TI);    ##% Hinv is covariance of beta
	b0 = Hinv%*%(t(xstar)%*%ystar + as.numeric(sige)*TIc);   ##% b0 is mean of beta
	bhat = norm_rnd(Hinv) + b0;            ##% bhat is simulated beta; norm_rnd is MVN, mean 0, covar Hinv  ##code for norm_rnd
	xb = x%*%bhat;


	##   % update sige (here we take vi into account)
	Bx = (diag(n) - lambda*W2)%*%x; ##use sparse
	b = solve((t(Bx)%*%Bx),(t(Bx)%*%B%*%A%*%y));
	e = B%*%(A%*%y - x%*%b);
	ev = P*e;
	d1 = 2*d0 + t(ev)%*%ev;
	chi = chis_rnd(1,nu1); 	##code for chis_rnd
	sige = d1/chi;

	##   % update vi (based on e, without taking vi into account)
	if(novi_flag==0)
		{
		chiv = chis_rnd(n,rval+1);
		vi = (((e*e)/as.numeric(sige)) + In*rval)/chiv;
		P = In/vi; 
		}


	##   % update lambda using metropolis-hastings
	##          % numerical integration is too slow here
	xb = x%*%bhat;
	rhox = c_lambda_sac(rho,lambda,y,x,bhat,sige,W1,W2,detval2,P,a1,a2);    	##c_lambda
	accept = 0;
	lambda2 = lambda + cc1*rnorm(1);		##cc1 matrix??
	while (accept == 0)
		{
		if ((lambda2 > lmin) & (lambda2 < lmax)) 
			{
			accept = 1; 
			} 
		else 
			{
			lambda2 = lambda + cc1*rnorm(1);
			}
		} 
	rhoy = c_lambda_sac(rho,lambda2,y,x,bhat,sige,W1,W2,detval2,P,a1,a2);
	ru = runif(1);			##unif_rnd
	if ((rhoy - rhox) > exp(1))
		{
		pp = 1;
		}
	else 
		{         
		ratio = exp(rhoy-rhox);
		pp = min(ratio); 				##ratio??
		}
	if (ru < pp)
		{
		lambda = lambda2;
		acc1 = acc1 + 1;
		}
	acc_rate1[iter,1] = acc1/iter;		##acc_rate1 ??
	###% update cc based on std of rho draws
	if (acc_rate1[iter,1] < 0.4)
		{
		cc1 = cc1/1.1;
		}
	if (acc_rate1[iter,1] > 0.6)
		{
		cc1 = cc1*1.1;
		}

	B = diag(n) - lambda*W2;


	##% update rho using metropolis-hastings
	##   % numerical integration is too slow here
	xb = x%*%bhat;
	rhox = c_rho_sac(rho,lambda,y,x,bhat,sige,W1,W2,detval1,P,a1,a2);
	accept = 0;
	rho2 = rho + cc2*rnorm(1);
	while (accept == 0)
		{
		if ((rho2 > lmin) & (rho2 < lmax))
			{
			accept = 1;  
			}
		else
			{
			rho2 = rho + cc2*rnorm(1);
			}
		}	
	rhoy = c_rho_sac(rho2,lambda,y,x,bhat,sige,W1,W2,detval1,P,a1,a2);
	ru = runif(1);
	if ((rhoy - rhox) > exp(1))
		{
		 pp = 1;
		}  
	else
		{  
		ratio = exp(rhoy-rhox);
		pp = min(1,ratio);
		}
	if (ru < pp)
		{
		rho = rho2;
		acc2 = acc2 + 1;
		}
	acc_rate2[iter,1] = acc2/iter;			##acc_rate2
	##% update cc based on std of rho draws
	if (acc_rate2[iter,1] < 0.4)       cc2 = cc2/1.1;

	if (acc_rate2[iter,1] > 0.6)       cc2 = cc2*1.1;


	A = diag(n) - rho*W1;


	if (iter > nomit)
		{ ##% if we are past burn-in, save the draws
		bsave[iter-nomit,1:k] = t(bhat);
		ssave[iter-nomit,1] = sige;
		psave[iter-nomit,1] = rho;
		lsave[iter-nomit,1] = lambda;
		vmean = vmean + vi;        
		}


	iter = iter + 1; 
	##waitbar(iter/ndraw);         
	}
	
uiter=50;
maxorderu=100;
nobs = n;
rv=matrix(rnorm(nobs*uiter),nobs,uiter);
tracew=matrix(rep(0,maxorderu),maxorderu,1);
wjjju=rv;
for (jjj in 1:maxorderu)
	{
    wjjju=W1%*%wjjju;
    tracew[jjj]=mean(apply(rv*wjjju,2,mean));
    
	}

traces=tracew;
traces[1]=0;
traces[2]=sum(sum(t(W1)*W1))/nobs;#####mark###
trs=rbind(1,traces);
ntrs=length(trs);
trbig=t(trs);
                 
        if (cflag == 1){
        bdraws = bsave[,2:length(bsave)];}
        if (cflag == 0){
        bdraws = bsave;
        } 
        pdraws = psave;

        ree = 0:(ntrs-1);
		
        rmat = matrix(rep(0,ntrs),1,ntrs);####three dimentional matrix in R??
        total = array(0,c(ndraw-nomit,p,ntrs)) ;
        direct = array(0,c(ndraw-nomit,p,ntrs));	## 3D matrix in R
        indirect = array(0,c(ndraw-nomit,p,ntrs));
        
for (i in 1:ndraw-nomit){
    rmat = pdraws[i,1]^ree;
    for (j in 1:p){
            beta = bdraws[i,j];
            total[i,j,] = beta*rmat;
    direct[i,j,] = (beta%*%trbig)*rmat;
    indirect[i,j,] = total[i,j,] - direct[i,j,];
    }

}

#% compute posterior means 
vmean = vmean/(ndraw-nomit);
bmean = apply(bsave,2,mean);
bmean = t(bmean);
rho = apply(psave,2,mean);
lam = apply(lsave,2,mean);
sige = apply(ssave,2,mean);


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
results$bmean = c_beta;
results$bstd  = sqrt(diag(T));
results$rsqr  = rsqr;
results$rbar = 1 - (rsqr1/rsqr2); #% rbar-squared
results$sige = sige;
results$nobs  = n;
results$nvar  = nvar;
results$ndraw = ndraw;
results$nomit = nomit;
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
































