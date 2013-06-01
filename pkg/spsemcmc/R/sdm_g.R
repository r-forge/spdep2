#sdm_g.R
# list of external functions : draw_rho,sdm_parse,sdm_eigs,sdm_lndet,
#sdm_marginal,sdm_marginal2
#

sdm_g <- function(y,x,W,ndraw,nomit,prior){
##############################################
results <- list()
n = nrow(y)

#checking if the user handled the intercept term okay
if(sum(x[,1])!=n){
	tst=apply(x, 2, sum)
	ind=which(tst==n)
	if(length(ind)>0){
		stop("sdm_g: intercept term must be in first column of the x-matrix")
	}
	if(length(ind)==0){
		xsdm=cbind(x,W%*%x)
		cflag=0
		p=ncol(x)
		}
	}
if(sum(x[,1])==n){
	xsdm=cbind(x,W%*%x[,-1])
	cflag=1
	p=ncol(x)-1
}	

nobs=nrow(xsdm)
k=ncol(xsdm)

results$nobs  = n;
results$nvar  = k;
results$y = y; 
# ommitted: check no of parameters input

pprior=prior_parse(prior,k) 
attach(pprior)


n1=nrow(W)
n2=ncol(W)

results$cflag = cflag;
results$p = p;

if(n1!=n2)
{
        stop('sar_g: wrong size weight matrix W');
}
else
{
        if(n1!=n)
                stop('sar_g: wrong size weight matrix W');
}

#    [nchk junk] = size(y);
nchk <- nrow(y)
if(nchk!=n)
        stop('sar_g: wrong size y vector input');

results$order = order;
results$iter = iter;

checkk<-nrow(c_beta)
junk<-ncol(c_beta)

if((checkk!=k) | (junk!=1))
        stop("sdm_g: prior means are wrong")


checkk<-nrow(Tbeta)
junk<-ncol(Tbeta)

if((checkk!=k) | (junk!=k))
        stop("sdm_g: prior bcov is wrong")


out_temp = set_eigs(eflag,W,rmin,rmax,n);
rmin=out_temp$rmin
rmax=out_temp$rmax

detval = set_lndet(ldetflag,W,rmin,rmax,detval,order,iter);

bsave = matrix(0,ndraw-nomit,k);

if (mm!= 0)
	rsave = matrix(0,ndraw-nomit,1)

psave = matrix(0,ndraw-nomit,1)
ssave = matrix(0,ndraw-nomit,1)
vmean = matrix(0,n,1)
acc_rate = matrix(0,ndraw,1)

ntrs = 101;

#comment: initialisations

TI = solve(Tbeta);
TIc = TI%*%c_beta; ##name changed from c to c_beta for less ambiguity
iter = 1;

In = matrix(1,n,1);  ##name changed in to In
V = In;
vi = In;
Wy = W%*%y;
x = xsdm;

#commment: homosc. and heteroc. cases are handled in a single loop

iter = 1;
	while (iter <= ndraw){ #% start sampling;

          #% update beta   
          xs = matmul(sqrt(V),x); ## code for matmul in support functions
          ys = sqrt(V)*y;
          Wys = sqrt(V)*Wy
          AI = solve(t(xs)%*%xs + sige*TI)#,diag(k))         
          yss = ys - rho*Wys;          
          xpy = t(xs)%*%yss;
          b = t(xs)%*%yss + sige*TIc;
          b0 = solve(t(xs)%*%xs + sige*TI,b);
          bhat = mvrnorm(1, b0, sige*AI)
          xb = xs%*%bhat; 
                    
          #% update sige
          nu1 = n + 2*nu; 
          e = (yss - xb);
          d1 = 2*d0 + t(e)%*%e;
          chi = rchisq(1,nu1); ##code for chis_rnd
          sige = as.numeric(d1/chi); 
	  ## not present in homosc case.
	  if(novi_flag==0){
	  	#% update vi
          	ev = y - rho*Wy - x%*%bhat; 
          	#chiv = chis_rnd(n,rval+1);  
          	chiv = matrix(rchisq(n,rval+1),n,1); 
          	vi = ((ev*ev/sige) + In%*%rval)/chiv; 
          	V = In/vi; 
                        
          	#% update rval
          	if (mm != 0)
			rval = rgamma(1,shape=mm, rate=kk);  
          }
          
      #% we use griddy Gibbs to perform rho-draw
          b0 = solve((t(xs)%*%xs + sige*TI ),(t(xs)%*%ys + sige*TIc));
          bd = solve((t(xs)%*%xs + sige*TI),(t(xs)%*%Wys + sige*TIc));
          e0 = ys - xs%*%b0;
          ed = Wys - xs%*%bd;
          epe0 = t(e0)%*%e0;
          eped = t(ed)%*%ed;
          epe0d = t(ed)%*%e0;
          rho = draw_rho(detval,epe0,eped,epe0d,n,k,rho,a1,a2);  #function later


#	print(t(bhat))
#	print(sige)
#	print(rho)
          
    if (iter > nomit){ ## if we are past burn-in, save the draws
        bsave[iter-nomit,1:k] = t(bhat);
        ssave[iter-nomit,1] = sige;
        psave[iter-nomit,1] = rho;
        vmean = vmean + vi; 
	if(novi_flag==0 & mm != 0)   
			rsave[iter-nomit,1] = rval
	
   	         
    }
                    
iter = iter + 1; 
    
}
uiter=50;
maxorderu=100;
nobs = n;
rv=matrix(rnorm(nobs*uiter),nobs,uiter);
tracew=matrix(0,maxorderu,1);
wjjju=rv;
for (jjj in 1:maxorderu){
    wjjju=W%*%wjjju;
    tracew[jjj]=mean(apply(rv*wjjju, 2, mean));
}

traces=tracew;
traces[1,1]=0;
traces[2,1]=sum(apply(t(W)*W, 2, sum))/nobs;
trs=matrix(c(1,traces), ncol=1);
ntrs=length(trs);
trbig=t(trs);
trbig2=matrix(c(trbig[1,2:ncol(trbig)],trbig[1,ncol(trbig)]))
trmat=rbind(trbig,t(trbig2))

if (cflag == 1){
		bdraws = as.matrix(bsave[,2:ncol(bsave)]);
}
if(cflag == 0){
		bdraws = bsave;
}
pdraws = psave;
ree = 0:(ntrs-1);
rmat = matrix(0,1,ntrs);
total = array(0,c(ndraw-nomit,p,ntrs)) ;
direct = array(0,c(ndraw-nomit,p,ntrs));	## 3D matrix in R
indirect = array(0,c(ndraw-nomit,p,ntrs));

for (i in 1:(ndraw-nomit)){
    rmat = pdraws[i,1]^ree;
    for (j in 1:p){
            bbeta = c(bdraws[i,j], bdraws[i, j+p])
            total[i,j,] = sum(bbeta)*rmat; ##used beta instead of beta[1,1]
    direct[i,j,] = (bbeta%*%trmat)*rmat; 
    indirect[i,j,] = total[i,j,] - direct[i,j,];
    }

}

bmean = apply(bsave, 2, mean);
bbeta = t(bmean);
rho = apply(psave, 2, mean);
sige = apply(ssave, 2, mean);
results$sige = sige;
vmean = vmean/(ndraw-nomit);
V = In/vmean;

### calculate log marginal likelihood

nobs=nrow(x)
nvar=ncol(x)
xs = matmul(x,sqrt(V));
ys = sqrt(V)*y;
Wys = sqrt(V)*Wy
AI = solve(t(xs)%*%xs + sige*TI);
b0 = AI%*%(t(xs)%*%ys + sige*TIc);
bd = AI%*%(t(xs)%*%Wys + sige*TIc);
e0 = ys - xs%*%b0;
ed = Wys - xs%*%bd;
epe0 = t(e0)%*%e0;
eped = t(ed)%*%ed;
epe0d = t(ed)%*%e0;
logdetx = log(det(t(xs)%*%xs + sige*TI));
  if (inform_flag == 0){
		mlike = rho_marginal(detval,e0,ed,epe0,eped,epe0d,nobs,nvar,logdetx,a1,a2);}
  if(inform_flag == 1){
#   mlike = sar_marginal2(detval,e0,ed,epe0,eped,epe0d,nobs,nvar,a1,a2,c_beta,TI,xs,ys,sige,W);
		mlike = rho_marginal2(detval,e0,ed,epe0,eped,epe0d,nobs,nvar,logdetx,a1,a2,c_beta,TI,xs,ys,sige);#VIRGILIO: W is not needed
  }
  
yhat = solve((diag(nobs) - rho*W),(xs%*%bmean));
e = y - yhat; 

## compute R squared

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

results$meth  = 'sar_g';
results$total = total;
results$direct = direct;
results$indirect = indirect;
results$beta = bbeta;
results$rho = rho;
results$bdraw = bsave;
results$pdraw = psave;
results$sdraw = ssave;
results$beta_std = t(apply(bsave,2,sd));
results$sige_std = apply(ssave,2,sd);
results$rho_std = apply(psave,2,sd);
results$mlike = mlike;
results$vmean = vmean;
results$yhat  = yhat;
results$resid = e;
results$bmean = c_beta;
results$bstd  = sqrt(diag(Tbeta));
results$ndraw = ndraw;
results$nomit = nomit;
#results$time  = time;
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
results$cflag=cflag;
results$p=p;

if (mm!= 0){
	results$rdraw = rsave;
	results$m     = mm;
	results$k     = kk;
}
else
{
	results$r     = rval;
	results$rdraw = 0;
}

detach(pprior)
return(results)
}













