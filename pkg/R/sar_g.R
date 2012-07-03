## just to test sar_g with out function encapsulations
sar_g <- function(y,x,W,ndraw,nomit,prior){

n=dim(y)[1]
junk = dim(y)[2]
n1=dim(x)[1]
k=dim(x)[2]
n2=dim(W)[1]
n4=dim(W)[2]

results <-list()

results$nobs  = n;
results$nvar  = k;
results$y = y; 
temp=prior_parse(prior,k)
attach(temp)
##check for intercept
n=length(y)
if(sum(x[,1])!=n){
	tst=sum(x)
	ind=which(tst==n)
	if(length(ind)>0){
		print("sar_g: intercept term must be in first column of the x-matrix")
	}
	if(length(ind)==0){
	cflag=0
	p=ncol(x)
	}
	if(sum(x[,1]==n)){
	cflag=1
	p=ncol(x)-1
	}	
	results$cflag=cflag
	results$p=p
	}
	
out_temp = set_eigs(eflag,W,rmin,rmax,n)
rmin=out_temp$rmin
rmax=out_temp$rmax
detval = set_lndet(ldetflag, W, rmin, rmax, detval, order, iter)
results$order = order;
results$iter = iter;
#storage for draws
bsave = matrix(rep(0,(ndraw-nomit)*k),ndraw-nomit,k)
if (mm !=0) rsave = matrix(rep(0,ndraw-nomit),ndraw-nomit,1)
psave = matrix(rep(0,ndraw-nomit),ndraw-nomit,1)
ssave = matrix(rep(0,ndraw-nomit),ndraw-nomit,1)
vmean= matrix(rep(0,n),n,1)

#% ====== initializations
#% compute this stuff once to save time

TI = solve(T);
TIc = TI%*%c_beta;

In = matrix(rep(1,n),n,1);
V = In;
vi = In;
Wy = W%*%y;



#hwait = waitbar(0,'sar: MCMC sampling ...');

#t0 = clock;                  
iter = 1;
          while (iter <= ndraw){ #% start sampling;
                  
          #% update beta   
          xs = matmul(sqrt(V),x); ## code for matmul in support functions
          ys = sqrt(V)*y;
          Wys = sqrt(V)*Wy
          AI = solve((t(xs)%*%xs + as.numeric(sige)*TI),diag(k))         
          yss = ys - rho*Wys;          
          xpy = t(xs)%*%yss;
          b = t(xs)%*%yss + as.numeric(sige)*TIc;
          b0 = solve((t(xs)%*%xs + as.numeric(sige)*TI),b);
          bhat = norm_rnd(as.numeric(sige)*AI) + b0;  
          xb = xs%*%bhat; 
                    
          #% update sige
          nu1 = n + 2*nu; 
          e = (yss - xb);
          d1 = 2*d0 + t(e)%*%e;
          chi = chis_rnd(1,nu1); ##code for chis_rnd
          sige = d1/chi; 
	  ## not present in homosc case.
	  if(novi_flag==0){
	  	#% update vi
          	ev = y - rho*Wy - x%*%bhat; 
          	#chiv = chis_rnd(n,rval+1);  
          	chiv = matrix(rchisq(n,rval+1),n,1); 
          	vi = ((ev*ev/as.numeric(sige)) + In%*%rval)/chiv; 
          	V = In/vi; 
                        
          	#% update rval
          	if (mm != 0)   rval = rgamma(1,mm,1/kk);  
          }
          
      #% we use griddy Gibbs to perform rho-draw
          b0 = solve((t(xs)%*%xs + as.numeric(sige)*TI ),(t(xs)%*%ys + as.numeric(sige)*TIc));
          bd = solve((t(xs)%*%xs + as.numeric(sige)*TI),(t(xs)%*%Wys + as.numeric(sige)*TIc));
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
	if(novi_flag==0){
		if (mm != 0)   rsave[iter-nomit,1] = rval
	}
   	         
    }
                    
iter = iter + 1; 
    
}   ##% end of sampling loop


##% pre-calculate traces for the x-impacts calculations
uiter=50;
maxorderu=100;
nobs = n;
rv=matrix(rnorm(nobs*uiter),nobs,uiter);
tracew=matrix(rep(0,maxorderu),maxorderu,1);
wjjju=rv;
for (jjj in 1:maxorderu){
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
        bdraws = bsave[,2:length(bsave)];
        }
	if(cflag == 0){
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
            total[i,j,] = beta*rmat; ##used beta instead of beta[1,1]
    direct[i,j,] = (beta*trbig)%*%rmat; ##%*%??
    indirect[i,j,] = total[i,j,] - direct[i,j,];
    }

}


##% compute posterior means and log marginal likelihood for return arguments
bmean = mean(bsave);
beta = t(bmean);
rho = mean(psave);
sige = mean(ssave);
vmean = vmean/(ndraw-nomit);
V = In/vmean;


results$sige = sige;
nobs=dim(x)[1]
nvar=dim(x)[2]
xs = matmul(x,sqrt(V));
ys = sqrt(V)*y;
Wys = W%*%ys;  ## sige %*% or * ??
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
   mlike = sar_marginal2(detval,e0,ed,epe0,eped,epe0d,nobs,nvar,logdetx,a1,a2,c,TI,xs,ys,sige,W);
  }
 yhat = solve((diag(nobs) - rho*W),(x*as.numeric(beta)));
 ####ERROR!! : dim(yhat)=dim(x) but yhat should be similar to y.
# e = y - yhat; 

####% compute R-squared
#epe = t(e)%*%e;
#sige = epe/(n-k);
#results$sigma = sige;
#ym = y - mean(y);
#rsqr1 = epe;
#rsqr2 = t(ym)%*%ym;
#results$rsqr = 1- rsqr1/rsqr2; ##% r-squared
#rsqr1 = rsqr1/(nobs-nvar);
#rsqr2 = rsqr2/(nobs-1.0);
#results$rbar = 1 - (rsqr1/rsqr2); ###% rbar-squared


results$meth  = 'sar_g';
results$total = total;
results$direct = direct;
results$indirect = indirect;
results$beta_std = t(apply(bsave,2,sd));
results$sige_std = apply(ssave,2,sd);
results$rho_std = apply(psave,2,sd);
results$beta = beta;
results$rho = rho;
results$bdraw = bsave;
results$pdraw = psave;
results$sdraw = ssave;
results$mlike = mlike;
results$vmean = vmean;
results$yhat  = yhat;
#results$resid = e;
results$bmean = c_beta;
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
return(results)
}#### end of sar_g


