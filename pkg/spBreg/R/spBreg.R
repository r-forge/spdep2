spBreg_lag <- function(formula, data = list(), listw, na.action, type="lag",
    zero.policy=NULL, control=list()) {
    timings <- list()
    .ptime_start <- proc.time()
# FIXME
#control
    con <- list(tol.opt=.Machine$double.eps^0.5, ldet_method="SE_classic",
        Imult=2, cheb_q=5, MC_p=16L, MC_m=30L, super=NULL, spamPivot="MMD",
        in_coef=0.1, type="MC", correct=TRUE, trunc=TRUE,
        SE_method="LU", nrho=200, interpn=2000, SElndet=NULL, LU_order=FALSE,
        pre_eig=NULL, interval=c(-1, 1))
    nmsC <- names(con)
    con[(namc <- names(control))] <- control
    if (length(noNms <- namc[!namc %in% nmsC])) 
        warning("unknown names in control: ", paste(noNms, collapse = ", "))
    if (is.null(zero.policy))
        zero.policy <- get.ZeroPolicyOption()
    stopifnot(is.logical(zero.policy))
    if (class(formula) != "formula") formula <- as.formula(formula)
    mt <- terms(formula, data = data)
    mf <- lm(formula, data, na.action=na.action,  method="model.frame")
    na.act <- attr(mf, "na.action")
    if (!inherits(listw, "listw")) stop("No neighbourhood list")
    can.sim <- FALSE
    if (listw$style %in% c("W", "S")) can.sim <- can.be.simmed(listw)
    if (!is.null(na.act)) {
        subset <- !(1:length(listw$neighbours) %in% na.act)
        listw <- subset(listw, subset, zero.policy=zero.policy)
    }
    y <- model.extract(mf, "response")
    x <- model.matrix(mt, mf)
    n <- nrow(x)
    if (n != length(listw$neighbours))
        stop("Input data and weights have different dimensions")
    xcolnames <- colnames(x)
    wy <- lag.listw(listw, y, zero.policy=zero.policy)
    if (anyNA(wy)) stop("NAs in lagged dependent variable")
#create_WX
    if (type == "Durbin") {
        WX <- create_WX(x, listw, zero.policy=zero.policy, prefix="lag")
        x <- cbind(x, WX)
        rm(WX)
    } else if (type != "lag") stop("No such type:", type)
    m <- ncol(x)
    lm.base <- lm(y ~ x - 1)
    aliased <- is.na(coefficients(lm.base))
    cn <- names(aliased)
    names(aliased) <- substr(cn, 2, nchar(cn))
    if (any(aliased)) {
        nacoef <- which(aliased)
        x <- x[,-nacoef]
    }
    m <- ncol(x)
    timings[["set_up"]] <- proc.time() - .ptime_start
    .ptime_start <- proc.time()

    env <- new.env()
    assign("can.sim", can.sim, envir=env)
    assign("listw", listw, envir=env)
    assign("similar", FALSE, envir=env)
    assign("n", n, envir=env)
    interval <- spdep::jacobianSetup(con$ldet_method, env, con,
        pre_eig=con$pre_eig, interval=con$interval)
    nmsC <- names(con)
    con[(namc <- names(control))] <- control
    if (length(noNms <- namc[!namc %in% nmsC])) 
        warning("unknown names in control: ", paste(noNms, collapse = ", "))
    if (is.null(zero.policy))
        zero.policy <- get.ZeroPolicyOption()
    stopifnot(is.logical(zero.policy))
    if (class(formula) != "formula") formula <- as.formula(formula)
    mt <- terms(formula, data = data)
    mf <- lm(formula, data, na.action=na.action,  method="model.frame")
    na.act <- attr(mf, "na.action")
    if (!inherits(listw, "listw")) stop("No neighbourhood list")
    can.sim <- FALSE
    if (listw$style %in% c("W", "S")) can.sim <- can.be.simmed(listw)
    if (!is.nullcon$interval=interval)
    assign("interval", interval, envir=env)

    nm <- paste(con$ldet_method, "set_up", sep="_")
    timings[[nm]] <- proc.time() - .ptime_start
    .ptime_start <- proc.time()



#bounds
#ldets handling

#homoskedastic
#dflag MH/integration

#run loop(s)

#sar_g <- function(y,x,W,ndraw,nomit,prior){


#% error checking on inputs
n=nrow(y)
junk = ncol(y)
n1=nrow(x)
k=ncol(x)
n2=nrow(W)
n4=ncol(W)


results <-list()
results$nobs  = n;
results$nvar  = k;
results$y = y; 
  
    
# FIXME pprior=prior_parse(prior,k)
# FIXME attach(pprior)
#n = nrow(y)
     
results$cflag = cflag;
results$p = p;
    

results$order = order;
results$iter = iter;

#timet = clock; % start the timer

out_temp = set_eigs(eflag,W,rmin,rmax,n)
rmin=out_temp$rmin
rmax=out_temp$rmax
#time1=out_temp$time1
#results$time1 = time1;

detval = set_lndet(ldetflag,W,rmin,rmax,detval,order,iter)
#detval=out_temp$detval
#time2=out_temp$time2
#results$time2 = time2;

results$order = order;
results$iter = iter;


#% storage for draws
bsave = matrix(0,ndraw-nomit,k)
if (mm !=0) 
	rsave = matrix(0,ndraw-nomit,1)

psave = matrix(0, ndraw-nomit,1)
ssave = matrix(0, ndraw-nomit,1)
vmean= matrix(0, n,1)


#% ====== initializations
#% compute this stuff once to save time

TI = solve(Tbeta); # see eq 5.29, p. 140
TIc = TI%*%c_beta;

In = matrix(1,n,1);
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
          AI = solve((t(xs)%*%xs + sige*TI))#,diag(rep(1,k)));         
          yss = ys - rho*Wys;          
          xpy = t(xs)%*%yss;
          b = t(xs)%*%yss + sige*TIc;
          b0 = solve((t(xs)%*%xs + sige*TI),b); # see eq 5.29, p. 140
          #bhat = norm_rnd(sige*AI) + b0; ##code for norm_rnd ?? check  
          bhat = mvrnorm(1, b0, sige*AI)##From MASS package
          xb = xs%*%bhat; 
                    
          #% update sige
          nu1 = n + 2*nu; 
          e = (yss - xb);
          d1 = 2*d0 + t(e)%*%e;
          chi = rchisq(1,nu1);
          sige = as.numeric(d1/chi); # see eq 5.30, p. 141

	  #% update vi
          ev = y - rho*Wy - x%*%bhat; 
          #chiv = chis_rnd(n,rval+1);  
          chiv = matrix(rchisq(n,rval+1),n,1);
          vi = (ev*ev/sige + In%*%rval)/chiv;
          V = In/vi; 
                        
          #% update rval
          if (mm != 0)   
		rval = rgamma(1, shape=mm, rate=kk)#gamm_rnd(1,1,mm,kk);  
          
      #% we use griddy Gibbs to perform rho-draw
#          b0 = (t(xs)%*%xs + sige%*%TI )\(t(xs)%*%ys + sige%*%TIc);
#          bd = (t(xs)%*%xs + sige%*%TI)\(t(xs)%*%Wys + sige%*%TIc);
          b0 = solve( (t(xs)%*%xs + sige*TI ), (t(xs)%*%ys + sige*TIc) );
          bd = solve( (t(xs)%*%xs + sige*TI), (t(xs)%*%Wys + sige*TIc) );
          e0 = ys - xs%*%b0;
          ed = Wys - xs%*%bd;
          epe0 = t(e0)%*%e0;
          eped = t(ed)%*%ed;
          epe0d = t(ed)%*%e0;
          rho = draw_rho(detval,epe0,eped,epe0d,n,k,rho,a1,a2);  #function later

          
    if (iter > nomit){ ## if we are past burn-in, save the draws
        bsave[iter-nomit,1:k] = as.vector(bhat);
        ssave[iter-nomit,1] = as.vector(sige);
        psave[iter-nomit,1] = as.vector(rho);
        vmean = vmean + vi; 
    	if (mm != 0)
		rsave[iter-nomit,1] = as.vector(rval)
    	         
    }

	iter = iter + 1; 
	}
                    
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
          AI = solve((xpx + sige*TI))#,diag(rep(1,k)));        
          ys = y - rho*Wy;          
          b = t(x)%*%ys + sige*TIc;
          b0 = solve((xpx + sige*TI),b);# see eq 5.29, p. 140
          bhat = mvrnorm(1, b0, sige*AI); #norm_rnd(sige*AI) + b0;  
          xb = x%*%bhat;
          
          ##% update sige
          nu1 = n + 2*nu; 
          ##%e = e0 - rho*ed;
          e = (ys - xb);
          d1 = 2*d0 + t(e)%*%e;
          #chi = chis_rnd(1,nu1);
          chi = rchisq(1,nu1);
          sige = as.numeric(d1/chi); # see eq 5.30, p. 141
          
          ###% update rho using griddy Gibbs
          AI = solve((xpx + sige*TI))#,diag(rep(1,k)));
          b0 = solve((xpx + sige*TI),(xpy + sige*TIc));
          bd = solve((xpx + sige*TI),(xpWy + sige*TIc));
          e0 = y - x%*%b0;
          ed = Wy - x%*%bd;
          epe0 = t(e0)%*%e0;
          eped = t(ed)%*%ed;
          epe0d = t(ed)%*%e0;
          rho = draw_rho(detval,epe0,eped,epe0d,n,k,rho,a1,a2);

    if (iter > nomit){ ###% if we are past burn-in, save the draws
    bsave[iter-nomit,1:k] = as.vector(bhat);
    ssave[iter-nomit,1] = as.vector(sige);
    psave[iter-nomit,1] = as.vector(rho);
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
tracew=matrix(0, maxorderu,1);
wjjju=rv;

for (jjj in 1:maxorderu){
    wjjju=W%*%wjjju;
#    tracew[jjj]=mean(mean(rv*wjjju));
    tracew[jjj, 1]=mean(apply(rv*wjjju, 2, mean));
    
}

traces=tracew;
traces[1,1]=0;
traces[2,1]=sum(t(W)*W)/nobs;
trs=matrix(c(1,traces), ncol=1);
ntrs=length(trs);
trbig=t(trs);
                 
        if (cflag == 1){
        bdraws = bsave[,-1];#2:length(bsave)];
        }
	if(cflag == 0){
        bdraws = bsave;
        }
        pdraws = psave;

        ree = 0:(ntrs-1);

        rmat = matrix(0, 1,ntrs);####three dimentional matrix in R??
#        total = matrix(rep(0,(ndraw-nomit)*p*ntrs)),ndraw-nomit,p,ntrs);
        total = array(0, dim=c(ndraw-nomit,p,ntrs));
#        direct = zeros(ndraw-nomit,p,ntrs);
        direct = array(0,dim=c(ndraw-nomit,p,ntrs));
#        indirect = zeros(ndraw-nomit,p,ntrs);
        indirect = array(0,dim=c(ndraw-nomit,p,ntrs));
       
for (i in 1:(ndraw-nomit))
{
    rmat = pdraws[i,1]^ree;#FIXME: Check this
    for (j in 1:p)
	{
            bbeta = bdraws[i,j];#VIRGILIO: Check this
            total[i,j,] = bbeta*rmat;
    direct[i,j,] = (bbeta*trbig)*rmat;
    indirect[i,j,] = total[i,j,] - direct[i,j,];
    }

}
#### ?? ###
#time4 = etime(clock,t0);
#results.time4 = time4;


#% compute posterior means and log marginal likelihood for return arguments
bmean = apply(bsave, 2, mean);
beta = matrix(bmean, ncol=1);
rho = mean(psave);
sige = mean(ssave);
vmean = vmean/(ndraw-nomit);
V = In/vmean;


results$sige = sige;
nobs=nrow(x)
nvar=ncol(x)
#[nobs,nvar] = size(x);
          xs = matmul(x,sqrt(V));
          ys = sqrt(V)*y;
          Wys = W%*%ys;
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
   mlike = rho_marginal2(detval,e0,ed,epe0,eped,epe0d,nobs,nvar,logdetx,a1,a2,c_beta,TI,xs,ys,sige,W);
  }
 yhat = solve((diag(nobs) - rho*W),(x%*%beta));
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
results$beta_std = apply(bsave, 2, sd);
results$sige_std = apply(ssave, 2, sd);
results$rho_std = apply(psave, 2, sd);
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
results$bstd  = sqrt(diag(Tbeta));
results$ndraw = ndraw;
results$nomit = nomit;
#results$time  = time;
results$nu = nu;
results$d0 = d0;
results$a1 = a1;
results$a2 = a2;
results$tflag = "plevel";
results$rmax = rmax; 
results$rmin = rmin;
results$lflag = ldetflag;
results$lndet = detval;
results$novi  = novi_flag;
results$priorb = inform_flag;

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


# FIXME detach(pprior)
# FIXME return(results)
#### end of sar_g


#output mcmc object
#do impacts externally
#do summary externally
}
