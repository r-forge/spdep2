#################### Start far_g ####
######
far_g <-function(y,W,ndraw,nomit,prior){

#List to store results
results<-list()
### prior inputs ##

temp=prior_parse(prior)
attach(temp)

#VIRGILIO: Commented out all this because of attach(temp)
#nu=temp[1]
#d0=temp[2]
#rval=temp[3]
#mm=temp[4]
#kk=temp[5]
#rho=temp[6]
#sige=temp[7]
#rmin=temp[8]
#rmax=temp[9]
#detval=temp[10]
#ldetflag=temp[11]
#eflag=temp[12]
#order=temp[13]
#iter=temp[14]
#novi_flag=temp[15]
#cc=temp[16]
#metflag=temp[17]
#a1=temp[18]
#a2=temp[19]
#### end of parse input ##

results$order = prior$order
results$iter = prior$iter

n<-nrow(y)

V = matrix(rep(1,n),n,1)
In = matrix(rep(1,n),n,1)
ys <- y *sqrt(V)#VIRGILIO:Added this which was missing
vi = In

#### Allocate storage for results ####
psave = matrix(rep(0,ndraw-nomit),ndraw-nomit,1)
ssave = matrix(rep(0,ndraw-nomit),ndraw-nomit,1)
vmean = matrix(rep(0,n),n,1)
yhat = matrix(rep(0,n),n,1)
acc_rate = matrix(rep(0,ndraw),ndraw,1)

#### storage for draw on rvalue ####
if (mm!=0) rsave = matrix(rep(0,ndraw-nomit),ndraw-nomit,1)

#tmp = far_eigs(eflag,W,rmin,rmax,n)
tmp = set_eigs(eflag,W,rmin,rmax,n)
rmin = tmp$rmin
rmax = tmp$rmax
#time2 = tmp$time2

#tmp2 = far_lndet(ldetflag, W, rmin, rmax, detval, order, iter)
detval = set_lndet(ldetflag, W, rmin, rmax, detval, order, iter)
#detval = tmp2$detval
#time1 = tmp2$time1

iter=1
#time3 = 0

#### The sampler starts ####

Wy=W%*%y
acc = 0
cc = 0.1

#### For novi_flag=0 or Heteroscedastic Model
if(novi_flag==0){ ##novi_flag_if

#start sampling
while(iter<=ndraw){
	#update sige
	nu1=n+nu
	Wys = sqrt(V)*Wy
	e = ys-rho*Wys
	d1 = d0 +t(e)%*%e
	chi = rchisq(1,nu1)
	t2 = chi/d1
	sige = as.vector(1/t2)#VIRGILIO: FIXED issue with dimension
	
	###update vi
	e = y-rho*Wy
	chiv = rchisq(n,rval+1)
	vi = ((e*e/sige)+In*rval)/chiv
	V = In/vi
	ys = y*sqrt(V)
	
	##update rval
	if(mm!=0){
		rval = rgamma(1,mm,1/kk)
		}
	if(metflag==1){
		#metropolis step to get rho update
		rhox = c_far(rho,y,sige,W,detval,In,a1,a2)
		accept = 0
		rho2 = rho+cc*rnorm(1)
		while(accept==0){
			if((rho2 > rmin)&(rho2 < rmax)){ 
				accept=1
			}
			else rho2 = rho+cc*rnorm(1)
		}
		rhoy = c_far(rho2,y,sige,W,detval,In,a1,a2)
		ru = runif(1)
		if((rhoy-rhox)>exp(1)){
			p=1
		}
		else {
			ratio = exp(rhoy-rhox)
			p = min(ratio)
		}
		if(ru < p){
			rho = rho2
			acc = acc+1
		}
		rtmp[iter,1]=rho
		acc_rate[iter,1]=acc/iter

#### update cc based on sd of rho draws

		if(acc_rate[iter,1]<0.4) cc=cc/1.1
		if(acc_rate[iter,1]>0.6) cc=cc*1.1
	}
	if(metflag==0){
	### update rho using numerical integration ##
	e0 = ys
	ed = Wys
	epe0 = t(e0)%*%e0
	eped = t(ed)%*%ed
	epe0d = t(ed)%*%e0
	rho = draw_rho(detval,epe0,eped,epe0d,n,1,rho,a1,a2)
	}	
	if(iter > nomit){
	ssave[iter-nomit,1] = sige
	psave[iter-nomit,1] = rho
	vmean = vmean+vi
	if(mm!=0) rsave[iter-nomit,1] = rval
	}
iter =iter +1
}
##time??

##compute posterior means and log marginal likelihood for return arguments 

rho = mean(psave)	
sigm = mean(ssave)
vmean = vmean/(ndraw-nomit)
V = In/vmean

ys = sqrt(V)*y
Wys = sqrt(V)*Wy
e0 = ys
ed = Wys
epe0 = t(e0)%*%e0
eped = t(ed)%*%ed
epe0d = t(ed)%*%e0		
e = (e0 - rho*ed)
yhat = y-e
sige = (1/n)*t(e0-rho*ed)%*%(e0-rho*ed)
#mlike = far_marginal(detval,e0,ed,eped,epe0d,n,1,a1,a2)
mlike = rho_marginal(detval,e0,ed, epe0, eped,epe0d,n,1, logdetx=0, a1,a2)

####results

results$y = y      
results$nobs = n
results$nvar = 1   
results$meth = 'far_g'
results$pdraw = psave
results$sdraw = ssave
results$vmean = vmean
results$yhat = yhat
results$resid = e
results$tflag = 'plevel'
results$lflag = ldetflag
results$dflag = metflag
results$nobs  = n
results$ndraw = ndraw
results$nomit = nomit
results$y = y
results$nvar = 1
results$mlike = mlike
results$sige = sige
results$rho = rho
results$lndet = detval
results$acc = acc_rate
results$novi = novi_flag

if(mm != 0){
results$rdraw = rsave
results$m     = mm
results$k     = kk
}
else {
results$r     = rval
results$rdraw = 0
}

#results$time = etime(clock,timet)
#results$time1 = time1
#results$time2 = time2
#results$time3 = time3
results$lndet = detval
results$rmax = rmax 
results$rmin = rmin

} ##end novi_flag_if

if(novi_flag==1){ ##strt novi_flag_if
nu1 = n + nu 
e = y - rho*Wy
d1 = d0 + t(e)%*%e
chi = rchisq(1,nu1)
t2 = chi/d1
sige = 1/t2

	if(metflag==1){
		#metropolis step to get rho update
		rhox = c_far(rho,y,sige,W,detval,In,a1,a2)
		accept = 0
		rho2 = rho+cc*rnorm(1)
		while(accept==0){
			if((rho2 > rmin)&(rho2 < rmax)){ 
				accept=1
			}
			else rho2 = rho+cc*rnorm(1)
		}
		rhoy = c_far(rho2,y,sige,W,detval,In,a1,a2)
		ru = runif(1)
		if((rhoy-rhox)>exp(1)){
			p=1
		}
		else {
			ratio = exp(rhoy-rhox)
			p = min(ratio)
		}
		if(ru < p){
			rho = rho2
			acc = acc+1
		}
		rtmp[iter,1]=rho
		acc_rate[iter,1]=acc/iter
		#### update cc based on sd of rho draws

		if(acc_rate[iter,1]<0.4) cc=cc/1.1
		if(acc_rate[iter,1]>0.6) cc=cc*1.1
	} ###end of metflag==1

	if(metflag==0){
	### update rho using numerical integration ##
	e0 = ys
	ed = Wy
	epe0 = t(e0)%*%e0
	eped = t(ed)%*%ed
	epe0d = t(ed)%*%e0
	rho = draw_rho(detval,epe0,eped,epe0d,n,1,rho,a1,a2)
	} ###end of metflag==0	
	if(iter > nomit){
	ssave[iter-nomit,1] = sige
	psave[iter-nomit,1] = rho
	
	if(mm!=0) rsave[iter-nomit,1] = rval
	}
iter =iter +1
#}


##time??

##compute posterior means and log marginal likelihood for return arguments 

rho = mean(psave)	

e0 = y
ed = Wy
epe0 = t(e0)%*%e0
eped = t(ed)%*%ed
epe0d = t(ed)%*%e0		
e = (e0 - rho*ed)
yhat = y-e
sige = (1/(n-1))*t(e0-rho*ed)%*%(e0-rho*ed)
#mlike = far_marginal(detval,e0,ed,eped,epe0d,n,1,a1,a2)
mlike = rho_marginal(detval,e0,ed,epe0, eped,epe0d,n,1,logdetx=0, a1,a2)#VIRGILIO:FIXED

####results

results$y = y      
results$nobs = n
results$nvar = 1   
results$meth = 'far_g'
results$pdraw = psave
results$sdraw = ssave
results$vmean = vmean
results$yhat = yhat
results$resid = e
results$tflag = 'plevel'
results$lflag = ldetflag
results$dflag = metflag
results$nobs  = n
results$ndraw = ndraw
results$nomit = nomit
results$y = y
results$nvar = 1
results$mlike = mlike
results$sige = sige
results$rho = rho
results$lndet = detval
results$acc = acc_rate
results$novi = novi_flag

if(mm != 0){
results$rdraw = rsave
results$m     = mm
results$k     = kk
}
else {
results$r     = rval
results$rdraw = 0
}

#results$time = etime(clock,timet)
#results$time1 = time1
#results$time2 = time2
#results$time3 = time3
results$lndet = detval
results$rmax = rmax 
results$rmin = rmin

} ##end novi_flag_if

#return (list(results.meth,results.pdraw,results.sdraw,results.vmean,results.rdraw,results.nu,results.d0,results.a1,results.a2,results.r,results.m,results.k,results.nobs,results.ndraw,results.nomit,results.y,results.yhat,results.time,results.time1,results.time2,results.time3,results.rmax,results.rmin,results.tflag,results.lflag,results.dflag,results.iter,results.order,results.limit,results.lndet,results.acc,results.mlike))

detach()
return(results)

}###end of far_g

