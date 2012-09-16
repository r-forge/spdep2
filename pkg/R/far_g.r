#################### Start far_g ####
######
far_g <-function(y,W,ndraw,nomit,prior){

results<-list()
### prior inputs ##

temp=prior_parse(prior)
attach(temp)

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
rtmp = matrix(rep(0,ndraw-nomit),ndraw-nomit,1)
vmean = matrix(rep(0,n),n,1)
yhat = matrix(rep(0,n),n,1)
acc_rate = matrix(rep(0,ndraw),ndraw,1)

#### storage for draw on rvalue ####
if (mm!=0) rsave = matrix(rep(0,ndraw-nomit),ndraw-nomit,1)

#tmp = far_eigs(eflag,W,rmin,rmax,n)
tmp = set_eigs(eflag,W,rmin,rmax,n)
rmin = tmp$rmin
rmax = tmp$rmax
detval = set_lndet(ldetflag, W, rmin, rmax, detval, order, iter)


iter=1

Wy=W%*%y
Wys = sqrt(V)*Wy
acc = 0
cc = 0.1

while(iter<=ndraw){

	#update sige
	nu1=n+nu
	Wys = sqrt(V)*Wy
	e = ys-rho*Wys
	d1 = d0 +t(e)%*%e
	chi = rchisq(1,nu1)
	t2 = chi/d1
	sige = as.vector(1/t2)#VIRGILIO: FIXED issue with dimension
	
	e = y-rho*Wy

	###update vi  and rval IF novi_flag==0
	if(novi_flag==0)
	{
		chiv = rchisq(n,rval+1)
		vi = ((e*e/sige)+In*rval)/chiv
		V = In/vi
		ys = y*sqrt(V)
	##update rval
	if(mm!=0){
		rval = rgamma(1,mm,1/kk)
		}
	}


	if(metflag==1){
		#metropolis step to get rho update
		rhox = c_rho(rho,y,sige,W,detval,In,a1,a2)
		accept = 0
		rho2 = rho+cc*rnorm(1)
		while(accept==0){
			if((rho2 > rmin)&(rho2 < rmax)){ 
				accept=1
			}
			else rho2 = rho+cc*rnorm(1)
		}
		rhoy = c_rho(rho2,y,sige,W,detval,In,a1,a2)

		#print(c(rho, rhox, rho2, rhoy))

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
	
	#print(c(rho, iter, nomit))
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

results$lndet = detval
results$rmax = rmax 
results$rmin = rmin


detach()
return(results)

}###end of far_g

