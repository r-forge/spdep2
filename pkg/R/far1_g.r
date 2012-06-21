# An implementation of Bayesian FAR model using Metropolis within gibbs method

W <-read.table('wmat.dat')
far1_g <-function(W,ndraw=1100,nomit=100){
n = length(W[1,])
nadj <-ndraw-nomit
IN <-diag(rep(1,n))
In <-matrix(rep(1,n),n,1)
weig_tmp <-eigen(W)
weig <-weig_tmp$values
#bounds on rho
lmin <-1/min(weig)
lmax <-1/max(weig)
rho <-0.7
y <-solve(IN-rho*W)%*%matrix(rnorm(n),n,1)
ydev <-y-mean(y)
Wy <-W*ydev

#### set starting values ####
rho <-0.5
sige <-10.0
c <-0.5
rsave <-matrix(rep(0,nadj),nadj,1)
ssave <-matrix(rep(0,nadj),nadj,1)
rtmp <-matrix(rep(0,nadj),nomit,1)
iter <-1
cnt <-0
while (iter <= ndraw){
	e <-ydev-rho*Wy
	ssr <-t(e)%*%e
	chi <-rchisq(1,df=n)
	sige <-ssr/chi

#### Metropolis step to get rho update ####

	rhox <-c_rho(rho,sige,ydev,W)
	rho2 <-rho +c*rnorm(1)
	accept <-0
	while (accept==0){
		if((rho2>lmin)&(rho2<lmax)) accept <-1
		rho2 <- rho+c*rnorm(1)
		cnt <-cnt+1
	}
	rhoy <-c_rho(rho2,sige,ydev,W)
	ru <-runif(1,0,1)
	ratio <-rhoy/rhox
	p <-min(1,ratio)
	if(ru < p) rho=rho2
	rtmp[iter,1] <-rho
	if(iter >=nomit){
		if(iter ==nomit) c <-2*sd(rtmp[1:nomit,1])
	}
	ssave[iter-nomit+1,1] <- sige
	rsave[iter-nomit+1,1] <-rho

iter <-iter+1
}
#### Results ####
rho_mean <-apply(rsave,2,mean)
rho_sd <-apply(rsave,2,sd)
sig_mean <-apply(ssave,2,mean)
sig_sd <-apply(ssave,2,sd)
hit_rate <-ndraw/cnt

output <- list(hit_rate=hit_rate,rho_mean=rho_mean,rho_sd=rho_sd,sig_mean=sig_mean,sig_sd=sig_sd)
output
}











