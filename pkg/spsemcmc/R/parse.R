prior_parse <- function(prior,k=2){
if(is.null(prior)){
		prior<-list(
	    eflag = 0,    # % default to not computing eigenvalues
	    ldetflag = 1,  #% default to 1999 Pace and Barry MC determinant approx
	    mflag = 1,     #% default to compute log marginal likelihood
	    order = 50,    #% there are parameters used by the MC det approx
	    iter = 30,     #% defaults based on Pace and Barry recommendation
	    rmin = -1,     #% use -1,1 rho interval as default
	    rmax = 1,
	    lmin=-1,
	    lmax=1,
	    detval = 0,    #% just a flag
	    detval1=1,
	    detval2=1,
	    rho = 0.5,
	    lambda=0.2,
	    sige = 1.0,
	    rval = 4,
	    mm = 1,#changed gamma prior defaults
	    kk = 1,
	    nu = 0,
	    d0 = 0,
	    a1 = 1.01,
	    a2 = 1.01,
	    cc = 0.2,
	    novi_flag = 0,#% do vi-estimates
	    metflag = 0, #% use integration instead of M-H sampling for rho
	    nsample=1,
	    #for sar
	    c_beta=matrix(rep(0,k),k,1), #c is changed to c_beta for c() is a fun.
	    Tbeta=diag(k)*1e+12,
	    inform_flag=0,
	    #for sac
	    lambda=0.2,
	    lmin=-1,
	    lmax=1,
	    detval1=0,
	    detval2=0,
	    cc1=0.2,
	    cc2=0.2,
	    #for sem
	    mlog=1,	    
	    #for sart
	    cflag=0,
	    vflag=0.0
    )
	}
else {
		if(length(prior$novi_flag)!=1)	prior$novi_flag=0
		if(length(prior$eflag)!=1) prior$eflag=0
		if(length(prior$ldetflag)!=1) prior$ldetflag=1
		if(length(prior$mflag)!=1) prior$mflag=1
		if(length(prior$nu)!=1) prior$nu=0
		if(length(prior$d0)!=1) prior$d0=0
		if(length(prior$rval)!=1) prior$rval=4
		if(length(prior$dflag)!=1) prior$dflag=0
		if(length(prior$rho)!=1) prior$rho=0.5
		if(length(prior$lambda)!=1) prior$lambda=0.2
		if(length(prior$sige)!=1) prior$sige=1.0
		if(length(prior$a1)!=1) prior$a1=1.01
		if(length(prior$a2)!=1) prior$a2=1.01
		if(length(prior$cc)!=1) prior$cc=0.2
		if(length(prior$cc1)!=1) prior$cc1=0.2
		if(length(prior$cc2)!=1) prior$cc2=0.2
		if(length(prior$metflag)!=1) prior$metflag=0
		if(length(prior$limit)!=1) prior$vflag=0.0
		if(length(prior$nsample)!=1) prior$nsample=1
		if(length(prior$m)!=1) prior$mm=1
		if(length(prior$k)!=1) prior$kk=1
		if(is.null(prior$c_beta))
			{ 
			prior$c_beta=matrix(rep(0,k),k,1)
	    		Tbeta=diag(k)*1e+12
			}
			else 
			{
				prior$c_beta=prior$c_beta	
	    			prior$Tbeta=prior$Tbeta
			}
		if(prior$k!=0){
			prior$rval=rgamma(1,prior$mm,1/prior$kk)
			}
		if(length(prior$rmin)!=1) {
			prior$rmin=-1	#eflag =0 ?
			prior$eflag=0
			}
		if(length(prior$rmax)!=1) {
			prior$rmax=1
			prior$eflag=0
			}
		if(length(prior$lmin)!=1) {
			prior$lmin=-1	#eflag =0 ?
			prior$eflag=0
			}
		if(length(prior$lmax)!=1) {
			prior$lmax=1
			prior$eflag=0
			}	
#VIRGILIO: What exactly is this???
#		if(length(prior$beta)!=1) {
#			prior$c_beta=prior$beta
#			prior$inform_flag=1
#			}
#		if(length(prior$bcov)!=1) {
#			prior$Tbeta=prior$bcov
#			prior$inform_flag=1
#			}
		if(length(dim(prior$lndet))==2) {
			prior$detval = prior$lndet;
			prior$ldetflag = -1;
			prior$eflag = 0;
			prior$rmin = prior$detval[1,1]; ##detval[1,1]
			nr = length(prior$detval);
			prior$rmax = prior$detval[nr,1];
			}
		if(length(dim(prior$lndet1))==2){
			prior$detval1 = prior$lndet1;
			prior$ldetflag = -1;
			prior$eflag = 0;
			prior$rmin = prior$detval1[1,1]; ##detval[1,1]
			nr = length(prior$detval1);
			prior$rmax = prior$detval1[nr,1];
			}
		if(length(dim(prior$lndet2))==2){
			prior$detval2 = prior$lndet2;
			prior$ldetflag = -1;
			prior$eflag = 0;
			prior$lmin = prior$detval2[1,1]; ##detval[1,1]
			nr = length(prior$detval2);
			prior$lmax = prior$detval[nr,1];
			}	
		## ?: lflag info not getting passed ??
			if(length(prior$order)!=1) prior$order=50
			if(length(prior$iter)!=1) prior$iter=30
			if(length(prior$eflag)!=1) prior$eflag=0

		

}

	return(prior)
}	
