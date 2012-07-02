
#
#Parse values in prior. If prior=NULL return list of default values
#
#prior: list of values
#model: Model. One of "sar", "far", etc.
prior_parse <-function(prior, model,k){

if(model=='far'){
    #% set defaults
	    if(is.null(prior))
	    {
		    prior<-list(
		    eflag = 0,    # % default to not computing eigenvalues
		    ldetflag = 1,  #% default to 1999 Pace and Barry MC determinant approx
		    mflag = 1,     #% default to compute log marginal likelihood
		    order = 50,    #% there are parameters used by the MC det approx
		    iter = 30,     #% defaults based on Pace and Barry recommendation
		    rmin = -1,     #% use -1,1 rho interval as default
		    rmax = 1,
		    detval = 0,    #% just a flag
		    rho = 0.5,
		    sige = 1.0,
		    rval = 4,
		    mm = 0,
		    kk = 0,
		    nu = 0,
		    d0 = 0,
		    a1 = 1.01,
		    a2 = 1.01,
		    cc = 0.2,
		    novi_flag = 0,#% do vi-estimates
		    metflag = 0 #% use integration instead of M-H sampling for rho
    )

		    #VIRGILIO: Add if(model=="sar"), for example, to add specific parameter
		    return(prior)
    }

    # check with user input
    #VIRGILIO: Abhirup, could you modify the code so that 'prior' elements
    #are checked accordingly. As it is now, it seems thatthey are not.
    #You can check for a missing argument using, for example, 
    #if(length(prior$novi)==0) and then set it to the default.
    #Alternatively, you can check for an argument having the wrong lengthas
    #in the previous example:
    #  

    #This is an example:
    #Original code: if(length(prior$novi)!=0) novi_flag=prior$novi
    if(length(prior$novi_flag)!=1)
	    prior$novi_falg<-0



    if(length(prior$nu)!=1) prior$nu=0
    if(length(prior$d0)!=1) prior$d0=0
    if(length(prior$rval)!=1) prior$rval=4
    if(length(prior$dflag)!=1) prior$dflag=0
    if(length(prior$a1)!=1) prior$a1=1.01
    if(length(prior$a2)!=1) prior$a2=1.01
    if(length(prior$m)!=0){ 
	    prior$mm=prior$m
	    prior$kk=prior$k
	    prior$rval=rgamma(1,prior$mm,1/prior$kk)
	    }
    if(length(prior$rmin)!=0) {
	    prior$rmin=prior$rmin
	    prior$eflag=0
	    }
    if(length(prior$rmax)!=0) {
	    prior$rmin=prior$rmax
	    prior$eflag=0
	    }
    if(length(prior$lndet)!=0){
	    prior$detval = prior$lndet;
	    prior$ldetflag = -1;
	    prior$eflag = 0;
	    prior$rmin = prior$detval; ##detval[1,1]
	    prior$nr = length(prior$detval);
	    prior$rmax = prior$detval[nr]; ###detval[nr,1]
	    }
    if(length(prior$lflag)!=0){
	    tst = prior$lflag;
	    if (tst == 0)        prior$ldetflag = 0; 
	    if (tst == 1)        prior$ldetflag = 1; 
	    if (tst == 2)        prior$ldetflag = 2;
	    }
    if(length(prior$order)!=1) prior$order=50
    if(length(prior$iter)!=1) prior$iter=30
    if(length(prior$eig)!=1) prior$eflag=0
    if(length(prior$mlog)!=1) prior$mlog=0


    #return (c(nu,d0,rval,mm,kk,rho,sige,rmin,rmax,detval,ldetflag,eflag,order,iter,novi_flag,cc,metflag,a1,a2))

    #Virgilio: Return the (modified) 'prior' list
    return(prior)
    }
if(model=='sar') {
    if(is.null(prior))
		{
    	prior <-list(	
	eflag = 0,    # % default to not computing eigenvalues
	ldetflag = 1,  #% default to 1999 Pace and Barry MC determinant approx
	mflag = 1,     #% default to compute log marginal likelihood
	order = 50,    #% there are parameters used by the MC det approx
	iter = 30,     #% defaults based on Pace and Barry recommendation
	rmin = -1,     #% use -1,1 rho interval as default
	rmax = 1,
	detval = 0,    #% just a flag
	rho = 0.5,
	sige = 1.0,
	rval = 4,
	mm = 0,
	kk = 0,
	nu = 0,
	d0 = 0,
	a1 = 1.01,
	a2 = 1.01,
	c = matrix(rep(0,k),k,1), #defuse prior for beta
	T = diag(k)*1e+12, ## Abhirup: ??
	prior_beta = 0,   #% flag for diffuse prior on beta
	novi_flag = 0, #% do vi-estimates
	inform_flag = 0)
	return(prior)
	}
	# check with user input
	if(length(prior$novi)!=1) novi_flag=0
	if(length(prior$nu)!=1) prior$nu=0
	if(length(prior$d0)!=1) prior$d0=0
	if(length(prior$rval)!=1) prior$rval=4
	if(length(prior$a1)!=1) prior$a1=1.01
	if(length(prior$a2)!=1) prior$a2=1.01
	if(length(prior$m)!=0){ 
		prior$mm=prior$m
		prior$kk=prior$k
		prior$rval=rgamma(1,prior$mm,1/prior$kk)
		}
	if(length(prior$m)!=0){ 
		prior$c=prior$beta
		prior$inform_flag=1
		}
	if(length(prior$rmin)!=0) {
		prior$rmin=prior$rmin
		prior$eflag=0
		}
	if(length(prior$rmax)!=0) {
		rmin=prior$rmax
		prior$eflag=0
		}
	if(length(prior$lndet)!=0){
		prior$detval = prior$lndet;
		prior$ldetflag = -1;
		prior$eflag = 0;
		prior$rmin = prior$detval; ##detval[1,1]
		prior$nr = length(prior$detval);
		prior$rmax = prior$detval[nr]; ###detval[nr,1]
		}
	if(length(prior$lflag)!=0){
		tst = prior$lflag;
		if (tst == 0)        prior$ldetflag = 0; 
		if (tst == 1)        prior$ldetflag = 1; 
		if (tst == 2)        prior$ldetflag = 2;
		}
	if(length(prior$order)!=1) prior$order=50
	if(length(prior$iter)!=1) prior$iter=30
	if(length(prior$eig)!=1) prior$eflag=prior$eig
	}
	    
    return(prior)
}
