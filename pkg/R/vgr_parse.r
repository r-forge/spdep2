
#
#Parse values in prior. If prior=NULL return list of default values
#
#prior: list of values
#model: Model. One of "sar", "far", etc.
prior_parse <-function(prior, model){

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

		#VIRGILIO: Add if(model=="sar"), for example, to add specific parameters


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



if(length(prior$nu)!=0) nu=prior$nu
if(length(prior$d0)!=0) d0=prior$d0
if(length(prior$rval)!=0) rval=prior$rval
if(length(prior$dflag)!=0) metflag=prior$dflag
if(length(prior$nu)!=0) nu=prior$nu
if(length(prior$a1)!=0) a1=prior$a1
if(length(prior$a2)!=0) a2=prior$a2
if(length(prior$m)!=0){ 
	mm=prior$m
	kk=prior$k
	rval=rgamma(1,mm,1/kk)
	}
if(length(prior$rmin)!=0) {
	rmin=prior$rmin
	eflag=0
	}
if(length(prior$rmax)!=0) {
	rmin=prior$rmax
	eflag=0
	}
if(length(prior$lndet)!=0){
	detval = prior$lndet;
	ldetflag = -1;
	eflag = 0;
	rmin = detval; ##detval[1,1]
	nr = length(detval);
	rmax = detval[nr]; ###detval[nr,1]
	}
if(length(prior$lflag)!=0){
	tst = prior$lflag;
        if (tst == 0)        ldetflag = 0; 
        if (tst == 1)        ldetflag = 1; 
        if (tst == 2)        ldetflag = 2;
	}
if(length(prior$order)!=0) order=prior$order
if(length(prior$iter)!=0) iter=prior$iter
if(length(prior$eig)!=0) eflag=prior$eig
if(length(prior$mlog)!=0) mlog=prior$mlog


#return (c(nu,d0,rval,mm,kk,rho,sige,rmin,rmax,detval,ldetflag,eflag,order,iter,novi_flag,cc,metflag,a1,a2))

#Virgilio: Return the (modified) 'prior' list
return(prior)
}
