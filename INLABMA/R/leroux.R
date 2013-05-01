#
#Code to fit some the model by Leroux et al. (1999)
#


#Simultaneous Autorregresive Model
#
#Return an INLA model
#Formula is for the FIXED effects ONLY
#'...' are passed to the INLA call
#areaid is the column to be used as area index in the random spatial effect
#

#W is a binary adjacency matrix
#lambda is the parameter in the mixture of precision matrices
#

leroux.inla<-function(formula, d, W, lambda, improve=TRUE, impacts=FALSE, fhyper=NULL, probit=FALSE,...)
{

	W2<-diag(apply(W, 1, sum))-W
	Q<-lambda*diag(nrow(W))+(1-lambda)*W2

	environment(formula)<-environment()

	if(is.null(fhyper))
	{
	formula<-update(formula, 
	  . ~ . + f(idx, model="generic0", Cmatrix =Q) )
	}
	else
	{
	formula<-update(formula, 
	  . ~ . + f(idx, model="generic0", Cmatrix =Q, hyper=fhyper) )
	}


	res<-inla(formula, data=d, ...)

	if(improve)
		res<-inla.hyperpar(res, diff.logdens=20)

	#Compute log-determinat to correct the marginal-loglikelihood
	res$logdet<-as.numeric(determinant(Q)$modulus)
	res$mlik<-res$mlik+res$logdet/2

	return(res)
}
