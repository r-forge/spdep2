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


library(spdep)
library(INLA)
library(multicore)


#Load data
data(boston)
boston.c$idx<-1:nrow(boston.c)#Create index for INLA

boston.mat<-nb2mat(boston.soi)
boston.matB<-listw2mat(nb2listw(boston.soi, style="B"))
bmspB<-as(boston.matB, "CsparseMatrix")

#Some stuff we need later
#fhyper=list(prec=list(prior="loggamma", param=c(1,1), initial=log(1), fixed=TRUE))
fhyper=list(prec=list(prior="loggamma", param=c(1,1), initial=log(1), fixed=FALSE))

zero.variance = list(prec=list(initial = 25, fixed=TRUE))

form<-log(CMEDV) ~CRIM+ZN + INDUS + CHAS + I(NOX^2)+
   I(RM^2) +  AGE + log(DIS) + log(RAD) + TAX + PTRATIO + B + log(LSTAT)


#Run models
rlambda<-seq(.01, .2, length.out=10)

lerouxmodels = mclapply(rlambda,
        function(lambda) {
                leroux.inla(form, d=boston.c, W=bmspB, lambda=lambda,
                        fhyper=fhyper,
                        family = "gaussian",
   control.family = list(hyper = zero.variance),
                        control.predictor=list(compute=TRUE),
                        control.compute=list(dic=TRUE, cpo=TRUE),
#control.inla=list(print.joint.hyper=TRUE, strategy="laplace"), verbose=FALSE
                        control.inla=list(print.joint.hyper=TRUE), verbose=FALSE
                )
        })


plot(rlambda, unlist(lapply(lerouxmodels, function(X){X$mlik[1,1]})), type="l")

source("BMArho.R")
bmaleroux<-BMA2(lerouxmodels, rlambda, 0, impacts=FALSE)


#Run with CARBAYES
library(CARBayes)

attach(boston.c)
lcarbayes<- gaussian.lerouxCAR (form, W=as.matrix(bmspB), burnin=500, n.sample=1500)
detach(boston.c)


#Compare results between INLA and CARBayes
plot(density(lcarbayes$samples.rho))
lines(1-bmaleroux$rho$marginal[,1], bmaleroux$rho$marginal[,2], col="red")



