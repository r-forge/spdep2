##demo file

library(Matrix)
source('parse.R')
source('utils.R')
source('far_g.R')
source('sac_g.R')
source('sar_g.R')
source('sem_g.R')
source('sdm_g.R')

#########initialisations######
W<-structure(c(0, .5, 0.5, 0.5, 0, 0.5, 0.5, 0.5, 0), .Dim = c(3L, 
3L))

W1<-structure(c(0, .6, 0.4, 0.5, 0, 0.5, 0.8, 0.2, 0), .Dim = c(3L, 
3L))
W2<-structure(c(0, .5, 0.5, 0.7, 0, 0.3, 0.8, 0.2, 0), .Dim = c(3L, 
3L))

prior = NULL

y<-matrix(rnorm(3), ncol=1)
y<-solve(diag(3)-0.9*W1)%*%y
x<-matrix(rnorm(6),ncol=2)

ndraw=10
nomit=2
#############

demo_set <- function(method){
if(method=="far")
	{
	results<-far_g(y, W, ndraw, nomit, prior)
	}
if(method=="sar")
	{
	results<-sar_g(y,x,W,ndraw,nomit,prior)
	}
if(method=="sac")
	{
	results<-sac_g(y, x, W1,W2, ndraw,nomit, prior)
	}
if(method=="sdm")
	{
	results<-sdm_g(y,x,W,ndraw,nomit,prior)
	}
if(method=="sem")
	{
	results<-sem_g(y, x, W, ndraw, nomit, prior)
	}
	
plot(results$pdraw, type="l", main="rho")
summary(results$pdraw)	
}

