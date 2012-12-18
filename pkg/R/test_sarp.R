## test file

library(Matrix)
library(MASS)
library(msm)
library(spdep)
source('parse.R')
source('utils.R')

#W<-structure(c(0, .5, 0.5, 0.5, 0, 0.5, 0.5, 0.5, 0), .Dim = c(3L, 3L))

prior = NULL
#prior<-prior_parse(NULL,2)
#prior$ldetflag<-0
#prior$rmin<-0.1
#prior$rmax<-0.99#If this is set to 1 then problem when interpolating log-det's

#prior$metflag<-0#M-H sampler

nsim<-100
xx<-matrix(rnorm(nsim*2), ncol=2)
W<-nb2mat(knn2nb(knearneigh(xx, 2)))

err<-matrix(1*rnorm(nsim), ncol=1)
x<-cbind(1, matrix(rnorm(2*nsim),ncol=2))
#xWx<-cbind(x, W%*%x)
ystar<-solve(diag(nsim)-0.2*W, x%*%matrix(c(1,5, 10), ncol=1)+err)
y<-matrix(as.numeric(ystar>0), ncol=1)


ndraw=2000
nomit=1000

source('sarp_g.R')

prior<-prior_parse(NULL, k=ncol(x))
prior$rmin<-.01
prior$rmax<-.99
prior$novi_flag<-1

results=sarp_g(y,x,W,ndraw,nomit,prior)
plot(density(results$pdraw), type="l", main="rho")


#Fit model to latent variable
source('sar_g.R')
results2=sar_g(ystar,x,W,ndraw,nomit,prior)
plot(density(results2$pdraw), type="l", main="rho")



