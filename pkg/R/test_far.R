## test file

library(Matrix)
source('parse.R')
source('utils.R')

W<-structure(c(0, .5, 0.5, 0.5, 0, 0.5, 0.5, 0.5, 0), .Dim = c(3L, 
3L))

prior = NULL
#prior<-prior_parse(NULL,2)
#prior$ldetflag<-0
#prior$rmin<-0.1
#prior$rmax<-0.99#If this is set to 1 then problem when interpolating log-det's

#prior$metflag<-0#M-H sampler

nsim<-3
y<-3*matrix(rnorm(nsim), ncol=1)
y<-solve(diag(nsim)-0.9*W, y)


ndraw=1000
nomit=500

source('far_g.R')

prior<-prior_parse(NULL, k=0)
prior$rmin<-.01
prior$rmax<-.99
prior$novi_flag<-1

results=far_g(y,W,ndraw,nomit,prior)
plot(density(results$pdraw), type="l", main="rho")
summary(results$pdraw)
summary(sqrt(results$sdraw))



#Bigger example
library(spdep)
data(boston)

W<-nb2mat(boston.soi)
n<-nrow(W)


y<-3*matrix(rnorm(n), ncol=1)
y<-solve(diag(nrow(W))-0.4*W)%*%y


#Heteroc.  model
prior<-prior_parse(NULL, k=0)
prior$novi_flag<-0
results<-far_g(y, W, 1000, 500, prior)
plot(density(results$pdraw), type="l", main="rho")
summary(results$pdraw)
summary(sqrt(results$sdraw))



#Homoc. model
prior$novi_flag<-1
results1<-far_g(y, W, 1000, 500, prior)
plot(density(results1$pdraw), type="l", main="rho")
summary(results1$pdraw)
summary(sqrt(results1$sdraw))




#Example with Boston housing data
#Code from help(boston, package="spdep")
y<-matrix(log(boston.c$MEDV), ncol=1)

prior<-prior_parse(NULL, k=0)
prior$novi_flag<-1
prior$metflag<-1

resultsx<-far_g(y, W, 300, 100, prior)
#plot(results$pdraw, type="l", main="rho")
plot(density(resultsx$pdraw), type="l", main="rho")
summary(resultsx$pdraw)

#With lagged covariates
prior<-prior_parse(NULL, k=0)
prior$novi_flag<-0
prior$metflag<-1
resultsxlag<-far_g(y, W, 1000, 500, prior)
#plot(results$pdraw, type="l", main="rho")
plot(density(resultsxlag$pdraw), type="l", main="rho")
summary(resultsxlag$pdraw)


