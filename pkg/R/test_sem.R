## test file

library(Matrix)
source('parse.R')
source('utils.R')

W<-structure(c(0, .6, 0.4, 0.5, 0, 0.5, 0.8, 0.2, 0), .Dim = c(3L, 
3L))

prior = NULL
#prior<-prior_parse(NULL,2)
#prior$ldetflag<-0
#prior$rmin<-0.1
#prior$rmax<-0.99#If this is set to 1 then problem when interpolating log-det's

#prior$metflag<-0#M-H sampler

y<-matrix(rnorm(3), ncol=1)
y<-solve(diag(3)-0.9*W)%*%y
x<-matrix(rnorm(6),ncol=2)

#prior$novi_flag<-0
ndraw=150
nomit=30
source('sem_g.R')
results<-sem_g(y, x, W, 500, 10, prior)
plot(results$pdraw, type="l", main="rho")
summary(results$pdraw)



##Bigger example
library(spdep)
data(boston)

W<-nb2mat(boston.soi)
n<-nrow(W)
x<-cbind(rep(1, n), rnorm(n))


y<-cbind(x, W%*%x[,-1])%*%matrix(c(0,5,10), ncol=1)
y<-matrix(y+rnorm(n), ncol=1)
y<-solve(diag(nrow(W))-0.4*W)%*%y


#Heteroc.  model
prior<-prior_parse(NULL, k=1)
prior$novi_flag<-0
results<-sem_g(y, x, W, 150, 50, prior)
plot(results$pdraw, type="l", main="rho")
summary(results$pdraw)



##Homoc. model
prior$novi_flag<-1
results1<-sem_g(y, x, W, 150, 50, prior)
plot(results1$pdraw, type="l", main="rho")
summary(results1$pdraw)


