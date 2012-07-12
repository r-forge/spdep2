## test file

library(Matrix)
source('parse.R')
source('utils.R')

W1<-structure(c(0, .6, 0.4, 0.5, 0, 0.5, 0.8, 0.2, 0), .Dim = c(3L, 
3L))
W2<-structure(c(0, .5, 0.5, 0.7, 0, 0.3, 0.8, 0.2, 0), .Dim = c(3L, 
3L))

prior = NULL
#prior<-prior_parse(NULL,2)
#prior$ldetflag<-0
#prior$rmin<-0.1
#prior$rmax<-0.99#If this is set to 1 then problem when interpolating log-det's

#prior$metflag<-0#M-H sampler

y<-matrix(rnorm(3), ncol=1)
y<-solve(diag(3)-0.9*W1)%*%y
x<-matrix(rnorm(6),ncol=2)

#prior$novi_flag<-0
ndraw=10
nomit=2
source('sac_g.R')
results<-sac_g(y, x, W1,W2, 10, 2, prior)
plot(results$pdraw, type="l", main="rho")
summary(results$pdraw)
