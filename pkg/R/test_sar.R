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

y<-matrix(rnorm(3), ncol=1)
y<-solve(diag(3)-0.9*W)%*%y
x<-matrix(rnorm(6),ncol=2)

#prior$novi_flag<-0
ndraw=1000
nomit=0
source('sar_g.R')
results=sar_g(y,x,W,ndraw,nomit,prior)
plot(results$pdraw)

