library(Matrix)
source('parse.R')
source('utils.R')

W<-structure(c(0, .5, 0.5, 0.5, 0, 0.5, 0.5, 0.5, 0), .Dim = c(3L, 
3L))

prior = NULL

y<-matrix(rnorm(3), ncol=1)
y<-solve(diag(3)-0.9*W)%*%y
x<-matrix(rnorm(6),ncol=2)

ndraw=100
nomit=20
source('semt_g.R')
results=semt_g(y,x,W,ndraw,nomit,prior)
plot(results$pdraw)
