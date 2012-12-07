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
x<-matrix(rnorm(6),ncol=2)
y<-solve(diag(3)-0.9*W)%*%(y+x%*%matrix(c(0,5), ncol=1))

#prior$novi_flag<-0
ndraw=150
nomit=30
source('sem_g.R')
results<-sem_g(y, x, W, ndraw, nomit, prior)
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
prior<-prior_parse(NULL, k=2)
prior$novi_flag<-0
results<-sem_g(y, x, W, 150, 50, prior)
plot(results$pdraw, type="l", main="rho")
summary(results$pdraw)



##Homoc. model
prior$novi_flag<-1
results1<-sem_g(y, x, W, 150, 50, prior)
plot(results1$pdraw, type="l", main="rho")
summary(results1$pdraw)




#Example with Boston housing data
#Code from help(boston, package="spdep")
hr0 <- lm(log(MEDV) ~ CRIM + ZN + INDUS + CHAS + I(NOX^2) + I(RM^2) +
      AGE + log(DIS) + log(RAD) + TAX + PTRATIO + B + log(LSTAT)
, data=boston.c)

y<-matrix(log(boston.c$MEDV), ncol=1)
x<-model.matrix(hr0)
xlag<-cbind(x, W%*%x[,-1])

prior<-prior_parse(NULL, k=ncol(x))
prior$novi_flag<-0
prior$metflag<-1

resultsx<-sem_g(y, x, W, 3000, 1000, prior)
#plot(results$pdraw, type="l", main="rho")
plot(density(resultsx$pdraw), type="l", main="rho")
summary(resultsx$pdraw)
cbind(apply(resultsx$bdraw, 2, mean), apply(resultsx$bdraw, 2, sd))

#With lagged covariates
prior<-prior_parse(NULL, k=ncol(xlag))
prior$novi_flag<-0
prior$metflag<-1
resultsxlag<-sem_g(y, xlag, W, 3000, 1000, prior)
#plot(results$pdraw, type="l", main="rho")
plot(density(resultsxlag$pdraw), type="l", main="rho")
summary(results$pdraw)
apply(resultsxlag$bdraw, 2, mean)
apply(resultsxlag$bdraw, 2, mean)




