#
#Some tests
#
library(Matrix)

source("vgr_parse.r")
source("vgr_utils.R")

#Define adjacency matrix
W<-structure(c(0, .5, 0.5, 0.5, 0, 0.5, 0.5, 0.5, 0), .Dim = c(3L, 
3L))

#Test for lndetint
Wdetint<-lndetint(W)

lndettest<-sapply(Wdetint$rho, function(rho){log(det(diag(3)-rho*W))})

sum((Wdetint$lndet - lndettest)^2)
plot( Wdetint$lndet, lndettest); abline(0,1)

#Test for lndetfull
Wdetfull<-lndetfull(W)
lndetfulltest<-sapply(Wdetfull$rho, function(rho){log(det(diag(3)-rho*W))})

sum((Wdetfull$lndet - lndetfulltest)^2)
plot( Wdetfull$lndet, lndetfulltest); abline(0,1)

#Test for lndetmc
Wdetmc<-lndetmc(40, 20, W)
lndetmctest<-sapply(Wdetmc$rho, function(rho){log(det(diag(3)-rho*W))})

sum((Wdetmc$lndet - lndetmctest)^2)
#plot( Wdetmc$lndet, lndetmctest, type="l"); abline(0,1)
plot( Wdetmc$rho, lndetmctest, type="l");# abline(0,1)
lines(Wdetmc$rho, Wdetmc$lndet, lty=1, col="red")
lines(Wdetmc$rho, Wdetmc$lo95, lty=2);
lines(Wdetmc$rho, Wdetmc$hi95, lty=2);

#Test for far_g
source("vgr_far_g.r")

prior<-prior_parse(NULL)
prior$ldetflag<-0
prior$rmin<-0.1
prior$rmax<-0.99#If this is set to 1 then problem when interpolating log-det's

prior$metflag<-0#M-H sampler

y<-matrix(rnorm(3), ncol=1)
y<-solve(diag(3)-0.9*W)%*%y


#Heteroc.  model
prior$novi_flag<-0
results<-far_g(y, W, 1000, 0, prior)
plot(results$pdraw, type="l", main="rho")

#Homoc. model
prior$novi_flag<-1
results1<-far_g(y, W, 1000, 0, prior)
plot(results1$pdraw, type="l", main="rho")



#Bigger example
library(spdep)
data(boston)

W<-nb2mat(boston.soi)
y<-matrix(rnorm(nrow(W)), ncol=1)
y<-solve(diag(nrow(W))-0.4*W)%*%y


#Heteroc.  model
prior$novi_flag<-0
results<-far_g(y, W, 1000, 0, prior)
plot(results$pdraw, type="l", main="rho")
summary(results$pdraw)



#Homoc. model
prior$novi_flag<-1
results1<-far_g(y, W, 1000, 0, prior)
plot(results1$pdraw, type="l", main="rho")
summary(results1$pdraw)

