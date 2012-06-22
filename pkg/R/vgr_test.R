#
#Some tests
#

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
Wdetmc<-lndetmc(10, 10, W)
lndetmctest<-sapply(Wdetmc$rho, function(rho){log(det(diag(3)-rho*W))})

sum((Wdetmc$lndet - lndetmctest)^2)
plot( Wdetmc$lndet, lndetmctest); abline(0,1)

