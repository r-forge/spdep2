
library(condMVNorm)
library(MCMCpack)
library(MASS)
library(spdep)
data(boston)

fm <- CMEDV~CRIM + AGE + RM
x <- cbind(1, boston.c$CRIM, boston.c$AGE,  boston.c$RM)
y <- boston.c$CMEDV
w <- listw2mat(nb2listw(boston.soi))
w <- Matrix(w)
Im <- Diagonal(ncol(w))
stsres <- stsls(fm, data = boston.c,listw = nb2listw(boston.soi))

lmbd <- coefficients(stsres)[1]
selmbd <- stsres$var[1,1]
stsres$sse/n
summary(stsres)
Mllag <- lagsarlm(fm, data = boston.c,listw = nb2listw(boston.soi))
lllmbd <- Mllag$rho
vl <- vcov(Mllag)[1,1]
# Generating the initial values 
# What is the best way of generating 

# Initial values for lambda, s2, and b
lambda_0 <- 0.5
s2_0 <- 0.9
b <- rep(0.5, 4) 

xpx <- crossprod(x)
vcx <- solve(xpx)

niter <- 1000
burn <- 150
mcmc <- niter - burn
B1 <- Matrix(0, mcmc, 4)
S2 <- Matrix(0, mcmc, 1)
Lambda <- Matrix(0, mcmc, 1)

i <- 1
while(i < niter){

#update beta
# using initial values for lambda and s2 
# compute bhat and  
yst <- y - lambda_0 * w %*% y
dtr <- det(Im - lambda_0 *w)
bhat <- vcx %*% crossprod(x,yst)
ehat <- yst - x %*% bhat

A <- crossprod(ehat) + t((bhat - b)) %*% xpx %*% (bhat - b)
b <- mvrnorm( n=1, mu = bhat, Sigma = (s2_0 * vcx)) 

if(i > burn) B1[(i-burn),] <- b

# update s2_0 
gamma <-  as.numeric(A /2)
alpha <- ((n+1)/2)-1
s2_1 <- rinvgamma(1, alpha, gamma)
# fl <- log(dtr) * (-1/(2*s2_0)) * A

if(i > burn) S2[(i-burn),] <- s2_1

# lambda_1 <- runif(1, 0, 2.5)
lambda_1 <- rnorm(1, lllmbd, sqrt(vl))
# print(lambda_1)
# print(lambda_1)
yst1 <- y - lambda_1*w %*% y
dtr1 <- det(Im - lambda_1 *w)
bhat1 <- vcx %*% crossprod(x,yst1)
ehat1 <- yst1 - x %*% bhat1
A1 <- crossprod(ehat1) + t((bhat1 - b)) %*% xpx %*% (bhat1 - b)
fl <- log(dtr) * (-1/(2*s2_1)) * A1
fl1 <- log(dtr1) * (-1/(2*s2_1)) * A1

# print(fl1/fl)
ara <- c(as.numeric(exp(fl1/fl)), 1)
alpha <- min(ara)
# print(alpha)

if(i > burn){
	
if(alpha == 1) Lambda[(i-burn),] <- lambda_1

else{
	
ru <- runif(1, 0, 1)	

if(ru > alpha) 	Lambda[(i-burn),] <- lambda_0
else Lambda[(i-burn),] <- lambda_1
}
}

lambda_0 <- lambda_1
s2_0 <- s2_1

i = i+1	
}

colSums(Lambda)/ mcmc
colSums(B1)/ mcmc
colSums(S2)/ mcmc

plot(density(as.numeric(Lambda)))


