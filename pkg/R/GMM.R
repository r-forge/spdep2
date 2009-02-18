# vv <- .kpwuwu(listw, residuals(ols), zero.policy=zero.policy)
# optres <- nlminb(pars, .kpgm, v=vv, verbose=verbose, control=control)
#
# Copyright 2004 by Luc Anselin
# Kelejian-Prucha generalized moments equations
# helper function to provide function to nonlinear optimizer
# must have parameter vector first for nlm
# Usage:
#    kpgm(par,v)
# Arguments:
#    par: 2x1 parameter vector rho,sig2
#    v: list containing bigG and litg as computed by kpwuwu
# Details:
#    sets up the equation as squared residuals
# Value:
#    value: evaluated nonlinear least squares for parameter value

.kpgm <- function(rhopar,v,verbose=FALSE) {
  vv <- v$bigG %*% c(rhopar[1],rhopar[1]^2,rhopar[2]) - v$litg
  value <- sum(vv^2)
  if (verbose)
    cat("function:", value, "lambda:", rhopar[1], "sig2:", rhopar[2], "\n")
  value
  
}

# Kelejian-Prucha generalized moments equations 2005-2008

.kpgm2 <- function(rho, v, verbose=FALSE) {
  vv <- v$litg - v$bigG %*% c(rho,rho^2)
  value <- c(t(vv) %*% v$Ups %*% vv)
  if (verbose)
    cat("function:", value, "lambda:", rho, "\n")
  value
}

.kpwuwu2 <- function(W, u, UPS=1) {
  n <- length(u)
  ubar <- W %*% u
  ubarbar <- W %*% ubar
  A1 <- (t(W) %*% W) - .symDiagonal(n, x=apply(W, 2, crossprod))
  litG <- c(as.matrix(crossprod(u, (A1 %*% u))),
    as.matrix(crossprod(u, ubar)))/n
  bigG <- matrix(0,2,2)
  bigG[1,1] <- as.matrix(crossprod(ubar, (A1 %*% u)))/2*n
  bigG[1,2] <- as.matrix(crossprod(ubar, (A1 %*% ubar)))/n
  bigG[2,1] <- as.matrix(crossprod(ubar, ((W + t(W)) %*% u)))/n
  bigG[2,2] <- as.matrix(crossprod(ubar, ubarbar))/n
  Ups <- diag(2)
  if (UPS == 2) Ups <- (1/(1 + ((sum(diag(W %*% t(W))))/n)^2)) * Ups
  list(litG=litG, bigG=bigG, Ups=Ups)
}

# Copyright 2004 by Luc Anselin
# Kelejian-Prucha generalized moments equations
# helper function
# Usage:
#    kpwuwu(listw,u)
# Arguments:
#    listw: spatial weights file as listw object
#    u: OLS residual vector
#    zero.policy: allow no-neighbour observations if TRUE
# Details:
#    sets up the bigG matrix and littleg vector needed
#    for the nonlinear least squares in the GM estimator
#    see Kelejian-Prucha(1999) p. 515
# Value:
# a list with two elements
#    bigG: the 3x3 G matrix
#    litg: the 3x1 g vector

.kpwuwu <- function(W, u, zero.policy=FALSE) {
	n <- length(u)
# Gianfranco Piras 081119 
        trwpw <- sum(unlist(W$weights)^2)
#	tt <- matrix(0,n,1)
#	for (i in 1:n) {tt[i] <- sum(W$weights[[i]]^2) }
#	trwpw <- sum(tt)
	wu <- lag.listw(W, u, zero.policy=zero.policy)
	wwu <- lag.listw(W, wu, zero.policy=zero.policy)
    	uu <- crossprod(u,u)
    	uwu <- crossprod(u,wu)
    	uwpuw <- crossprod(wu,wu)
    	uwwu <- crossprod(u,wwu)
    	wwupwu <- crossprod(wwu,wu)
    	wwupwwu <- crossprod(wwu,wwu)
    	bigG <- matrix(0,3,3)
    	bigG[,1] <- c(2*uwu,2*wwupwu,(uwwu+uwpuw))/n
    	bigG[,2] <- - c(uwpuw,wwupwwu,wwupwu) / n
    	bigG[,3] <- c(1,trwpw/n,0)
    	litg <- c(uu,uwpuw,uwu) / n
    	list(bigG=bigG,litg=litg)
}

