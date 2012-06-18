# Chebyshev approximation setup and run functions
get_trT <- function(W, q=5) {
# W a CSparseMatrix object
# q order
    n <- nrow(W)
    IM <- .symDiagonal(n)
    T <- vector(mode="list", length=(q+1))
    T[[1]] <- IM
    T[[2]] <- W
    trT <- numeric(q+1)
    trT[1] <- n
    trT[2] <- 0
    if (q > 1) {
        for (k in 3:(q+1)) {
            T[[k]] <- 2*(W %*% T[[(k-1)]]) - T[[(k-2)]]
            trT[k] <- sum(diag(T[[k]]))
        }
    }
    trT
}

lndetcheb_in <- function(trT, alpha) {
# trT output from get_trT()
# alpha spatial coefficient
    q <- length(trT)-1
    n <- trT[1]
    C1 <- cheb_in(alpha, j=1, q)
    x <- 0.0
    for (j in 1:(q+1)) {
        x <- x + (cheb_in(alpha, j=j, q)*trT[j])
    }
    x <- x - (n/2)*C1
    x
}
cheb_in <- function(alpha, j, q) {
    res <- (2/(q+1))
    x <- 0.0
    for (k in 1:(q+1)) {
        x <- x + log(((1 - (alpha*cos((pi*(k - 0.5))/(q + 1)))))) * 
            cos((pi*(j - 1)*(k - 0.5))/(q + 1))
    }
    res <- res * x
    res
}

# MC approximation setup and run functions
mcdet_prep <- function(W, p=16, m=30) {
# W a CSparseMatrix object
# p, m given in papers
	n <- dim(W)[1]
	x <- matrix(rnorm(n*p), nrow=n, ncol=p)
        int1  <- vector(mode="list", length=m)
	xx <- x
        for (k in 1:m) {
            xx <- W %*% xx
            int1[[k]] <- apply(x * as.matrix(xx), 2, sum)
        }
        int2 <- apply(x * x, 2, sum)
        list(m=m, p=p, n=n, int1=int1, int2=int2)
}

mcdet_in <- function(alpha, clx) {
# clx output from mcdet_prep()
# alpha spatial coefficient
	vk <- numeric(length=clx$p)
	for (k in 1:clx$m) {
		vk <- clx$n*(alpha^k)*(clx$int1[[k]]/k) + vk
	}
	v <- c(as.matrix(vk/clx$int2))
	v
}

