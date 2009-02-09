ml_env_setup <- function(formula, data, listw, verbose=TRUE) {
    if (!is.logical(verbose)) verbose=FALSE
    mt <- terms(formula, data=data)
    mf <- lm(formula=formula, data = data, method = "model.frame")
    y <- model.response(mf, "numeric")
    n <- length(y)
    x <- model.matrix(mt, mf)
    W <- as(as_dgRMatrix_listw(listw), "CsparseMatrix")
    wy <- as.matrix(W %*% y)
    WX <- as.matrix(W %*% x)
    env <- new.env(parent=globalenv())
    assign("n", n, envir=env)
    assign("verbose", verbose, envir=env)
    assign("y", y, envir=env)
    assign("x", x, envir=env)
    assign("wy", wy, envir=env)
    assign("WX", WX, envir=env)
    I <- .symDiagonal(n)
    slw <- similar.listw(listw)
    dsT <- as_dsTMatrix_listw(slw)
    dWd <- as(dsT, "dsCMatrix")
    ndWd <- -dWd
    C1p <- Cholesky(ndWd, Imult = 2)
    C1n <- Cholesky(dWd, Imult = 2)
    a <- 0 - .Machine$double.eps
    b <- 0 + .Machine$double.eps
    assign("a", a, envir=env)
    assign("b", b, envir=env)
    assign("dWd", dWd, envir=env)
    assign("ndWd", ndWd, envir=env)
    assign("C1p", C1p, envir=env)
    assign("C1n", C1n, envir=env)

    env
}

do_LL <- function(val, env, interp=c(FALSE, FALSE)) {
    n <- get("n", envir=env)
    verbose <- get("verbose", envir=env)
    if (interp[1]) SSE <- sse_fn(env, val)
    else SSE <- .Call("R_ml_sse_env", env, val, PACKAGE="spdep2")
    s2 <- SSE/n
    if (interp[2]) Jacobian <- J_fn(env, val)
    else Jacobian <- .Call("R_ml_Jac_env", env, val, PACKAGE="spdep2")
    loglik <- Jacobian - ((n/2) * log(2 * pi)) -
        (n/2) * log(s2) - (1/(2 * (s2))) * SSE
    if (verbose) cat("val:", val, "SSE:", SSE, "Jacobian:", Jacobian,
    "LL:", loglik, "\n")
    loglik
}

sse_fn <- function(env, lambda) {
    yl <- get("y", envir=env) - lambda * get("wy", envir=env)
    xl <- get("x", envir=env) - lambda * get("WX", envir=env)
    xlQR <- qr(xl)
    xl.q <- qr.Q(xlQR)
    xl.q.yl <- crossprod(xl.q, yl)
    SSE <- crossprod(yl) - crossprod(xl.q.yl)
    SSE
}


J_fn <- function(env, lambda) {
    n <- get("n", envir=env)
    a <- get("a", envir=env)
    b <- get("b", envir=env)
    dWd <- get("dWd", envir=env)
    ndWd <- get("ndWd", envir=env)
    C1p <- get("C1p", envir=env)
    C1n <- get("C1n", envir=env)
    Jacobian <- ifelse(lambda > b, n * log(lambda) +
        c(determinant(update(C1p, ndWd, 1/lambda))$modulus),
        ifelse(lambda < a, n* log(-(lambda)) + c(determinant(update(C1n,
        dWd, 1/(-lambda)))$modulus), 0.0))
    Jacobian
}

