ml_env_setup <- function(formula, data, listw, weights,
    na.action=na.fail, verbose=TRUE, similar=TRUE, imult=2, CAR=FALSE) {
    if (!is.logical(verbose)) verbose=FALSE
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "weights", "na.action"), names(mf), 0)
    mf <- mf[c(1, m)]
    mf$drop.unused.levels <- TRUE
    mf[[1]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame())
    mt <- attr(mf, "terms")
    na.act <- attr(mf, "na.action")
    if (!is.null(na.act)) {
        subset <- !(1:length(listw$neighbours) %in% na.act)
        listw <- subset(listw, subset, zero.policy=zero.policy)
    }

    dcar <- ifelse(CAR, 0.5, 1.0)
    
    y <- model.extract(mf, "response")
    if (any(is.na(y))) stop("NAs in dependent variable")
    n <- length(y)
    x <- model.matrix(mt, mf)
    if (any(is.na(x))) stop("NAs in independent variable")
    weights <- as.vector(model.extract(mf, "weights"))
    if (!is.null(weights) && !is.numeric(weights)) 
        stop("'weights' must be a numeric vector")
    if (is.null(weights)) weights <- rep(as.numeric(1), n)
    if (any(is.na(weights))) stop("NAs in weights")
    if (any(weights < 0)) stop("negative weights")
    sw <- sqrt(weights)
    sum_lw <- sum(log(weights))

    W <- as(as_dgRMatrix_listw(listw), "CsparseMatrix")
    wy <- as.matrix(W %*% y) * sw
    y <- y * sw
    WX <- as.matrix(W %*% x) * sw
    x <- x * sw
    env <- new.env(parent=globalenv())
    assign("n", n, envir=env)
    assign("sum_lw", sum_lw, envir=env)
    assign("dcar", dcar, envir=env)
    assign("verbose", verbose, envir=env)
    assign("y", y, envir=env)
    assign("x", x, envir=env)
    assign("wy", wy, envir=env)
    assign("WX", WX, envir=env)
    I <- .symDiagonal(n)
    if (similar) listw <- similar.listw(listw)
    dsT <- as_dsTMatrix_listw(listw)
    dWd <- as(dsT, "dsCMatrix")
    ndWd <- -dWd
    C1p <- Cholesky(ndWd, Imult = imult)
    C1n <- Cholesky(dWd, Imult = imult)
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
    sum_lw <- get("sum_lw", envir=env)
    dcar <- get("dcar", envir=env)
    verbose <- get("verbose", envir=env)
    if (dcar < 1.0) {
      SSE <- sse_fn_car(env, val)
    } else {
      if (interp[1]) SSE <- sse_fn(env, val)
      else SSE <- .Call("R_ml_sse_env", env, val, PACKAGE="spdep2")
    }
    s2 <- SSE/n
    if (interp[2]) Jacobian <- J_fn(env, val)
    else Jacobian <- .Call("R_ml_Jac_env", env, val, PACKAGE="spdep2")
    loglik <- Jacobian*dcar + (1/2)*sum_lw - ((n/2) * log(2 * pi)) -
        ((n/2) * log(s2)) - (1/(2 * (s2))) * SSE
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

sse_fn_car <- function(env, lambda) {
    n <- get("n", envir=env)
    a <- get("a", envir=env)
    b <- get("b", envir=env)
    y <- get("y", envir=env)
    x <- get("x", envir=env)
    if (lambda > b) Z <- as(update(get("C1p", envir=env), get("ndWd", envir=env), 1/lambda), "sparseMatrix")
    else (lambda < a) Z <- as(update(get("C1n", envir=env), get("dWd", envir=env), 1/(-lambda)), "sparseMatrix")
    else Z <- .symDiagonal(n)
#    Z <- ifelse(lambda > b, update(get("C1p", envir=env), get("ndWd", envir=env), 1/lambda), ifelse(lambda < a, update(get("C1n", envir=env), get("dWd", envir=env), 1/(-lambda)), .symDiagonal(n)))
    yl <- Z %*% y
    xl <- Z %*% x
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

