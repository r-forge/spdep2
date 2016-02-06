spBreg_lag <- function(formula, data = list(), listw, na.action, type="lag",
    zero.policy=NULL, control=list()) {
    timings <- list()
    .ptime_start <- proc.time()
#control
    con <- list(tol.opt=.Machine$double.eps^0.5, ldet_method="SE_classic",
        Imult=2, cheb_q=5, MC_p=16L, MC_m=30L, super=NULL, spamPivot="MMD",
        in_coef=0.1, type="MC", correct=TRUE, trunc=TRUE,
        SE_method="LU", nrho=200, interpn=2000, SElndet=NULL, LU_order=FALSE,
        pre_eig=NULL, interval=c(-1, 1), ndraw=2500L, nomit=500L, thin=1L,
        verbose=FALSE, detval=NULL, prior=list(Tbeta=NULL, c_beta=NULL,
        rho=0.5, sige=1, nu=0, d0=0, a1 = 1.01, a2 = 1.01))
    nmsC <- names(con)
    con[(namc <- names(control))] <- control
    if (length(noNms <- namc[!namc %in% nmsC])) 
        warning("unknown names in control: ", paste(noNms, collapse = ", "))
    stopifnot(is.logical(con$verbose))
    stopifnot(is.integer(con$ndraw))
    stopifnot(is.integer(con$nomit))
    stopifnot(is.integer(con$thin))

    if (is.null(zero.policy))
        zero.policy <- get.ZeroPolicyOption()
    stopifnot(is.logical(zero.policy))
    if (class(formula) != "formula") formula <- as.formula(formula)
    mt <- terms(formula, data = data)
    mf <- lm(formula, data, na.action=na.action,  method="model.frame")
    na.act <- attr(mf, "na.action")
    if (!inherits(listw, "listw")) stop("No neighbourhood list")
    can.sim <- FALSE
    if (listw$style %in% c("W", "S")) can.sim <- spdep::can.be.simmed(listw)
    if (!is.null(na.act)) {
        subset <- !(1:length(listw$neighbours) %in% na.act)
        listw <- subset(listw, subset, zero.policy=zero.policy)
    }
    y <- model.extract(mf, "response")
    x <- model.matrix(mt, mf)
    n <- nrow(x)
    if (n != length(listw$neighbours))
        stop("Input data and weights have different dimensions")
    xcolnames <- colnames(x)
    K <- ifelse(xcolnames[1] == "(Intercept)", 2, 1)
    wy <- lag.listw(listw, y, zero.policy=zero.policy)
    if (anyNA(wy)) stop("NAs in lagged dependent variable")
#create_WX
    if (type == "Durbin") {
        WX <- create_WX(x, listw, zero.policy=zero.policy, prefix="lag")
        x <- cbind(x, WX)
        rm(WX)
    } else if (type != "lag") stop("No such type:", type)
    m <- ncol(x)
    lm.base <- lm(y ~ x - 1)
    aliased <- is.na(coefficients(lm.base))
    cn <- names(aliased)
    names(aliased) <- substr(cn, 2, nchar(cn))
    if (any(aliased)) {
        nacoef <- which(aliased)
        x <- x[,-nacoef]
    }
    m <- ncol(x)
    timings[["set_up"]] <- proc.time() - .ptime_start
    .ptime_start <- proc.time()

    env <- new.env()
    assign("can.sim", can.sim, envir=env)
    assign("listw", listw, envir=env)
    assign("similar", FALSE, envir=env)
    assign("n", n, envir=env)
    assign("verbose", con$verbose, envir=env)
    assign("family", "SAR", envir=env)
    W <- as(listw, "CsparseMatrix")
    assign("W", W, envir=env)

    interval <- spdep::jacobianSetup(con$ldet_method, env, con,
        pre_eig=con$pre_eig, interval=con$interval)
    bprior <-  dbeta(get("detval1", envir=env)[,1], con$prior$a1, con$prior$a2)
    assign("bprior", bprior, envir=env)

    nm <- paste(con$ldet_method, "set_up", sep="_")
    timings[[nm]] <- proc.time() - .ptime_start
    .ptime_start <- proc.time()

    k <- m 

    if (is.null(con$prior$c_beta))
        con$prior$c_beta <- rep(0, k)
    else 
        stopifnot(length(con$prior$c_beta) == k)

    if (is.null(con$prior$Tbeta))
        con$prior$Tbeta <- diag(k)*1e+12
    else
        stopifnot(nrow(con$prior$Tbeta) == k && ncol(con$prior$Tbeta) == k)

    sige <- con$prior$sige
    rho <- con$prior$rho

#% storage for draws
    bsave <- matrix(0, nrow=con$ndraw, ncol=k)
    psave <- numeric(con$ndraw)
    ssave <- numeric(con$ndraw)

#% ====== initializations
#% compute this stuff once to save time

    TI = solve(con$prior$Tbeta); # see eq 5.29, p. 140
    TIc = TI%*%con$prior$c_beta;
           
    iter = 1;
    xpx = crossprod(x)
    xpy = crossprod(x, y)
    Wy = W %*% y
    xpWy = crossprod(x, Wy)
    nu1 = n + 2*con$prior$nu 

    timings[["complete_setup"]] <- proc.time() - .ptime_start
    .ptime_start <- proc.time()


    while (iter <= con$ndraw) { #% start sampling;
                  
          ##% update beta   
        AI = solve((xpx + sige*TI))#,diag(rep(1,k)));        
        ys = y - rho*Wy          
        b = crossprod(x, ys) + sige*TIc
        b0 = solve((xpx + sige*TI), b)# see eq 5.29, p. 140
        bhat = MASS::mvrnorm(1, b0, sige*AI) #norm_rnd(sige*AI) + b0;  
        xb = x%*%bhat
        bsave[iter, 1:k] = as.vector(bhat)
          
          ##% update sige
        e = (ys - xb)
        d1 = 2*con$prior$d0 + crossprod(e)
        chi = rchisq(1, nu1) #chi = chis_rnd(1,nu1);
        sige = as.numeric(d1/chi) # see eq 5.30, p. 141
        ssave[iter] = as.vector(sige)
          
          ###% update rho using griddy Gibbs
        AI = solve((xpx + sige*TI))
        b0 = solve((xpx + sige*TI),(xpy + sige*TIc))
        bd = solve((xpx + sige*TI),(xpWy + sige*TIc))
        e0 = y - x%*%b0
        ed = Wy - x%*%bd
        epe0 = as.vector(crossprod(e0))
        eped = as.vector(crossprod(ed))
        epe0d = as.vector(crossprod(ed, e0))
        rho = draw_rho_env(env, epe0, eped, epe0d, n, k, rho)
        psave[iter] = as.vector(rho)

        iter = iter + 1

    }
### % end of sampling loop
    timings[["sampling"]] <- proc.time() - .ptime_start
    .ptime_start <- proc.time()

    mat <- cbind(bsave, psave, ssave)
    colnames(mat) <- c(colnames(x), "rho", "sige")
    res <- coda::mcmc(mat, start=con$nomit+1, end=con$ndraw, thin=con$thin)
    timings[["finalise"]] <- proc.time() - .ptime_start
    attr(res, "timings") <- do.call("rbind", timings)
    attr(res, "control") <- con
    res

#output mcmc object
}
