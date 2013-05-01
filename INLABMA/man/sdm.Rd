\name{sdm.inla}
\alias{sdm.inla}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
sdm.inla(formula, d, W, rho, mmatrix = NULL, intercept = TRUE, impacts = FALSE, improve = TRUE, fhyper = NULL, probit = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{
%%     ~~Describe \code{formula} here~~
}
  \item{d}{
%%     ~~Describe \code{d} here~~
}
  \item{W}{
%%     ~~Describe \code{W} here~~
}
  \item{rho}{
%%     ~~Describe \code{rho} here~~
}
  \item{mmatrix}{
%%     ~~Describe \code{mmatrix} here~~
}
  \item{intercept}{
%%     ~~Describe \code{intercept} here~~
}
  \item{impacts}{
%%     ~~Describe \code{impacts} here~~
}
  \item{improve}{
%%     ~~Describe \code{improve} here~~
}
  \item{fhyper}{
%%     ~~Describe \code{fhyper} here~~
}
  \item{probit}{
%%     ~~Describe \code{probit} here~~
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (formula, d, W, rho, mmatrix = NULL, intercept = TRUE, 
    impacts = FALSE, improve = TRUE, fhyper = NULL, probit = FALSE, 
    ...) 
{
    IrhoW <- diag(nrow(W)) - rho * W
    IrhoW2 <- IrhoW \%*\% t(IrhoW)
    if (is.null(mmatrix)) {
        mmatrix <- model.matrix(formula, d)
        if (intercept) 
            mmW <- W \%*\% mmatrix[, -1]
        else mmW <- W \%*\% mmatrix
        mmatrix <- cbind(mmatrix, as.matrix(mmW))
    }
    mm <- as.data.frame(as.matrix(solve(IrhoW) \%*\% mmatrix))
    names(mm) <- paste("x", 1:ncol(mm), sep = "")
    xnam <- names(mm)
    d2 <- cbind(d, mm)
    fmla <- paste(as.character(formula)[2], "~ -1+", paste(xnam, 
        collapse = "+"))
    if (is.null(fhyper)) 
        fmla <- paste(fmla, "+f(idx, model=\"generic0\", Cmatrix=IrhoW2)", 
            sep = "")
    else fmla <- paste(fmla, "+f(idx, model=\"generic0\", Cmatrix=IrhoW2, hyper=fhyper)", 
        sep = "")
    fmla <- as.formula(fmla)
    if (impacts) {
        if (!probit) {
            wtotal <- rep(1/(1 - rho), 2)
            wdirect <- c(trIrhoWinv(W, rho), trIrhoWinv(W, rho, 
                1))/nrow(W)
        }
        else {
            error("Cannot compute impacts because we need to run the model twice")
            Df <- dnorm(res$summary.linear.predictor[, 1])
            wtotal <- mean(Df) * rep(1/(1 - rho), 2)
            wdirect <- c(trIrhoWinv(W, rho, Df = diag(Df)), trIrhoWinv(W, 
                rho, 1, Df = diag(Df)))/nrow(W)
        }
        windirect <- wtotal - wdirect
        ncov <- (ncol(mm) - 1)/2
        lc.total.impacts <- lapply(1:ncov, function(X) {
            lc <- list(a = wtotal[1], b = wtotal[2])
            names(lc) <- xnam[1 + c(X, X + ncov)]
            inla.make.lincomb(lc)
        })
        lc.total.impacts <- do.call(c, lc.total.impacts)
        names(lc.total.impacts) <- paste("totalimp.", xnam[1 + 
            1:ncov], sep = "")
        lc.direct.impacts <- lapply(1:ncov, function(X) {
            lc <- list(a = wdirect[1], b = wdirect[2])
            names(lc) <- xnam[1 + c(X, X + ncov)]
            inla.make.lincomb(lc)
        })
        lc.direct.impacts <- do.call(c, lc.direct.impacts)
        names(lc.direct.impacts) <- paste("directimp.", xnam[1 + 
            1:ncov], sep = "")
        lc.indirect.impacts <- lapply(1:ncov, function(X) {
            lc <- list(a = windirect[1], b = windirect[2])
            names(lc) <- xnam[1 + c(X, X + ncov)]
            inla.make.lincomb(lc)
        })
        lc.indirect.impacts <- do.call(c, lc.indirect.impacts)
        names(lc.indirect.impacts) <- paste("indirectimp.", xnam[1 + 
            1:ncov], sep = "")
        lc.impacts <- c(lc.total.impacts, lc.direct.impacts, 
            lc.indirect.impacts)
        res <- inla(fmla, data = d2, lincomb = lc.impacts, ...)
    }
    else {
        res <- inla(fmla, data = d2, ...)
    }
    if (improve) 
        res <- inla.hyperpar(res, diff.logdens = 20)
    res$logdet <- as.numeric(determinant(IrhoW2)$modulus)
    res$mlik <- res$mlik + res$logdet/2
    res$impacts <- FALSE
    if (impacts) {
        res$impacts <- TRUE
        idx <- substr(rownames(res$summary.lincomb.derived), 
            1, 3)
        res$summary.total.impacts <- res$summary.lincomb.derived[idx == 
            "tot", ]
        res$summary.direct.impacts <- res$summary.lincomb.derived[idx == 
            "dir", ]
        res$summary.indirect.impacts <- res$summary.lincomb.derived[idx == 
            "ind", ]
        idx <- substr(names(res$marginals.lincomb.derived), 1, 
            3)
        res$marginals.total.impacts <- res$marginals.lincomb.derived[idx == 
            "tot"]
        res$marginals.direct.impacts <- res$marginals.lincomb.derived[idx == 
            "dir"]
        res$marginals.indirect.impacts <- res$marginals.lincomb.derived[idx == 
            "ind"]
    }
    return(res)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
