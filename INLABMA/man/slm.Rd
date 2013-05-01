\name{slm.inla}
\alias{slm.inla}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
slm.inla(formula, d, W, rho, mmatrix = NULL, improve = TRUE, impacts = FALSE, fhyper = NULL, probit = FALSE, ...)
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
  \item{improve}{
%%     ~~Describe \code{improve} here~~
}
  \item{impacts}{
%%     ~~Describe \code{impacts} here~~
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
function (formula, d, W, rho, mmatrix = NULL, improve = TRUE, 
    impacts = FALSE, fhyper = NULL, probit = FALSE, ...) 
{
    IrhoW <- diag(nrow(W)) - rho * W
    IrhoW2 <- IrhoW \%*\% t(IrhoW)
    if (is.null(mmatrix)) 
        mmatrix <- model.matrix(formula, d)
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
    res <- inla(fmla, data = d2, ...)
    if (improve) 
        res <- inla.hyperpar(res, diff.logdens = 20)
    res$logdet <- as.numeric(determinant(IrhoW2)$modulus)
    res$mlik <- res$mlik + res$logdet/2
    res$impacts <- FALSE
    if (impacts) {
        res$impacts <- TRUE
        if (!probit) {
            wtotal <- 1/(1 - rho)
            wdirect <- trIrhoWinv(W, rho)/nrow(W)
        }
        else {
            Df <- dnorm(res$summary.linear.predictor[, 1])
            wtotal <- mean(Df) * 1/(1 - rho)
            wdirect <- trIrhoWinv(W, rho, Df = diag(Df))/nrow(W)
        }
        windirect <- wtotal - wdirect
        res$summary.total.impacts <- wtotal * res$summary.fixed[rownames(res$summary.fixed) != 
            "(Intercept)", 1:2]
        res$summary.direct.impacts <- wdirect * res$summary.fixed[rownames(res$summary.fixed) != 
            "(Intercept)", 1:2]
        res$summary.indirect.impacts <- windirect * res$summary.fixed[rownames(res$summary.fixed) != 
            "(Intercept)", 1:2]
        res$marginals.total.impacts <- res$marginals.fixed[names(res$marginals.fixed) != 
            "(Intercept)"]
        res$marginals.direct.impacts <- res$marginals.fixed[names(res$marginals.fixed) != 
            "(Intercept)"]
        res$marginals.indirect.impacts <- res$marginals.fixed[names(res$marginals.fixed) != 
            "(Intercept)"]
        for (i in 1:length(res$marginals.total.impacts)) {
            xx <- res$marginals.total.impacts[[i]]
            res$marginals.total.impacts[[i]] <- rescalemarg(xx, 
                wtotal)
            res$marginals.direct.impacts[[i]] <- rescalemarg(xx, 
                wdirect)
            res$marginals.indirect.impacts[[i]] <- rescalemarg(xx, 
                windirect)
        }
    }
    return(res)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
