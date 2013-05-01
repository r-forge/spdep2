\name{BMA2}
\alias{BMA2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
BMA2(models, rho, logrhoprior = rep(1, length(rho)), impacts = FALSE, usenormal = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{models}{
%%     ~~Describe \code{models} here~~
}
  \item{rho}{
%%     ~~Describe \code{rho} here~~
}
  \item{logrhoprior}{
%%     ~~Describe \code{logrhoprior} here~~
}
  \item{impacts}{
%%     ~~Describe \code{impacts} here~~
}
  \item{usenormal}{
%%     ~~Describe \code{usenormal} here~~
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
function (models, rho, logrhoprior = rep(1, length(rho)), impacts = FALSE, 
    usenormal = FALSE) 
{
    mlik <- unlist(lapply(models, function(X) {
        X$mlik[1]
    }))
    post.func <- fitmarg(rho, mlik, logrhoprior, usenormal)
    ws <- (post.func(rho))
    ws <- ws/sum(ws)
    mfit <- list(rho = list())
    mfit$rho$marginal <- data.frame(x = seq(min(rho), max(rho), 
        len = 100))
    mfit$rho$marginal$y <- post.func(mfit$rho$marginal$x)
    mfit$rho$mean <- sum(mfit$rho$marginal$x * mfit$rho$marginal$y) * 
        diff(mfit$rho$marginal$x)[1]
    mfit$rho$sd <- sqrt(sum(((mfit$rho$marginal$x - mfit$rho$mean))^2 * 
        mfit$rho$marginal$y) * diff(mfit$rho$marginal$x)[1])
    fquant <- function(x, ff, qtile) {
        integrate(ff, 0, x)$value - qtile
    }
    qvect <- c(0.025, 0.5, 0.975)
    mfit$rho$quantiles <- sapply(qvect, function(qtile) {
        uniroot(fquant, interval = c(min(rho), max(rho)), ff = post.func, 
            qtile = qtile)$root
    })
    names(mfit$rho$quantiles) <- as.character(qvect)
    mateff <- c("summary.fixed", "summary.lincomb", "summary.linear.predictor", 
        "summary.fitted.values")
    lmat <- mclapply(mateff, function(X) {
        fitmatrixBMA(models, ws, X)
    })
    names(lmat) <- mateff
    mfit <- c(mfit, lmat)
    listeff <- c("dic", "cpo")
    leff <- mclapply(listeff, function(X) {
        fitlistBMA(models, ws, X)
    })
    names(leff) <- listeff
    mfit <- c(mfit, leff)
    listmarg <- c("marginals.fixed", "marginals.lincomb", "marginals.lincomb.derived", 
        "marginals.linear.predictor", "marginals.hyperpar", "marginals.spde2.blc")
    margeff <- mclapply(listmarg, function(X) {
        fitmargBMA2(models, ws, X)
    })
    names(margeff) <- listmarg
    mfit <- c(mfit, margeff)
    mfit$impacts <- FALSE
    if (impacts) {
        mfit$impacts <- TRUE
        summimp <- c("summary.total.impacts", "summary.direct.impacts", 
            "summary.indirect.impacts")
        matsummimp <- mclapply(summimp, function(X) {
            fitmatrixBMA(models, ws, X)
        })
        names(matsummimp) <- summimp
        mfit <- c(mfit, matsummimp)
        margimp <- c("marginals.total.impacts", "marginals.direct.impacts", 
            "marginals.indirect.impacts")
        lmargimp <- mclapply(margimp, function(X) {
            fitmargBMA2(models, ws, X)
        })
        names(lmargimp) <- margimp
        mfit <- c(mfit, lmargimp)
        mfit <- recompute.impacts(mfit)
    }
    return(mfit)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
