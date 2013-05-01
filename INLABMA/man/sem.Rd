\name{sem.inla}
\alias{sem.inla}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
sem.inla(formula, d, W, rho, improve = TRUE, impacts = FALSE, fhyper = NULL, probit = FALSE, ...)
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
function (formula, d, W, rho, improve = TRUE, impacts = FALSE, 
    fhyper = NULL, probit = FALSE, ...) 
{
    IrhoW <- diag(nrow(W)) - rho * W
    IrhoW2 <- IrhoW \%*\% t(IrhoW)
    environment(formula) <- environment()
    if (is.null(fhyper)) {
        formula <- update(formula, . ~ . + f(idx, model = "generic0", 
            Cmatrix = IrhoW2))
    }
    else {
        formula <- update(formula, . ~ . + f(idx, model = "generic0", 
            Cmatrix = IrhoW2, hyper = fhyper))
    }
    res <- inla(formula, data = d, ...)
    if (improve) 
        res <- inla.hyperpar(res, diff.logdens = 20)
    res$logdet <- as.numeric(determinant(IrhoW2)$modulus)
    res$mlik <- res$mlik + res$logdet/2
    res$impacts <- FALSE
    if (impacts) {
        res$impacts <- TRUE
        res$impacts$summary.impacts <- res$summary.fixed[rownames(res$summary.fixed) != 
            "(Intercept)", ]
        res$impacts$marginals.impacts <- res$marginals.fixed[names(res$marginals.fixed) != 
            "(Intercept)"]
        if (probit) {
            Dfm <- mean(dnorm(res$summary.linear.predictor[, 
                1]))
            mm <- lapply(res$impacts$marginals.impacts, function(X) {
                inla.marginal.transform(function(x) {
                  Dfm * x
                }, X)
            })
            names(mm) <- names(res$impacts$marginals.impacts)
            res$impacts$marginals.impacts <- mm
            res$impacts$summary.impacts <- Dfm * res$impacts$summary.impacts
        }
    }
    return(res)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
