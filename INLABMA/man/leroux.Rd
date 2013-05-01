\name{leroux.inla}
\alias{leroux.inla}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
leroux.inla(formula, d, W, lambda, improve = TRUE, impacts = FALSE, fhyper = NULL, probit = FALSE, ...)
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
  \item{lambda}{
%%     ~~Describe \code{lambda} here~~
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
function (formula, d, W, lambda, improve = TRUE, impacts = FALSE, 
    fhyper = NULL, probit = FALSE, ...) 
{
    W2 <- diag(apply(W, 1, sum)) - W
    Q <- lambda * diag(nrow(W)) + (1 - lambda) * W2
    environment(formula) <- environment()
    if (is.null(fhyper)) {
        formula <- update(formula, . ~ . + f(idx, model = "generic0", 
            Cmatrix = Q))
    }
    else {
        formula <- update(formula, . ~ . + f(idx, model = "generic0", 
            Cmatrix = Q, hyper = fhyper))
    }
    res <- inla(formula, data = d, ...)
    if (improve) 
        res <- inla.hyperpar(res, diff.logdens = 20)
    res$logdet <- as.numeric(determinant(Q)$modulus)
    res$mlik <- res$mlik + res$logdet/2
    return(res)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
