\name{trIrhoWinv}
\alias{trIrhoWinv}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
trIrhoWinv(W, rho, offset = 0, order = 20, direct = TRUE, Df = diag(nrow(W)))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{W}{
%%     ~~Describe \code{W} here~~
}
  \item{rho}{
%%     ~~Describe \code{rho} here~~
}
  \item{offset}{
%%     ~~Describe \code{offset} here~~
}
  \item{order}{
%%     ~~Describe \code{order} here~~
}
  \item{direct}{
%%     ~~Describe \code{direct} here~~
}
  \item{Df}{
%%     ~~Describe \code{Df} here~~
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
function (W, rho, offset = 0, order = 20, direct = TRUE, Df = diag(nrow(W))) 
{
    if (!direct) {
        lambdas <- eigen(W)$values
        tr <- 0
        for (i in 0:order) tr <- tr + sum(lambdas^(i + offset)) * 
            (rho^i)
    }
    else {
        tr <- 0
        WW <- diag(nrow(W))
        if (offset > 0) {
            for (i in 1:(offset)) WW <- WW \%*\% W
        }
        WW = Df \%*\% WW
        for (i in 0:order) {
            tr <- tr + sum(diag(WW)) * (rho^i)
            WW <- WW \%*\% W
        }
    }
    return(as.real(tr))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
