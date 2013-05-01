\name{fitmargBMA}
\alias{fitmargBMA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
fitmargBMA(margs, ws, len = 100)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{margs}{
%%     ~~Describe \code{margs} here~~
}
  \item{ws}{
%%     ~~Describe \code{ws} here~~
}
  \item{len}{
%%     ~~Describe \code{len} here~~
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
function (margs, ws, len = 100) 
{
    ws <- ws/sum(ws)
    xmin <- max(unlist(lapply(margs, function(X) {
        min(X[, 1])
    })))
    xmax <- min(unlist(lapply(margs, function(X) {
        max(X[, 1])
    })))
    xx <- seq(xmin, xmax, len = len)
    margsws <- lapply(1:length(margs), function(i) {
        func <- fitmarg(margs[[i]][, 1], margs[[i]][, 2])
        ws[i] * func(xx)
    })
    margsws <- do.call(cbind, margsws)
    d <- data.frame(x = xx, y = apply(margsws, 1, sum))
    names(d) <- c("x", "y")
    return(d)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
