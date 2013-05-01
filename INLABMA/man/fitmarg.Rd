\name{fitmarg}
\alias{fitmarg}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
fitmarg(x, logy, logp = 0, usenormal = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
}
  \item{logy}{
%%     ~~Describe \code{logy} here~~
}
  \item{logp}{
%%     ~~Describe \code{logp} here~~
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
function (x, logy, logp = 0, usenormal = FALSE) 
{
    if (!usenormal) {
        func = splinefun(x, logy - max(logy))
        post = exp(logp + func(x))
        post.func = splinefun(x, post)
        xx = seq(min(x), max(x), len = 1000)
        z = integrate(post.func, min(x), max(x))$value
        post.func = mysplinefun(x = xx, y = (post.func(xx)/z))
    }
    else {
        xx = seq(min(x), max(x), len = 1000)
        meany = sum(x * exp(logp + logy))/sum(exp(logp + logy))
        sdy = sqrt(sum(((x - meany)^2) * exp(logp + logy))/sum(exp(logp + 
            logy)))
        post.func = function(x) {
            dnorm(x, mean = meany, sd = sdy)
        }
    }
    return(post.func)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
