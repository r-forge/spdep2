\name{recompute.impacts}
\alias{recompute.impacts}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
recompute.impacts(obj, impacts = c("total", "direct", "indirect"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obj}{
%%     ~~Describe \code{obj} here~~
}
  \item{impacts}{
%%     ~~Describe \code{impacts} here~~
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
function (obj, impacts = c("total", "direct", "indirect")) 
{
    if (!obj$impacts) 
        return(obj)
    for (ii in impacts) {
        sumname <- paste("summary.", ii, ".impacts", sep = "")
        margname <- paste("marginals.", ii, ".impacts", sep = "")
        stab <- obj[[sumname]]
        for (i in 1:nrow(stab)) {
            xmean <- inla.emarginal(function(x) {
                x
            }, obj[[margname]][[i]])
            xsd <- sqrt(inla.emarginal(function(x) {
                (x - xmean)^2
            }, obj[[margname]][[i]]))
            stab[i, 1:2] <- c(xmean, xsd)
        }
        obj[[sumname]] <- stab[, 1:2]
    }
    return(obj)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
