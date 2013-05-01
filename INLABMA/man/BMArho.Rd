\name{BMArho}
\alias{BMArho}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
BMArho(models, rho, logrhoprior = rep(1, length(rho)))
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
function (models, rho, logrhoprior = rep(1, length(rho))) 
{
    mlik <- unlist(lapply(models, function(X) {
        X$mlik[1]
    }))
    post.func <- fitmarg(rho, mlik, logrhoprior)
    ws <- (post.func(rho))
    ws <- ws/sum(ws)
    fvalues <- mclapply(1:length(models), function(X) {
        ws[X] * models[[X]]$summary.fitted.values[, 1]
    })
    fvalues <- data.frame(fvalues)
    fvalues <- apply(fvalues, 1, sum)
    return(fvalues)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
