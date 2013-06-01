\name{fitmargBMA2}
\alias{fitmargBMA2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{fitmargBMA2
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
fitmargBMA2(models, ws, item)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{models}{
%%     ~~Describe \code{models} here~~
}
  \item{ws}{
%%     ~~Describe \code{ws} here~~
}
  \item{item}{
%%     ~~Describe \code{item} here~~
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
function (models, ws, item) 
{
    if (is.null(models[[1]][[item]])) 
        return(NULL)
    nlist <- names(models[[1]][[item]])
    auxlist <- as.list(rep(NA, length(nlist)))
    names(auxlist) <- nlist
    for (ele in nlist) {
        lmatrix <- lapply(models, function(X) {
            X[[item]][[ele]]
        })
        xxr <- c(NA, NA)
        xxr[1] <- sum(ws * unlist(lapply(lmatrix, function(X) {
            min(X[, 1])
        })))
        xxr[2] <- sum(ws * unlist(lapply(lmatrix, function(X) {
            max(X[, 1])
        })))
        xx <- seq(xxr[1], xxr[2], length.out = 81)
        auxbma <- rep(0, length(xx))
        for (i in 1:length(lmatrix)) {
            auxspl <- mysplinefun(x = lmatrix[[i]][, 1], y = lmatrix[[i]][, 
                2])
            auxbma <- auxbma + ws[i] * auxspl(xx)
        }
        auxlist[[ele]] <- cbind(xx, auxbma)
    }
    return(auxlist)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
