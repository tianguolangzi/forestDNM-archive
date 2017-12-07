\name{predDNM}
\alias{predDNM}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
A function to predict de novo SNVs from features produced by dnmFeat().
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
predDNM(X, rf, cutoff = 0.2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
%%     ~~Describe \code{X} here~~
}
  \item{rf}{
%%     ~~Describe \code{rf} here~~
}
  \item{cutoff}{
%%     ~~Describe \code{cutoff} here~~
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
function (X, rf, cutoff = 0.2) 
{
    info = X[, colnames(X) \%in\% c("chr", "pos", "par_allele", 
        "mut_allele")]
    X = X[, !colnames(X) \%in\% c("chr", "pos", "par_allele", "mut_allele")]
    info$prd = predict(rf, X, type = "prob")[, 1]
    info = info[info$prd >= cutoff, ]
    return(info)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
