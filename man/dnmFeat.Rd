\name{dnmFeat}
\alias{dnmFeat}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
A wrapper to create features from all chromosomes in a (human) VCF file.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
dnmFeat(filename, child.sex, pat.col = 1, mat.col = 2, child.col = 3, genome = "hg18", chrom.conv = "chr")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{filename}{
%%     ~~Describe \code{filename} here~~
}
  \item{child.sex}{
%%     ~~Describe \code{child.sex} here~~
}
  \item{pat.col}{
%%     ~~Describe \code{pat.col} here~~
}
  \item{mat.col}{
%%     ~~Describe \code{mat.col} here~~
}
  \item{child.col}{
%%     ~~Describe \code{child.col} here~~
}
  \item{genome}{
%%     ~~Describe \code{genome} here~~
}
  \item{chrom.conv}{
%%     ~~Describe \code{chrom.conv} here~~
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
function (filename, child.sex, pat.col = 1, mat.col = 2, child.col = 3, 
    genome = "hg18", chrom.conv = "chr") 
{
    require(VariantAnnotation)
    if (missing(child.sex) || !child.sex \%in\% c("M", "F")) 
        stop("offspring's sex must be specified (M or F)")
    chr = paste(chrom.conv, c(1:22), sep = "")
    X = list()
    for (i in 1:length(chr)) {
        rg = RangedData(IRanges(start = 1, end = 3e+08), space = chr[i])
        param <- ScanVcfParam(which = rg, geno = c("AD", "DP", 
            "GQ", "GT", "PL"), info = c("VQSLOD", "SB", "DB", 
            "BaseQRankSum", "Dels", "FS", "HaplotypeScore", "MQ", 
            "MQRankSum", "QD", "ReadPosRankSum", "AN"))
        vcf = readVcf(filename, genome, param)
        X[[i]] = pDNM(vcf, pat.col, mat.col, child.col)
        print(paste("done with", chr[i]))
    }
    if (child.sex == "M") {
        rg = RangedData(IRanges(start = 1, end = 3e+08), space = paste(chrom.conv, 
            "X", sep = ""))
        param <- ScanVcfParam(which = rg, geno = c("AD", "DP", 
            "GQ", "GT", "PL"), info = c("VQSLOD", "SB", "DB", 
            "BaseQRankSum", "Dels", "FS", "HaplotypeScore", "MQ", 
            "MQRankSum", "QD", "ReadPosRankSum", "AN"))
        vcf = readVcf(filename, genome, param)
        X[[23]] = pDNM_X_male(vcf, pat.col, mat.col, child.col)
        print("done with chrX")
        rg = RangedData(IRanges(start = 1, end = 3e+08), space = paste(chrom.conv, 
            "Y", sep = ""))
        param <- ScanVcfParam(which = rg, geno = c("AD", "DP", 
            "GQ", "GT", "PL"), info = c("VQSLOD", "SB", "DB", 
            "BaseQRankSum", "Dels", "FS", "HaplotypeScore", "MQ", 
            "MQRankSum", "QD", "ReadPosRankSum", "AN"))
        vcf = readVcf(filename, genome, param)
        X[[24]] = pDNM_Y(vcf, pat.col, child.col)
        X = do.call("rbind", X)
        print("done with chrY")
    }
    else {
        rg = RangedData(IRanges(start = 1, end = 3e+08), space = paste(chrom.conv, 
            "X", sep = ""))
        param <- ScanVcfParam(which = rg, geno = c("AD", "DP", 
            "GQ", "GT", "PL"), info = c("VQSLOD", "SB", "DB", 
            "BaseQRankSum", "Dels", "FS", "HaplotypeScore", "MQ", 
            "MQRankSum", "QD", "ReadPosRankSum", "AN"))
        vcf = readVcf(filename, genome, param)
        X[[23]] = pDNM_X_female(vcf, pat.col, mat.col, child.col)
        X = do.call("rbind", X)
        print("done with chrX")
    }
    X = X[!duplicated(paste(X$chr, X$pos, sep = ":")), ]
    return(X)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
