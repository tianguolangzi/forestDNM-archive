\name{pDNM_X_male}
\alias{pDNM_X_male}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
A function to create features from a male's X chromosome in a (human) VCF file.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
pDNM_X_male(vcf, pat.col = 1, mat.col = 2, child.col = 3)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{vcf}{
%%     ~~Describe \code{vcf} here~~
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
function (vcf, pat.col = 1, mat.col = 2, child.col = 3) 
{
    sn = seqnames(vcf)
    if (!all(sn \%in\% c("X", "chrX"))) 
        stop("this function intended for chrX only")
    st = start(vcf)
    g = geno(vcf)$GT
    p1 = (g[, pat.col] == g[, mat.col]) & g[, pat.col] != g[, 
        child.col]
    p2 = g[, pat.col] == "0/0" | g[, pat.col] == "1/1" | g[, 
        pat.col] == "2/2"
    c1 = g[, child.col] == "0/0" | g[, child.col] == "1/1" | 
        g[, child.col] == "2/2"
    ind = p1 & p2 & c1
    pa = as.integer(substr(g[, pat.col], 1, 1)) + 1
    ca = as.integer(substr(g[, child.col], 1, 1)) + 1
    mut = rep(NA, length(ca))
    mut[pa != ca] = ca[pa != ca]
    mut = mut[ind]
    pa = pa[ind]
    ca = ca[ind]
    sn = as.vector(sn[ind])
    st = st[ind]
    AD = geno(vcf)$AD[ind, ]
    PL = geno(vcf)$PL[ind, ]
    GQ = geno(vcf)$GQ[ind, ]
    GT = geno(vcf)$GT[ind, ]
    par_allele = rep("", length(REF))
    par_allele[pa == 1] = as.character(REF[pa == 1])
    par_allele[pa == 2] = sapply(ALT[pa == 2], function(x) as.character(x)[1])
    par_allele[pa == 3] = sapply(ALT[pa == 3], function(x) as.character(x)[2])
    mut_allele = rep("", length(REF))
    mut_allele[mut == 1] = as.character(REF[mut == 1])
    mut_allele[mut == 2] = sapply(ALT[mut == 2], function(x) as.character(x)[1])
    mut_allele[mut == 3] = sapply(ALT[mut == 3], function(x) as.character(x)[2])
    p_mut1 = mapply(function(x, y) x[y], AD[, pat.col], mut)
    p_mut2 = mapply(function(x, y) x[y], AD[, mat.col], mut)
    p_par1 = mapply(function(x, y) x[y], AD[, pat.col], pa) + 
        1
    p_par2 = mapply(function(x, y) x[y], AD[, mat.col], pa) + 
        1
    c_par = mapply(function(x, y) x[y], AD[, child.col], pa) + 
        1
    c_mut = mapply(function(x, y) x[y], AD[, child.col], mut)
    AD_X = data.frame(p_ar_max = pmax(p_mut1/p_par1, p_mut2/p_par2, 
        na.rm = T), p_ar_min = pmin(p_mut1/p_par1, p_mut2/p_par2, 
        na.rm = T), c_ar = c_mut/c_par)
    p_homo_PL1 = mapply(function(x, y) x[c(1, 3, 6)][y], PL[, 
        pat.col], pa)
    p_homo_PL2 = mapply(function(x, y) x[c(1, 3, 6)][y], PL[, 
        mat.col], pa)
    c_homo_PL = mapply(function(x, y) x[c(1, 3, 6)][y], PL[, 
        child.col], pa)
    c_gt = rep(NA, nrow(GT))
    c_gt[GT[, child.col] == "0/0"] = 1
    c_gt[GT[, child.col] == "1/1"] = 2
    c_gt[GT[, child.col] == "2/2"] = 3
    p_homoc_PL1 = mapply(function(x, y) x[c(1, 3, 6)][y], PL[, 
        pat.col], mut)
    p_homoc_PL2 = mapply(function(x, y) x[c(1, 3, 6)][y], PL[, 
        mat.col], mut)
    c_homoc_PL = mapply(function(x, y) x[c(1, 3, 6)][y], PL[, 
        child.col], mut)
    PL_X = data.frame(p_cg_max = pmax(p_homoc_PL1, p_homoc_PL2, 
        na.rm = T), p_cg_min = pmin(p_homoc_PL1, p_homoc_PL2, 
        na.rm = T), c_cg = c_homoc_PL, p_pg_max = pmax(p_homo_PL1, 
        p_homo_PL2, na.rm = T), p_pg_min = pmin(p_homo_PL1, p_homo_PL2, 
        na.rm = T), c_pg = c_homo_PL)
    X = data.frame(chr = sn, pos = st, QL = fixed(vcf)$QUAL[ind], 
        VQ = info(vcf)$VQSLOD[ind], SB = info(vcf)$SB[ind], BZ = info(vcf)$BaseQRankSum[ind], 
        DL = info(vcf)$Dels[ind], FS = info(vcf)$FS[ind], HS = info(vcf)$HaplotypeScore[ind], 
        MQ = info(vcf)$MQ[ind], MZ = info(vcf)$MQRankSum[ind], 
        QD = info(vcf)$QD[ind], PZ = info(vcf)$ReadPosRankSum[ind], 
        GQ_p_min = pmin(geno(vcf)$GQ[ind, pat.col], geno(vcf)$GQ[ind, 
            mat.col], na.rm = T), GQ_p_max = pmax(geno(vcf)$GQ[ind, 
            pat.col], geno(vcf)$GQ[ind, mat.col], na.rm = T), 
        GQ_c = geno(vcf)$GQ[ind, child.col], N_alt = sapply(alt(vcf)[ind], 
            length), par_allele = par_allele, mut_allele = mut_allele)
    X = cbind(AD_X, PL_X, X)
    rownames(X) = paste(sn, st, sep = ":")
    return(X)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
