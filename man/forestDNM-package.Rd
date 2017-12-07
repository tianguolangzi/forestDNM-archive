\name{forestDNM-package}
\alias{forestDNM-package}
\alias{forestDNM}
\docType{package}
\title{
Predicts de novo SNVs when provided a VCF containing variants of a genotyped family.
}
\description{
forestDNM is an R package built around a classifier that was trained to predict true de novo germline mutations (DNMs), using features derived from family genotype data contained in a VCF.  The classifier was trained on 10 families with monozygotic twins, whose putative DNMs had undergone extensive experimental validation (the classifier was trained to predict validation status).  In an independent test set of held-out data from the 10 families, sensitivity was > 95% while maintaining an FDR < 10%.  The balance between sensitivity and FDR can be tuned by adjusting the threshold used on the classifier score (RF vote proportion).  The default cutoff is 0.2, and plots showing how FDR and sensitivity vary with classifier score are included in the vignette.  SNVs in this training/test set were genotyped using GATK 2.1-13 (unified genotyper and VQSR).   
}
\details{
\tabular{ll}{
Package: \tab forestDNM\cr
Type: \tab Package\cr
Version: \tab 0.1\cr
Date: \tab 2013-05-14\cr
License: \tab What license is it under?\cr
}

}
\author{
Jacob J. Michaelson

Maintainer: Jake Michaelson <jake.michaelson@gmail.com>

}
\references{
Michaelson, et al. Whole-genome sequencing in autism identifies hot spots for de novo germline mutation. Cell (2012).
}


\keyword{ package }
\seealso{
randomForest
}
\examples{

}
