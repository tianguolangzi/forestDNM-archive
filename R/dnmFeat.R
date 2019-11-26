dnmFeat <-
function(filename,child.sex,pat.col=1,mat.col=2,child.col=3,
  genome="hg18",chrom.conv="chr"){
 require(VariantAnnotation,quietly=TRUE)
 if(missing(child.sex) || !child.sex%in%c("M","F")) stop("offspring's sex must be specified (M or F)")

 ### set up progress reporter
 time.st = proc.time()[3]
 prop.done = log(((1:23)*0.07303)+1.06785) ## rough approximation
 prop.done[1] = 0.08
 prop.done[23] = 1

 ### get the autosomes
 chr = paste(chrom.conv,c(1:22),sep="")
 X = list()

 for(i in 1:length(chr)){
  rg = GRanges(seqnames=chr[i],IRanges(start=1,end=3e8))
  ### set params
  param <- ScanVcfParam(which=rg,geno=c("AD","DP","GQ","GT","PL"),
  info = c("VQSLOD","DB","BaseQRankSum","Dels","FS",
    "HaplotypeScore","MQ","MQRankSum",
    "QD","ReadPosRankSum","AN"))
  ### read in file
  vcf = readVcf(filename,genome,param)
  X[[i]] = pDNM(vcf,pat.col,mat.col,child.col)
  ## report progress
  elapsed = proc.time()[3] - time.st
  tot.est = elapsed/prop.done[i]
  remaining = round(tot.est-elapsed,0)
  cat(paste("estimated completion in",round(remaining/60,1),"minutes [",
  format(Sys.time()+remaining, "%a %b %d %X"),"] \n"))
 }

 ### do the sex chromosomes
 if(child.sex=="M"){
  rg = GRanges(seqnames=paste(chrom.conv,"X",sep=""),IRanges(start=1,end=3e8))
  ### set params
  param <- ScanVcfParam(which=rg,geno=c("AD","DP","GQ","GT","PL"),
  info = c("VQSLOD","DB","BaseQRankSum","Dels","FS",
    "HaplotypeScore","MQ","MQRankSum",
    "QD","ReadPosRankSum","AN"))
  ### read in file
  vcf = readVcf(filename,genome,param)
  X[[23]] = pDNM_X_male(vcf,mat.col,child.col)
  ## report progress
  elapsed = proc.time()[3] - time.st
  tot.est = elapsed/prop.done[23]
  remaining = round(tot.est-elapsed,0)
  cat(paste("estimated completion in",round(remaining/60,1),"minutes [",
  format(Sys.time()+remaining, "%a %b %d %X"),"] \n"))

  rg = GRanges(seqnames=paste(chrom.conv,"Y",sep=""),IRanges(start=1,end=3e8),)
  ### set params
  cat("working on Y chromosome...")
  param <- ScanVcfParam(which=rg,geno=c("AD","DP","GQ","GT","PL"),
  info = c("VQSLOD","DB","BaseQRankSum","Dels","FS",
    "HaplotypeScore","MQ","MQRankSum",
    "QD","ReadPosRankSum","AN"))
  ### read in file
  vcf = readVcf(filename,genome,param)
  X[[24]] = pDNM_Y(vcf,pat.col,child.col)
  cat("done\n")
  ##########################
  #save(X,file="test.Rdata")
  ##########################
  X = do.call("rbind",X)
 }else{
  rg = GRanges(seqnames=paste(chrom.conv,"X",sep=""),IRanges(start=1,end=3e8))
  ### set params
  param <- ScanVcfParam(which=rg,geno=c("AD","DP","GQ","GT","PL"),
  info = c("VQSLOD","DB","BaseQRankSum","Dels","FS",
    "HaplotypeScore","MQ","MQRankSum",
    "QD","ReadPosRankSum","AN"))
  ### read in file
  vcf = readVcf(filename,genome,param)
  X[[23]] = pDNM_X_female(vcf,pat.col,mat.col,child.col)
  X = do.call("rbind",X)
  ## report progress
  elapsed = proc.time()[3] - time.st
  tot.est = elapsed/prop.done[23]
  remaining = round(tot.est-elapsed,0)
  cat(paste("estimated completion in",round(remaining/60,1),"minutes [",
  format(Sys.time()+remaining, "%a %b %d %X"),"] \n"))

 }
 ## remove duplicates (if any)
 X = X[!duplicated(paste(X$chr,X$pos,sep=":")),]
 return(X)
}
