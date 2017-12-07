pDNM_X_female <-
function(vcf,pat.col=1,mat.col=2,child.col=3){
 sn = as.character(seqnames(vcf))
 if(!all(sn%in%c("X","chrX"))) stop("this function intended for chrX only") 
 st = start(vcf)
 g = geno(vcf)$GT
 p1 = g[,pat.col] == g[,mat.col] ## parents have same genotype
 p2 = g[,pat.col] == "0/0" | g[,pat.col] == "1/1" | g[,pat.col] == "2/2" ## parents are homozygous
 c1 = g[,child.col] == "0/1" | g[,child.col] == "0/2" | 
	g[,child.col] == "1/2" ## child is het
 ind = p1 & p2 & c1
 
 ## get the index of the parental allele
 ## because of the genotype './.', some entries will be coerced to NA
 ## during the as.integer call; we replace these with '0'
 pa = suppressWarnings(as.integer(substr(g[,pat.col],1,1))+1)

 ## ...and those of the offspring alleles
 ca1 = suppressWarnings(as.integer(substr(g[,child.col],1,1))+1)
 ca2 = suppressWarnings(as.integer(substr(g[,child.col],3,3))+1)

 ## replace NAs with '0'
 pa[is.na(pa)] = 0L
 ca1[is.na(ca1)] = 0L
 ca2[is.na(ca2)] = 0L

 ## the index of the mutant allele in the offspring
 mut = rep(NA,length(ca1))
 mut[pa!=ca1] = ca1[pa!=ca1]
 mut[pa!=ca2] = ca2[pa!=ca2]
 
 ## in case any entries in mut are '0', set them back to NA
 mut[mut<1] = NA

 ## sanity check: at least one of the offspring alleles 
 ## must originate from the parents
 ind = ind & (pa==ca1 | pa==ca2) & (pa >= 1L)

 ## from now on we'll only deal with the data that are
 ## putative DNMs
 REF = ref(vcf)[ind]
 ALT = alt(vcf)[ind]
 mut = mut[ind]
 pa = pa[ind]
 ca1 = ca1[ind]
 ca2 = ca2[ind]
 sn = as.vector(sn[ind])
 st = st[ind]
 AD = geno(vcf)$AD[ind,] ## 3 cols (of lists)
 PL = geno(vcf)$PL[ind,] ## 3 cols (of lists)
 GQ = geno(vcf)$GQ[ind,] ## 3 cols
 GT = geno(vcf)$GT[ind,] ## 3 cols

 ### get the actual base for the parental and mutant alleles
 par_allele = rep("",length(REF))
 par_allele[pa==1] = as.character(REF[pa==1])
 par_allele[pa==2] = sapply(ALT[pa==2],function(x) as.character(x)[1])
 par_allele[pa==3] = sapply(ALT[pa==3],function(x) as.character(x)[2])

 mut_allele = rep("",length(REF))
 mut_allele[mut==1] = as.character(REF[mut==1])
 mut_allele[mut==2] = sapply(ALT[mut==2],function(x) as.character(x)[1])
 mut_allele[mut==3] = sapply(ALT[mut==3],function(x) as.character(x)[2])

 ##################################################################
 ###### D E B U G #################################################
 if(is.list(mut_allele) && any(sapply(mut_allele,length)>1)){
  save(mut_allele,file="debug.Rdata")
  stop("something weird with mut_allele; dumped to debug.Rdata")
 }
 if(is.list(par_allele) && any(sapply(par_allele,length)>1)){
  save(par_allele,file="debug.Rdata")
  stop("something weird with mut_allele; dumped to debug.Rdata")
 }
 ##################################################################
 ##################################################################

 ## ALLELE COUNT FEATURES
 
 ## get the mutant allele counts in the parents
 p_mut1 = mapply(function(x,y) x[y],AD[,pat.col],mut)
 p_mut2 = mapply(function(x,y) x[y],AD[,mat.col],mut)
 p_par1 = mapply(function(x,y) x[y],AD[,pat.col],pa) + 1
 p_par2 = mapply(function(x,y) x[y],AD[,mat.col],pa) + 1

 ## get the mutant and parental allele counts in the child
 c_par = mapply(function(x,y) x[y],AD[,child.col],pa) + 1
 c_mut = mapply(function(x,y) x[y],AD[,child.col],mut)

 ## get median coverage for mother, father, and child
 cvg = geno(vcf)$AD
 m_cvg_med = median(sapply(cvg[,mat.col],sum,na.rm=T))
 p_cvg_med = median(sapply(cvg[,pat.col],sum,na.rm=T))
 c_cvg_med = median(sapply(cvg[,child.col],sum,na.rm=T))

 ## the log2 coverage ratio (relative to median)
 m_cvg = log2((p_par2+p_mut2)/(m_cvg_med+1))
 p_cvg = log2((p_par1+p_mut1)/(p_cvg_med+1))
 c_cvg = log2((c_par+c_mut)/(c_cvg_med+1))

 AD_X = data.frame(p_ar_max=pmax(p_mut1/p_par1,p_mut2/p_par2,na.rm=T),
			p_ar_min=pmin(p_mut1/p_par1,p_mut2/p_par2,na.rm=T),
			c_ar=c_mut/c_par,
			p_cvg_max=pmax(m_cvg,p_cvg),
			p_cvg_min=pmin(m_cvg,p_cvg),
			c_cvg=c_cvg)


 ## GENOTYPE LIKELIHOOD FEATURES

 ## format the Phred genotype likelihoods: homozygous parental
 p_homo_PL1 = mapply(function(x,y) x[c(1,3,6)][y], PL[,pat.col], pa) 
 p_homo_PL2 = mapply(function(x,y) x[c(1,3,6)][y], PL[,mat.col], pa) 
 c_homo_PL = mapply(function(x,y) x[c(1,3,6)][y], PL[,child.col], pa)

 ## format the Phred genotype likelihoods: het (child genotype)
 c_gt = rep(NA,nrow(GT))
 c_gt[GT[,child.col]=="0/1"] = 1
 c_gt[GT[,child.col]=="0/2"] = 2
 c_gt[GT[,child.col]=="1/2"] = 3
 
 p_het_PL1 = mapply(function(x,y) x[c(2,4,5)][y], PL[,pat.col], c_gt) 
 p_het_PL2 = mapply(function(x,y) x[c(2,4,5)][y], PL[,mat.col], c_gt) 
 c_het_PL = mapply(function(x,y) x[c(2,4,5)][y], PL[,child.col], c_gt)

 PL_X = data.frame(p_cg_max=pmax(p_het_PL1,p_het_PL2,na.rm=T),
		p_cg_min=pmin(p_het_PL1,p_het_PL2,na.rm=T),
		c_cg=c_het_PL,
		p_pg_max=pmax(p_homo_PL1,p_homo_PL2,na.rm=T),
		p_pg_min=pmin(p_homo_PL1,p_homo_PL2,na.rm=T),
		c_pg=c_homo_PL)
	

 ## OTHER FEATURES
 X = data.frame(chr=sn,pos=st,
		 FL = fixed(vcf)$FILTER[ind], ## character vector
		 QL = fixed(vcf)$QUAL[ind], ## vector
		 VQ = info(vcf)$VQSLOD[ind], ## vector
		 #SB = info(vcf)$SB[ind], ## vector
		 BZ = info(vcf)$BaseQRankSum[ind], ## vector
		 DL = info(vcf)$Dels[ind], ## vector
		 FS = info(vcf)$FS[ind], ## vector
		 HS = info(vcf)$HaplotypeScore[ind], ## vector
		 MQ = info(vcf)$MQ[ind], ## vector
		 MZ = info(vcf)$MQRankSum[ind], ## vector
		 QD = info(vcf)$QD[ind], ## vector
		 PZ = info(vcf)$ReadPosRankSum[ind],
		 GQ_p_min=pmin(geno(vcf)$GQ[ind,pat.col],
			geno(vcf)$GQ[ind,mat.col],na.rm=T),
		 GQ_p_max=pmax(geno(vcf)$GQ[ind,pat.col],
			geno(vcf)$GQ[ind,mat.col],na.rm=T),
		 GQ_c=geno(vcf)$GQ[ind,child.col],
		 N_alt=sapply(ALT,length),
		 par_allele=unlist(par_allele),
		 mut_allele=unlist(mut_allele)) ## vector)

 X = cbind(AD_X,PL_X,X)
 rownames(X) = paste(sn,st,sep=":")
 return(X)
}
