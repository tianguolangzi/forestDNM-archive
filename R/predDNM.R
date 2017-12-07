predDNM <-
function(X,rf,cutoff=0.6){
 ## remove observations with NA values in features
 ind = sapply(X,function(x) is.na(x))
 ind = !rowSums(ind)>0
 if(any(!ind)){
  X = X[ind,]
  warning(paste("removed",sum(!ind),"SNVs due to NAs"))
 }
 ## remove chr, pos, and allele columns
 info = X[,colnames(X)%in%c("chr","pos","par_allele","mut_allele")]
 X = X[,!colnames(X)%in%c("chr","pos","par_allele","mut_allele")]
 
 ## make predictions
 info$prd = predict(rf,X,type="prob")[,1]
 
 ## only return those at or above the cutoff
 info = info[info$prd >= cutoff,]

 return(info)
}
