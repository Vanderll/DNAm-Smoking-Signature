rm(list=ls())

getBeta <- function(a){
  b = 2^a/(1 + 2^a)  
  return(b)
}


#### load the signature CpGs and weights 

mccartney = read.csv(file="/home/vanderll/Norris/DNAm_smokingSignature/McCartney_CpG_weights.csv")	

reese = read.csv(file="/home/vanderll/Norris/DNAm_smokingSignature/Reese_CpG_weights.csv")	

want = unique(c(reese$CpG, mccartney$CpG))

#### load DNAm data 

load(file="/home/biostats_share/Norris/data/methylation/sesame450K.batchAdj.Mmatrix.Rdata")
beta.450 = getBeta(M.sesame.batch)
rm(M.sesame.batch)

load(file="/home/biostats_share/Norris/data/methylation/sesameEPIC.batchAdj.Mmatrix.Rdata")
beta.EPIC = getBeta(M.sesame.batch)
rm(M.sesame.batch)

load(file="/home/vanderll/Norris/DNAm_smokingSignature/MLscore_CpGs.Rdata")

#### grab the CpGs we want and put into smaller datasets to use;

##from 664,614 and 375,020 EPIC and 450K total
beta.450.rauchert = beta.450[which(rownames(beta.450) %in% ML_cpgs),]
beta.EPIC.rauchert= beta.EPIC[which(rownames(beta.EPIC) %in% ML_cpgs),]

save(beta.450.rauchert, beta.EPIC.rauchert, file = "/home/vanderll/Norris/DNAm_smokingSignature/Rauschert_DNAm_beta_datasets.Rdata") 
