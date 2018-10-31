## Read-in data

MAP_serology = read.table("Documents/Werk/Promotie/MAP/MAP_data_serologie_2.csv", sep= ";", header=T)

## Remove standardised data columns
MAP_serology[,c(6,8,10,12,13,15,17,19,21,22)] = NULL
row.names(MAP_serology) <- MAP_serology[,1]
MAP_serology[,c(1)] = NULL

## create covariates files
MAP.covariates = as.data.frame(MAP_serology[,c(1,2,3)])
MAP.covariates$Sex = ifelse(MAP.covariates$Sex == "Man", "1", "2")
MAP_serology[,c(1,2,3)] = NULL
MAP.covariates <- cbind(names = rownames(MAP.covariates), MAP.covariates)
colnames(MAP.covariates)[1] <- "SampleID"
MAP.covariates[,c(4)] = NULL

## Drop cases with missing covariates
MAP.covariates=MAP.covariates[c(-41,-53),]
MAP_serology=MAP_serology[c(-41,-53),]
MAP.covariates=MAP.covariates[,-1]

### inverse rank function:
invrank= function(x) {qnorm((rank(x,na.last="keep")-0.5)/sum(!is.na(x)))}
#MAP_serology_clean = as.data.frame(t(MAP_serology))
MAP_serology_logN = as.data.frame(apply(MAP_serology,2,invrank))


#### linear model correction 
MAP_serology_corrected = apply(MAP_serology_logN,2,function(x){
  #x.subset = x[!is.na(x)]
  #covariate.subset = MAP.covariates[!is.na(x),,drop = FALSE]
  #covariate.subset.matrix = data.matrix(covariate.subset)
  #if(ncol(covariate.subset)==ncol(covariate.subset.matrix)){
   # covariate.subset = covariate.subset[,apply(covariate.subset.matrix,2,sd) !=0,drop = FALSE]
  #}
  
  x.resid = resid(lm(x ~ .,data = MAP.covariates))
  x=x.resid+100
  #x[!is.na(x)] = x.resid+100
  #x[is.na(x)] = 0
  return(x)
})

## Creat annotation file
MAP_serology_corrected = as.data.frame(t(MAP_serology_corrected))
MAP_serology_corrected2 = data.frame(rownames(MAP_serology_corrected),MAP_serology_corrected)
colnames(MAP_serology_corrected2)[1] = "-"
annot = data.frame(platform = "RDP",
                   HT12v3.ArrayAddress = rownames(MAP_serology_corrected2),
                   Gene = rownames(MAP_serology_corrected2),
                   Chr = 4,
                   ChrStart = 1000,
                   ChrEnd = 1000)
colnames(annot)[2] = "HT12v3-ArrayAddress"
write.table(MAP_serology_corrected2, file = "~/Documents/Werk/Promotie/MAP/MAP_corrected.txt",sep="\t",row.names = F,quote = F)
write.table(annot,file = "~/Documents/Werk/Promotie/MAP/MAP_corrected.txt.annot",sep="\t",row.names=F,quote = F)
