library(PheWAS)

# Read combined results
setwd('phewasresults/')
results = read.csv('PheWas_All_snps_without_NA.csv', sep = ',', colClasses=c("phenotype"="character"))

# Supporting files for plot

pheinfo = read.csv("pheinfo.csv", sep = ";", colClasses=c("character",rep("character",4)))
phemap= read.csv("phemap.csv", sep = ";", colClasses = c("character", "character"))
phemap[3]=NULL
phemap[3]=NULL
phemap[3]=NULL
phemap[3]=NULL

# Plot results
pdf("PheWas_All_snps.pdf")
phewasManhattan(results, annotate.phenotype.description = T, title = "PheWas All snps")
dev.off()
jpeg("PheWas_All_snps.jpg")
phewasManhattan(results, annotate.phenotype.description = T, title = "PheWAS All snps")
dev.off()
