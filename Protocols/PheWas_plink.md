# Do phewas sort of manual. That is: perform several GWAS's for all different phenotypes and see whether results correlate well.

ml plink
wd=/groups/umcg-weersma/tmp04/Michiel/GSA-redo/phewas/plinkanalyses
inputfolder=/groups/umcg-weersma/tmp04/Michiel/GSA-redo/imputation/european/results/european_maf001/mergedplinkfiles
input=GSA_chr_1-22

cd $wd

plink --bfile $inputfolder/$input --update-sex /groups/umcg-weersma/tmp04/Michiel/GSA-redo/GSA.sex.info.postimp --out tmp --make-bed

# make covariate file including diagnosis and first 5 genetic PCA

plink --bfile tmp --pca 10 header tabs --out genetic.pca


ml R
R
# in R
gen.cov = read.table('genetic.pca.eigenvec', header=T)
diag.cov = read.table('covariates.raw.txt', header=T)
mergefile = merge(diag.cov,gen.cov,by=c('FID','IID')
write.table(mergefile,'covariates.diag.gen.txt',row.names=F,col.names=T,sep='\t',quote=F)

# distinguish binary phenotypes from quantitative phenotypes

# identify CD and UC/IBDU samples
awk -F"\t" '$3 == "TRUE" { print $1"\t"$2 }' covariates.diag.gen.txt > CD.samples
awk -F"\t" '$3 == "FALSE" { print $1"\t"$2 }' covariates.diag.gen.txt > UC.samples


#binary
awk 'BEGIN {OFS = "\t"} { print $1,$2,$3,$4,$5,$6,$7,$10,$11,$12,$13,$14,$15,$17,$18,$19,$20,$21,$22,$23,$24,$25,$26,$27,$28,$29,$30,$31,$34,$35,$36,$37,$38,$39,$40,$41,$42}' phenotypes.raw.txt > binary.phenotypes.txt

#quantitative
awk 'BEGIN {OFS = "\t"} { print $1,$2,$8,$9,$16,$32,$33,$43,$44,$45,$46,$47,$48,$49,$50,$51,$52,$53,$54,$55}' phenotypes.raw.txt > quantitative.phenotypes.txt

# do actual logistic/linear regression for all phenotypes with diagnosis and first 5 PC's as covariates. 
#binary/logistic (screen -r 21978.pts-9.calculon/screen -r 22545.pts-15.calculon )
plink --bfile tmp --covar covariates.diag.gen.txt --covar-number 4-9 --ci 0.95 --pheno binary.phenotypes.txt --all-pheno --out binary.hide.covar --logistic hide-covar --allow-no-sex
#quantitative/linear (screen -r 7773.pts-27.calculon/screen -r 24976.pts-38.calculon) 
plink --bfile tmp --covar covariates.diag.gen.txt --covar-number 4-9 --ci 0.95 --pheno quantitative.phenotypes.txt --all-pheno --out quantitative.hide.covar --linear hide-covar --allow-no-sex


# perform analyses seperately per diagnosis to see if indeed correlation between phenotype and covariates (diagnosis) causes trouble)
awk  'BEGIN{FS="\t";OFS="\t"}{ print $1,$2,$7,$8,$9,$10,$11 }' covariates.diag.gen.txt > covariates.gen.txt

#binary/logistic (screen -r 14976.pts-13.calculon)
plink --bfile tmp --covar covariates.gen.txt --keep CD.samples --covar-number 1-5 --ci 0.95 --pheno binary.phenotypes.txt --all-pheno --out binary.hide.covar.CD --logistic hide-covar --allow-no-sex
#quantitative/linear (screen -r 18258.pts-37.calculon) 
plink --bfile tmp --covar covariates.gen.txt --keep CD.samples --covar-number 1-5 --ci 0.95 --pheno quantitative.phenotypes.txt --all-pheno --out quantitative.hide.covar.CD --linear hide-covar --allow-no-sex

#binary/logistic (screen -r 21497.pts-47.calculon)
plink --bfile tmp --covar covariates.gen.txt --keep UC.samples --covar-number 1-5 --ci 0.95 --pheno binary.phenotypes.txt --all-pheno --out binary.hide.covar.UC --logistic hide-covar --allow-no-sex
#quantitative/linear (screen -r 23862.pts-48.calculo) 
plink --bfile tmp --covar covariates.gen.txt --keep UC.samples --covar-number 1-5 --ci 0.95 --pheno quantitative.phenotypes.txt --all-pheno --out quantitative.hide.covar.UC --linear hide-covar --allow-no-sex


# output files will look like this: {binary,quantitative}.hide.covar.P*.assoc.{logistic,linear}

# sort based on p-value



# The .assoc files contain some NA's. These are for SNPs that are either in the cases or controls monomorphic (usually cases of course)
# If all are NA then the model does not fit, so I have now run CD and UC separately, to circumvent this, but of course we half power... meta-??






# for cox regression on survival
# https://www.liverpool.ac.uk/translational-medicine/research/statistical-genetics/survival-gwas-sv/
