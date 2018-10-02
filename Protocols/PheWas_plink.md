# Author: Michiel Voskuil 

# Date: 2018/10/01


1. Do phewas first manually. That is: perform several GWAS's for all different phenotypes and see whether results correlate well. Step 1: prepare genotype files:
---------------------------------------------------
```
# Copy imputed genotype files to your workdirectory
ml plink
wd=/groups/umcg-weersma/tmp04/Michiel/GSA-redo/phewas/plink_per_phenotype
inputfolder=/groups/umcg-weersma/tmp04/Michiel/GSA-redo/imputation/european/results/european_maf001/mergedplinkfiles
input=GSA_chr_1-22

cd $wd
# Update sex
plink --bfile $inputfolder/$input --update-sex /groups/umcg-weersma/tmp04/Michiel/GSA-redo/GSA.sex.info.postimp --out tmp --make-bed

# Update map files: convert chr:bp to rsids.

# make covariate file including diagnosis and first 5 genetic PCA
plink --bfile tmp --pca 10 header tabs --out genetic.pca

ml R
R
# in R
gen.cov = read.table('genetic.pca.eigenvec', header=T)
diag.cov = read.table('covariates.raw.txt', header=T)
mergefile = merge(diag.cov,gen.cov,by=c('FID','IID')
write.table(mergefile,'covariates.diag.gen.txt',row.names=F,col.names=T,sep='\t',quote=F)

# make covariate files also without diagnosis
awk  'BEGIN{FS="\t";OFS="\t"}{ print $1,$2,$7,$8,$9,$10,$11 }' covariates.diag.gen.txt > covariates.gen.txt

# distinguish binary phenotypes from quantitative phenotypes

# identify CD and UC/IBDU samples
awk -F"\t" '$3 == "TRUE" { print $1"\t"$2 }' covariates.diag.gen.txt > CD.samples
awk -F"\t" '$3 == "FALSE" { print $1"\t"$2 }' covariates.diag.gen.txt > UC.samples

#update CHR:BP snp names to rsides using http://hgdownload.cse.ucsc.edu/goldenpath/hg19/database/snp150.txt.gz
wget http://hgdownload.cse.ucsc.edu/goldenpath/hg19/database/snp150.txt.gz > snp150.txt.gz 
zcat snp150.txt.gz | awk 'OFS="\t" {print $2,$4,$5}' | sed 's/^.\{,3\}//' | awk 'OFS="\t" {print $1 ":" $2,$3}' > rsid.list

# remove duplicate sites (tri-allelic), plink only handles bi-allelic sites
awk '{seen[$1]++; a[++count]=$0; key[count]=$1} END {for (i=1;i<=count;i++) if (seen[key[i]] == 1) print a[i]}' rsid.list > rsid.list2
```

2. Step 2: prepare phenotype files:
---------------------------------------------------
```cd /groups/umcg-weersma/tmp04/Michiel/GSA-redo/phewas/plink_per_phenotype
head binary.phenotypes.txt > binary.phenotypes.names -n1
head quantitative.phenotypes.txt > quantitative.phenotypes.names -n1

# Make separate binary phenotype files
for i in {Colectomy,Stenosing,Penetrating,PeriAnalDisease,Ileocaecal_resection,Smoking_EN,Smoking_CE,Complications,EIM_arthropathy,EIM_arthritis,Pouchitis,A1,A2,A3,PerianalDisease,E1,E2,E3,Azathioprine,Mercaptopurine,Immunomodulator,Mesalazine,PSC,Appendectomy,Pouch,Stoma,Uveitis,Erythema,Pyoderma,OralAphthae,AnalFissura,Skin,Eyes,TromboticEvents,EIM_BMD};
do plink --bfile tmp --pheno-name "$i" --pheno binary.phenotypes.txt --allow-no-sex --make-bed --out "$i"
rm "$i".bed
rm "$i".bim;
done

# Make separate quantitative phenotype files
for i in {CD_Time_to_Surgery,UC_Time_to_Surgery,AgeDiagnosis,HBImean,SCCAImean,Height,ASAT,AF,ALAT,BSE,CRP,GGT,Ht,Leuco,MCV,Creat,Thrombos,Hb};
do plink --bfile tmp --pheno-name "$i" --pheno quantitative.phenotypes.txt --allow-no-sex --make-bed --out "$i"
rm "$i".bed
rm "$i".bim;
done
```

Step 3: do GWAS's
---------------------------------------------------
```
# Do association test per binary phenotype 
for i in {Colectomy,Stenosing,Penetrating,PeriAnalDisease,Ileocaecal_resection,Smoking_EN,Smoking_CE,Complications,EIM_arthropathy,EIM_arthritis,Pouchitis,A1,A2,A3,PerianalDisease,E1,E2,E3,Azathioprine,Mercaptopurine,Immunomodulator,Mesalazine,PSC,Appendectomy,Pouch,Stoma,Uveitis,Erythema,Pyoderma,OralAphthae,AnalFissura,Skin,Eyes,TromboticEvents,EIM_BMD}; do
plink --bed tmp.bed --bim tmp.bim --fam "$i".fam --logistic --covar covariates.gen.txt --covar-number 1-5 -out "$i" --allow-no-sex;
done

# Do association testing per quantitative phenotype 
for i in {CD_Time_to_Surgery,UC_Time_to_Surgery,AgeDiagnosis,HBImean,SCCAImean,Height,ASAT,AF,ALAT,BSE,CRP,GGT,Ht,Leuco,MCV,Creat,Thrombos,Hb}; do
plink --bed tmp.bed --bim tmp.bim --fam "$i".fam --linear --covar covariates.gen.txt --covar-number 1-5 -out "$i" --allow-no-sex;
done
```

Step 4: analyse results: select only additive results from .assoc.{logistic,linear} and disregard PC1-5 assoc results. I have removed NA's here from the results as well. 
---------------------------------------------------
```
for i in {Colectomy,Stenosing,Penetrating,PeriAnalDisease,Ileocaecal_resection,Smoking_EN,Smoking_CE,Complications,EIM_arthropathy,EIM_arthritis,Pouchitis,A1,A2,A3,E1,E2,E3,Azathioprine,Mercaptopurine,Immunomodulator,Mesalazine,PSC,Appendectomy,Pouch,Stoma,Uveitis,Erythema,Pyoderma,OralAphthae,AnalFissura,Skin,Eyes,TromboticEvents,EIM_BMD};
do awk '(($5 == "ADD" || $5 == "TEST") && $9!="NA")' "$i".assoc.logistic > "$i".assoc.logistic_tmp1
sort -gk9 "$i".assoc.logistic_tmp1 > "$i".assoc.logistic_tmp2
cp *_tmp2 filtered_results; 
done

for i in {CD_Time_to_Surgery,UC_Time_to_Surgery,AgeDiagnosis,HBImean,SCCAImean,Height,ASAT,AF,ALAT,BSE,CRP,GGT,Ht,Leuco,MCV,Creat,Thrombos,Hb};
do awk '(($5 == "ADD" || $5 == "TEST") && $9!="NA")' "$i".assoc.linear > "$i".assoc.linear_tmp1
sort -gk9 "$i".assoc.linear_tmp1 > "$i".assoc.linear_tmp2
cp *_tmp2 filtered_results; 
done

bash create_manhattan.sh
bash create_manhattan_rscripts.sh

for i in {Colectomy,Stenosing,Penetrating,PeriAnalDisease,Ileocaecal_resection,Smoking_EN,Smoking_CE,Complications,EIM_arthropathy,EIM_arthritis,Pouchitis,A1,A2,A3,E1,E2,E3,Azathioprine,Mercaptopurine,Immunomodulator,Mesalazine,PSC,Appendectomy,Pouch,Stoma,Uveitis,Erythema,Pyoderma,OralAphthae,AnalFissura,Skin,Eyes,TromboticEvents,EIM_BMD};
do sbatch Manhattan_"$i".sh; done

for i in {CD_Time_to_Surgery,UC_Time_to_Surgery,AgeDiagnosis,HBImean,SCCAImean,Height,ASAT,AF,ALAT,BSE,CRP,GGT,Leuco,MCV,Creat,Thrombos,Hb};
do sbatch Manhattan_"$i".sh; done
```



Step (extra): to do
---------------------------------------------------
```
#For some reason Ht doesn't yield results. Have checked the reason yet though

#about results
Good way to visualize Manhattan way: https://www.r-graph-gallery.com/wp-content/uploads/2018/02/Manhattan_plot_in_R.html 

# for cox regression on survival
https://www.liverpool.ac.uk/translational-medicine/research/statistical-genetics/survival-gwas-sv/
```




deprecated
---------------------------------------------------
```
# binary
awk 'BEGIN {OFS = "\t"} { print $1,$2,$3,$4,$5,$6,$7,$10,$11,$12,$13,$14,$15,$17,$18,$19,$20,$21,$22,$23,$24,$25,$26,$27,$28,$29,$30,$31,$34,$35,$36,$37,$38,$39,$40,$41,$42}' phenotypes.raw.txt > binary.phenotypes.txt

# quantitative
awk 'BEGIN {OFS = "\t"} { print $1,$2,$8,$9,$16,$32,$33,$43,$44,$45,$46,$47,$48,$49,$50,$51,$52,$53,$54,$55}' phenotypes.raw.txt > quantitative.phenotypes.txt

# do actual logistic/linear regression for all phenotypes with diagnosis and first 5 PC's as covariates. 
plink --bfile tmp --covar covariates.diag.gen.txt --covar-number 4-9 --ci 0.95 --pheno binary.phenotypes.txt --all-pheno --out binary.hide.covar --logistic hide-covar --allow-no-sex

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


#######################################
#last updated: 20/09/2018

#Om het wat overzichtelijker te maken probeer ik eerst even allemaal losse GWAS met ook de titel van de phenotypes in de #filename. 

cd /groups/umcg-weersma/tmp04/Michiel/GSA-redo/phewas/plink_per_phenotype
head binary.phenotypes.txt > binary.phenotypes.names -n1
head quantitative.phenotypes.txt > quantitative.phenotypes.names -n1

# update CHR:BP snp names to rsides using http://hgdownload.cse.ucsc.edu/goldenpath/hg19/database/snp150.txt.gz

wget http://hgdownload.cse.ucsc.edu/goldenpath/hg19/database/snp150.txt.gz > snp150.txt.gz 

zcat snp150.txt.gz | awk 'OFS="\t" {print $2,$4,$5}' | sed 's/^.\{,3\}//' | awk 'OFS="\t" {print $1 ":" $2,$3}' > test.txt

sed 's/^.\{,3\}//' > test

| awk 'OFS="\t"  {print $1, $1":"$2 }' | awk 'OFS="\t" {print $4,$3}' | gzip -9c > 

# To update the plink files
gunzip 



gunzip -9c ChrBp_to_rsid.txt.gz
# Remove dupblicates first (that is, positions with multiple rsids -> we want to keep only bi-allelic sites)
awk '{seen[$1]++; a[++count]=$0; key[count]=$1} END {for (i=1;i<=count;i++) if (seen[key[i]] == 1) print a[i]}' ChrBp_to_rsid.txt > ChrBp_to_rsid_unique.txt
# Update map files
plink --bfile tmp --update-map ChrBp_to_rsid_uniqe-update-name 


# Make separate binary phenotype files
for i in {Colectomy,Stenosing,Penetrating,PeriAnalDisease,Ileocaecal_resection,Smoking_EN,Smoking_CE,Complications,EIM_arthropathy,EIM_arthritis,Pouchitis,A1,A2,A3,PerianalDisease,E1,E2,E3,Azathioprine,Mercaptopurine,Immunomodulator,Mesalazine,PSC,Appendectomy,Pouch,Stoma,Uveitis,Erythema,Pyoderma,OralAphthae,AnalFissura,Skin,Eyes,TromboticEvents,EIM_BMD};
do plink --bfile tmp --pheno-name "$i" --pheno binary.phenotypes.txt --allow-no-sex --make-bed --out "$i"
rm "$i".bed
rm "$i".bim;
done

# Make separate quantitative phenotype files
for i in {CD_Time_to_Surgery,UC_Time_to_Surgery,AgeDiagnosis,HBImean,SCCAImean,Height,ASAT,AF,ALAT,BSE,CRP,GGT,Ht,Leuco,MCV,Creat,Thrombos,Hb};
do plink --bfile tmp --pheno-name "$i" --pheno quantitative.phenotypes.txt --allow-no-sex --make-bed --out "$i"
rm "$i".bed
rm "$i".bim;
done

# Generate genetic covariate file
cp ../plinkanalyses/covariates.gen.txt .

# Do association test per binary phenotype 
for i in {Colectomy,Stenosing,Penetrating,PeriAnalDisease,Ileocaecal_resection,Smoking_EN,Smoking_CE,Complications,EIM_arthropathy,EIM_arthritis,Pouchitis,A1,A2,A3,PerianalDisease,E1,E2,E3,Azathioprine,Mercaptopurine,Immunomodulator,Mesalazine,PSC,Appendectomy,Pouch,Stoma,Uveitis,Erythema,Pyoderma,OralAphthae,AnalFissura,Skin,Eyes,TromboticEvents,EIM_BMD}; do
plink --bed tmp.bed --bim tmp.bim --fam "$i".fam --logistic --covar covariates.gen.txt --covar-number 1-5 -out "$i" --allow-no-sex;
done

# Do association testing per quantitative phenotype 
for i in {CD_Time_to_Surgery,UC_Time_to_Surgery,AgeDiagnosis,HBImean,SCCAImean,Height,ASAT,AF,ALAT,BSE,CRP,GGT,Ht,Leuco,MCV,Creat,Thrombos,Hb}; do
plink --bed tmp.bed --bim tmp.bim --fam "$i".fam --linear --covar covariates.gen.txt --covar-number 1-5 -out "$i" --allow-no-sex;
done

# Select only additive results from .assoc.{logistic,linear} and disregard PC1-5 assoc results
# I have removed NA's here from the results as well. 

for i in {Colectomy,Stenosing,Penetrating,PeriAnalDisease,Ileocaecal_resection,Smoking_EN,Smoking_CE,Complications,EIM_arthropathy,EIM_arthritis,Pouchitis,A1,A2,A3,E1,E2,E3,Azathioprine,Mercaptopurine,Immunomodulator,Mesalazine,PSC,Appendectomy,Pouch,Stoma,Uveitis,Erythema,Pyoderma,OralAphthae,AnalFissura,Skin,Eyes,TromboticEvents,EIM_BMD};
do awk '(($5 == "ADD" || $5 == "TEST") && $9!="NA")' "$i".assoc.logistic > "$i".assoc.logistic_tmp1
sort -gk9 "$i".assoc.logistic_tmp1 > "$i".assoc.logistic_tmp2
cp *_tmp2 filtered_results; 
done

for i in {CD_Time_to_Surgery,UC_Time_to_Surgery,AgeDiagnosis,HBImean,SCCAImean,Height,ASAT,AF,ALAT,BSE,CRP,GGT,Ht,Leuco,MCV,Creat,Thrombos,Hb};
do awk '(($5 == "ADD" || $5 == "TEST") && $9!="NA")' "$i".assoc.linear > "$i".assoc.linear_tmp1
sort -gk9 "$i".assoc.linear_tmp1 > "$i".assoc.linear_tmp2
cp *_tmp2 filtered_results; 
done

bash create_manhattan.sh
bash create_manhattan_rscripts.sh

for i in {Colectomy,Stenosing,Penetrating,PeriAnalDisease,Ileocaecal_resection,Smoking_EN,Smoking_CE,Complications,EIM_arthropathy,EIM_arthritis,Pouchitis,A1,A2,A3,E1,E2,E3,Azathioprine,Mercaptopurine,Immunomodulator,Mesalazine,PSC,Appendectomy,Pouch,Stoma,Uveitis,Erythema,Pyoderma,OralAphthae,AnalFissura,Skin,Eyes,TromboticEvents,EIM_BMD};
do sbatch Manhattan_"$i".sh; done

for i in {CD_Time_to_Surgery,UC_Time_to_Surgery,AgeDiagnosis,HBImean,SCCAImean,Height,ASAT,AF,ALAT,BSE,CRP,GGT,Leuco,MCV,Creat,Thrombos,Hb};
do sbatch Manhattan_"$i".sh; done


# For some reason Ht doesn't yield results. Have checked the reason yet though


# about results

Good way to visualize Manhattan way: https://www.r-graph-gallery.com/wp-content/uploads/2018/02/Manhattan_plot_in_R.html 



# for cox regression on survival
# https://www.liverpool.ac.uk/translational-medicine/research/statistical-genetics/survival-gwas-sv/
```
