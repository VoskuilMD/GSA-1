# Author: Michiel Voskuil 

Date: 2018/05/017

0. Necessary files and scripts
---------------------------------------------------
```
# Make sure you have the following scripts/files available in your workdirectory 

Scripts
create_phewas_rscripts.sh            <---- enter your work directory in here 
create_phewas_jobs.sh                <---- enter your work directory in here 
combined_phewas_manhatthanplot.r


Files
Binary plink files
List with multi-allelic sites [i.e. GSA_chr_1-22-merge.missnp]
Michiel_PSI_data_opgeschoond.dta
pheinfo.csv 
phemap.csv
PSI.imputation.key
```


0. Set variables
---------------------------------------------------
```
wd=[workingdirectory] #/groups/umcg-weersma/tmp03/Michiel/phewas/chunks
input=[prefix-to-plink-binary-file]  #GSA_chr_1-22
```

0. Prune SNPs
---------------------------------------------------
```
# # This should actually be done before the association testing
cd $wd

# Prune our dataset using with a 500kb window, 5 step, and 0.3 R2. 
module load plink

plink --bfile /groups/umcg-weersma/tmp04/Michiel/GSA-redo/imputation/european/results/european_maf001/mergedplinkfiles/$input --out $wd/$input-prune --indep-pairwise 500 5 0.3  

# 7737655 variants before pruning

plink --bfile /groups/umcg-weersma/tmp04/Michiel/GSA-redo/imputation/european/results/european_maf001/mergedplinkfiles/$input --extract $wd/$input-prune.prune.in --make-bed --out $wd/$input

# 1085628 variants after pruning
```

1. Prepare genotype files
---------------------------------------------------
```
module load plink
cd $wd

# Remove multi-allelic sites
plink --bfile $input --exclude GSA_chr_1-22-merge.missnp --out tmp01 --make-bed --allow-no-sex

# Cut into small chunks [here chunks of 10000 snps]
mkdir splits

split --lines=10000 tmp01.bim --numeric-suffixes -a 3 splits/chunk_

ls splits/* | while read line; 
    do plink --bfile tmp01 --extract <(awk '{print $2}' ${line}) --recode A --out ${line} --allow-no-sex;
done

# Remove unuseful stuff and move to correct directory
mv splits/*.raw .
rm splits/*.nosex
rm splits/*.log
rm tmp01.*

# create jobs for phewas
bash create_phewas_jobs.sh
mv splits/chunk_*.sh .

# create phewas r scripts
bash create_phewas_rscripts.sh
mv splits/chunk_*.r .
```

2. Perform phewas
---------------------------------------------------
```
cd $wd

# Run phewas on cluster 
# Make sure you have got all the necessary files in your workdirectory

for j in chunk_*.sh; do sbatch "$j"; done

## on boxy this took ~20 hours for all snps
```


3. Combine results in one results file
---------------------------------------------------
```
cd $wd

# Combine
head -1 phewasresults/PheWas_chunk_000.csv > phewasresults/PheWas_All_snps.csv; tail -n +2 -q phewasresults/PheWas_chunk*.csv >> phewasresults/PheWas_All_snps.csv

# Remove NAs from PheWas_All_snps.csv, as there are many NA results for phenotypes with < 20 cases or monomorphic snps

awk 'BEGIN{FS=OFS=","} $4!="NA"{print $0}' phewasresults/PheWas_All_snps.csv > phewasresults/PheWas_All_snps_without_NA.csv
awk 'BEGIN{FS=OFS=","} $4=="NA"{print $0}' phewasresults/PheWas_All_snps.csv > phewasresults/PheWas_All_snps_NA.csv


# Make combined Manhattan plot
Rscript combined_phewas_manhatthanplot.r

# Clean temporary stuff
rm chunk_*.raw
rm chunk_*.sh
rm chunk_*.r
rm chunk_*.err
rm chunk_*.out
rm -r splits
```
 
 
 99. GWAS cat snps
---------------------------------------------------
```
# Optionally, you can subset only GWAS catalogue snps

# module load VCFtools
# vcftools --gzvcf /groups/umcg-weersma/tmp04/Michiel/GSA-redo/imputation/european/results/european_maf001/merged/european_maf001.vcf.gz --positions /groups/umcg-weersma/tmp04/Michiel/GSA-redo/phewas/curated_snplist --recode --stdout | bgzip -c > /groups/umcg-weersma/tmp04/Michiel/GSA-redo/phewas/european_maf001_gwascatalogue.vcf.gz
# tabix -p vcf european_maf001_gwascatalogue.vcf.gz

# Here 35047 snps left.

# Make plink files from imputed VCFs (use prob > 0.4)

# sbatch /groups/umcg-weersma/tmp04/Michiel/GSA-redo/scripts/imputedVCFtoPlinkS_gwascat.sh

# module load plink
# plink --bfile /groups/umcg-weersma/tmp04/Michiel/GSA-redo/phewas/plink/GSA --recodeA --out /groups/umcg-weersma/tmp04/Michiel/GSA-redo/phewas/PheWASanalysis/GWAScatsnps/GSA.A
```
