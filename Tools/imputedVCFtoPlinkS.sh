#!/bin/bash
#SBATCH --job-name=vcfToPlink
#SBATCH --nodes=1
#SBATCH --time=23:59:59
#SBATCH --mem=60gb
#SBATCH --output=vcfToPlink.out
#SBATCH --error=vcfToPlink.err

module load tabix
module load GenotypeHarmonizer
module load plink

cohort="GSA"
input="/groups/umcg-weersma/tmp04/Michiel/GSA-redo/imputation/european/results/european_GRCh37_maf001/annotated/"
output="/groups/umcg-weersma/tmp04/Michiel/GSA-redo/imputation/european/results/european_GRCh37_maf001/annotated/plinkfiles"

mkdir -p $output

for chr in {1..22}
do

GenotypeHarmonizer.sh \
    -i "${input}chr_${chr}" \
    -I VCFFOLDER \
    -o "${output}${cohort}_chr_${chr}" \
    -O PLINK_BED 

done

echo "Finished converting VCF files into Plink BED files from ${cohort} cohort"


mergedOutPut="/groups/umcg-weersma/tmp04/Michiel/GSA-redo/imputation/european/results/european_GRCh37_maf001/annotated/mergedPlink/"


mkdir -p $mergedOutPut

cd $output
ls *.bed > bed.txt
ls *.bim > bim.txt
ls *.fam > fam.txt

paste bed.txt bim.txt fam.txt | column -s $'\t' -t > merge.list

plink --bfile "${output}${cohort}_chr_1" \
--merge-list merge.list \
--make-bed \
--out "${mergedOutPut}${cohort}"

echo "Finished concatenating all chromosome for ${cohort} cohort"
