

vcfFolder="/groups/umcg-weersma/tmp04/Michiel/GSA-redo/imputation/admixed/results/test/"
outputFolder="/groups/umcg-weersma/tmp04/Michiel/GSA-redo/imputation/admixed/results/test/admixed_maf001/"
snpDB="/groups/umcg-weersma/tmp04/Michiel/GSA-redo/imputation/annotation/GRCh37/All_20170710.vcf.gz"
cohort="GSA-admixed"


echo  "creating folders for VCFs"
mkdir -p $outputFolder
mkdir -p "${outputFolder}header"
mkdir -p "${outputFolder}filtered"
mkdir -p "${outputFolder}noID"
mkdir -p "${outputFolder}annotated"


for chr in {1..22}
do
    vcfChr=${chr}
    sbatch -J "$filterAnot_{chr}" -o "${outputFolder}/chr_${chr}.out" -e "${outputFolder}/chr_${chr}.err" \
        -v filterAndAnnotateMichigan.job \
        "$cohort" \
        $vcfChr \
        "${vcfFolder}/chr${chr}/chr${chr}.dose.vcf.gz" \
        "${outputFolder}" \
        $snpDB :q:
        #sbatch -J "${cohort}${chr}" -o "${outputFolder}/chr_${chr}.out" -e "${outputFolder}/chr_${chr}.err" \
        # 1- cohort
        # 2- chr
        # 3- input
        # 4- output
        # 4- snpDB reference
        sleep 0.4
done
