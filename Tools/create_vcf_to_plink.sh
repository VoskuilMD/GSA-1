# Script for converting imputed VCFs to plink files

#!/bin/bash

for i in {1..22};
do
echo "#!/bin/bash" >> chr"$i"_vcftoplink.sh
echo "#SBATCH --job-name=chr"$i".vcftoplink" >> chr"$i"_vcftoplink.sh
echo "#SBATCH --mem 5gb" >> chr"$i"_vcftoplink.sh
echo "#SBATCH --time=9:00:00" >> chr"$i"_vcftoplink.sh
echo "#SBATCH --nodes 1" >> chr"$i"_vcftoplink.sh
echo "#SBATCH --output=chr"$i"_vcftoplink.out" >> chr"$i"_vcftoplink.sh
echo "#SBATCH --error=chr"$i"_vcftoplink.err" >> chr"$i"_vcftoplink.sh

echo "module load plink" >> chr"$i"_vcftoplink.sh
echo "module load GenotypeHarmonizer" >> chr"$i"_vcftoplink.sh
echo "module load tabix " >> chr"$i"_vcftoplink.sh


# Set directory name for, input, output and r2 filter 

echo "mkdir -p /groups/umcg-weersma/tmp04/Michiel/GSA-redo/imputation/european/results/european_GRCh37_noMAF_noR2_filters_3/plinkfiles/" >> chr"$i"_vcftoplink.sh

echo "GenotypeHarmonizer.sh -i /groups/umcg-weersma/tmp04/Michiel/GSA-redo/imputation/european/results/european_GRCh37_noMAF_noR2_filters_3/annotated/chr_${i} -I VCFFOLDER -o /groups/umcg-weersma/tmp04/Michiel/GSA-redo/imputation/european/results/european_GRCh37_noMAF_noR2_filters_3/plinkfiles/GSA_chr_${i} -O PLINK_BED -mrf 0.3" >> chr"$i"_vcftoplink.sh;
done
