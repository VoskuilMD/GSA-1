# Script for 


#!/bin/bash

for i in {1..22};
do
echo "#!/bin/bash" >> chr"$i"_vcftoplink.sh
echo "#SBATCH --job-name=chr"$i".vcftoplink" >> chr"$i"_vcftoplink.sh
echo "#SBATCH --mem 5gb" >> chr"$i"_vcftoplink.sh
echo "#SBATCH --time=4:00:00" >> chr"$i"_vcftoplink.sh
echo "#SBATCH --nodes 1" >> chr"$i"_vcftoplink.sh
echo "#SBATCH --output=chr"$i"_vcftoplink.out" >> chr"$i"_vcftoplink.sh
echo "#SBATCH --error=chr"$i"_vcftoplink.err" >> chr"$i"_vcftoplink.sh

echo "module load plink" >> chr"$i"_vcftoplink.sh
echo "module load GenotypeHarmonizer" >> chr"$i"_vcftoplink.sh
echo "module tabix " >> chr"$i"_vcftoplink.sh


# Set directory name for output 

echo "mkdir -p /groups/umcg-weersma/tmp04/Michiel/GSA-redo/imputation/european/results/european_GRCh37_noMAF_noR2_filters_3/plinkfiles/" >> chr"$i"_vcftoplink.sh

echo "GenotypeHarmonizer.sh \ " >> chr"$i"_vcftoplink.sh

# Set directories for input and ouput 

echo "-i /groups/umcg-weersma/tmp04/Michiel/GSA-redo/imputation/european/results/european_GRCh37_noMAF_noR2_filters_3/filtered/chr_${i} \ " >> chr"$i"_vcftoplink.sh
echo "-I VCFFOLDER \ " >> chr"$i"_vcftoplink.sh
echo "-o /groups/umcg-weersma/tmp04/Michiel/GSA-redo/imputation/european/results/european_GRCh37_noMAF_noR2_filters_2/plinkfiles/GSA_chr_${i} \ " >> chr"$i"_vcftoplink.sh
echo "-O PLINK_BED " >> chr"$i"_vcftoplink.sh;
done
