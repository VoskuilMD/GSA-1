GENETIC_MAP_DIR=/apps/data/www.shapeit.fr/genetic_map_b37/
REF_DIR=/apps/data/gonl/imputationReference_phased_b37/
BED_DIR=/groups/umcg-weersma/tmp04/Michiel/GSA-redo/imputation/european/results/european_GRCh37_maf0001/plinkfiles_biallelic/
OUTPUT_DIR=/groups/umcg-weersma/tmp04/Michiel/GSA-redo/phasing/shapeit/results/
TEMP_DIR=/groups/umcg-weersma/tmp04/Michiel/GSA-redo/phasing/shapeit/

mkdir -p ${OUTPUT_DIR}
mkdir -p ${TEMP_DIR}

for chr in {1,2,6,7,10,11,12,13,14,15,16,19,22}
do
	echo '#!/bin/bash' >> ${TEMP_DIR}/chr${chr}.sh
	echo '#SBATCH --job-name=phase_beds'${chr} >> ${TEMP_DIR}/chr${chr}.sh
	echo '#SBATCH --output='${TEMP_DIR}/phase_beds_chr_$chr.out >> ${TEMP_DIR}/chr${chr}.sh
	echo '#SBATCH --error='${TEMP_DIR}/phase_beds_chr_$chr.err >> ${TEMP_DIR}/chr${chr}.sh
	echo '#SBATCH --time=23:59:00' >> ${TEMP_DIR}/chr${chr}.sh
	echo '#SBATCH --cpus-per-task=8' >> ${TEMP_DIR}/chr${chr}.sh
	echo '#SBATCH --mem-per-cpu=10gb' >> ${TEMP_DIR}/chr${chr}.sh
	echo '#SBATCH --nodes=1' >> ${TEMP_DIR}/chr${chr}.sh
	echo '#SBATCH --open-mode=truncate' >> ${TEMP_DIR}/chr${chr}.sh
	echo '#SBATCH --export=NONE' >> ${TEMP_DIR}/chr${chr}.sh
	echo ml shapeit >> ${TEMP_DIR}/chr${chr}.sh
	echo shapeit \
        -B ${BED_DIR}/GSA_chr_${chr}_biallelic.bed \
        -M ${GENETIC_MAP_DIR}genetic_map_chr${chr}_combined_b37.txt \
        --input-ref ${REF_DIR}chr${chr}.hap.gz \
        -O ${OUTPUT_DIR}/chr_${chr} \
        --output-log ${OUTPUT_DIR}/chr_${chr}.log \
        --thread 8 >> ${TEMP_DIR}/chr${chr}.sh

    sbatch ${TEMP_DIR}/chr${chr}.sh
    
done
