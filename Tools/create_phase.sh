#!/bin/bash

for i in {1..22};
do
echo "#!/bin/bash" >> chr"$i"_phase.sh
echo "#SBATCH --job-name=chr"$i".phase" >> chr"$i"_phase.sh
echo "#SBATCH --mem 8gb" >> chr"$i"_phase.sh
echo "#SBATCH --time=9:00:00" >> chr"$i"_phase.sh
echo "#SBATCH --nodes 1" >> chr"$i"_phase.sh
echo "#SBATCH --output=chr"$i"_phase.out" >> chr"$i"_phase.sh
echo "#SBATCH --error=chr"$i"_phase.err" >> chr"$i"_phase.sh
echo "module load Java" >> chr"$i"_phase.sh
echo "module load beagle/27Jul16.86a-Java-1.8.0_45" >> chr"$i"_phase.sh
echo "java -jar -Xmx8G /apps/software/beagle/27Jul16.86a-Java-1.8.0_45/beagle.27Jul16.86a.jar gt=chr_"$i"_filtered.vcf.gz map=plink.chr"$i".GRCh37.map out=chr"$i"_phased chrom="$i" impute=false"  >> chr"$i"_phase.sh;
done

#download b37 hapmap genetic map here: http://bochet.gcc.biostat.washington.edu/beagle/genetic_maps/plink.GRCh37.map.zip
