# Extract Lui loci from 2015 trans-ethnic meta-analysis summ stats:

cd /groups/umcg-weersma/tmp04/Michiel/PRSice/subphenotypes/BASE_FILES/raw/IBD
grep -f Lui_Sommeren_loci.txt IBD_trans_ethnic_association_summ_stats_b37.txt > IBD_trans_ethnic_association_summ_stats_b37_Lui_loci.txt
head -n1 IBD_trans_ethnic_association_summ_stats_b37.txt > IBD_trans_ethnic_association_summ_stats_b37_Lui_loci.assoc
cat IBD_trans_ethnic_association_summ_stats_b37_Lui_loci.txt >> IBD_trans_ethnic_association_summ_stats_b37_Lui_loci.assoc
sed -e 's/ [ ]*/\t/g' IBD_trans_ethnic_association_summ_stats_b37_Lui_loci.assoc > Lui_loci_b37.assoc 

# Extract Lui loci from GSA - GIEQ data
cd /groups/umcg-weersma/tmp04/Michiel/GIEQ/PRS.manually
ml plink
plink --bfile ../GIEQ --extract Lui_Sommeren_loci.txt --make-bed --out GIEQ.Lui.loci --allow-no-sex

# Align strand via PRSice
awk '{print $1,$2,$3,$4,$5,$10,$11,$12,$13,$14,$15}' Lui_loci_b37.assoc > BASE_Lui_loci_b37.tmp 
sed -e 's/ [ ]*/\t/g' BASE_Lui_loci_b37.tmp > BASE_Lui_loci_b37.assoc 
rm BASE_Lui_loci_b37.tmp

ml R
Rscript /groups/umcg-weersma/tmp04/Michiel/PRSice/PRSice.R --dir . --prsice /groups/umcg-weersma/tmp04/Michiel/PRSice/PRSice_linux \
--base /groups/umcg-weersma/tmp04/Michiel/GIEQ/PRS.manually/BASE_Lui_loci_b37.assoc --target GIEQ.Lui.loci \
--stat beta_EUR --beta --A1 A1_effect --se se_EUR --A2 A2_other --pvalue P_EUR --snp SNP --chr Chr --bp Pos \
--out TEMP --pheno-file GIEQ.diag.pheno --binary-target T \
--pheno-col pheno --lower 0.01 --upper 0.1 --interval 0.01 \
--all-score --print-snp  --no-regress

# We need to remove ambigious variants
Fortunately, PRSice does this automagically. See log file: "10 ambiguous variant(s) excluded 216 variant(s) included"
"Number of variant(s) after clumping : 213"

# We will only include the 213 clumped non-ambigious variants
awk '{print $2}' TEMP.snp > Lui_Sommeren_loci_non_amb.txt 
plink --bfile GIEQ.Lui.loci --extract Lui_Sommeren_loci_non_amb.txt --make-bed --out GIEQ.Lui.loci.non.amb --allow-no-sex





