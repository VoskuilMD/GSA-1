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

# We will manually add the NOD2:rs2066845 (flip), rs1505992(nonflip), rs7746082 (nonflip) here, since this is an ambigious snp but I have checked for strands
echo 'rs2066845' >> Lui_Sommeren_loci_non_amb.txt
echo 'rs1505992' >> Lui_Sommeren_loci_non_amb.txt
echo 'rs7746082' >> Lui_Sommeren_loci_non_amb.txt
plink --bfile GIEQ.Lui.loci --extract Lui_Sommeren_loci_non_amb.txt --make-bed --out GIEQ.Lui.loci.non.amb --allow-no-sex

# For these 216 I have checked manually whether the A1 .assoc allele corresponds with the A1 .bim allele. Please see excel: /GIEQ/GIEQ.snps.flips.xlsx
# I have used the .bim file as a lead. I.e. when A1.bim was A2.assoc, I converted the beta into the other direction. (similar to 1/OR).
# I have uploaded this info into the file GIEQ.snps.flips

awk '{print $1,$7,$9}' GIEQ.snps.flips > tmp
sed -e 's/ [ ]*/\t/g' tmp > snp.effectsizes
rm tmp

#Deprecated
# Now we are going to convert the binary plink files to an "additive (0/1/2) component file"
#plink --bfile GIEQ.Lui.loci.non.amb --recode A --out GIEQ.Lui.loci.non.amb.switch --allow-no-sex 

# A1 alleles are now counted and in the GIEQ.Lui.loci.non.amb.raw file.

# This is what we actually do!!
#Alternatively, we can also keep the orignal beta's and or's from the Lui paper. In this case, we have to count the .assoc A1 alleles in the plink files
# For this we use the recode-allele file
# Upload the file: recode.allels.according.to.base
plink --bfile  GIEQ.Lui.loci.non.amb --recode A --recode-allele recode.allels.according.to.base --out GIEQ.Lui.loci.non.amb.original

# as you can see, the snp.effectsizes files has corresponding BETA's/OR's for both the .switch and the .original genetic data. 
# Please check whether PRS scores will be the same eventually -> they are not, so stick with .original. 

# IN Rstudio:

#switch.raw = read.table("~/Documents/Werk/Promotie/GIEQ/GIEQ.Lui.loci.non.amb.switch.raw", header = F)
original.raw = read.table("~/Documents/Werk/Promotie/GIEQ/GIEQ.Lui.loci.non.amb.original.raw", header = F)
#switch.raw = switch.raw[,c(-1,-3,-4,-5,-6)]
original.raw = original.raw[,c(-1,-3,-4,-5,-6)]
#switch = as.data.frame(t(switch.raw))
original = as.data.frame(t(original.raw))

#names(switch) <- as.matrix(switch[1, ])
#switch <- switch[-1, ]
#switch[] <- lapply(switch, function(x) type.convert(as.character(x)))

names(original) <- as.matrix(original[1, ])
original <- original[-1, ]
original[] <- lapply(original, function(x) type.convert(as.character(x)))



