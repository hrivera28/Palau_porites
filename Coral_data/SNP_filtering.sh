#Using the TotalRawSNPs output file from dDocent
# Remove sites with where more than 50 % of samples are missing data
vcftools --vcf TotalRawSNPs.vcf --max-missing 0.5 --recode --recode-INFO-all --out mis05
# First remove PIPA samples
vcftools --vcf mis05.recode.vcf --remove pipa_to_rmv --recode --recode-INFO-all --out all_noPIPA
# 1,952,537 sites
vcftools --vcf all_noPIPA.recode.vcf --max-missing 0.7 --mac 3 --minQ 20 --minDP 8 --mac 0.1 --recode --recode-INFO-all --out mis07_mac3_Q20_DP8_maf1
# 295,000 Sites

# Remove individuals with more than 70% missing data
/docent/plut_ref/filt_script_puritz/filter_missing_ind.sh mis07_mac3_Q20_DP8_maf1.recode.vcf mis07_mac3_Q20_DP8_lowDPindvDrop
# removed 16 individual

docent/plut_ref/filt_script_puritz/dDocent_filters mis07_mac3_Q20_DP8_lowDPindvDrop.recode.vcf mis07_mac3_Q20_DP8_lowDPindvDrop_docentfilt

vcfallelicprimitives mis07_mac3_Q20_DP8_lowDPindvDrop_docentfilt.FIL.recode.vcf --keep-info --keep-geno > mis07_mac3_Q20_DP8_lowDPindvDrop_docentfilt_prim.vcf

# Remove indels
vcftools --vcf mis07_mac3_Q20_DP8_lowDPindvDrop_docentfilt_prim.vcf --remove-indels --recode --recode-INFO-all --out mis07_mac3_Q20_DP8_lowDPindvDrop_docentfilt_SNPS
# 155195 out of a possible 167699

# Biallele SNPs only
vcftools --vcf mis07_mac3_Q20_DP8_lowDPindvDrop_docentfilt_SNPS.recode.vcf --min-alleles 2 --max-alleles 2 --recode --recode-INFO-all --out mis07_mac3_Q20_DP8_lowDPindvDrop_docentfilt_SNPS_biall
#137854 out of a possible 155195

# Dropping samples with more than 50% missing data at sites
/docent/plut_ref/filt_script_puritz/filter_missing_ind.sh mis07_mac3_Q20_DP8_lowDPindvDrop_docentfilt_SNPS_biall.recode.vcf mis07_mac3_Q20_DP8_lowDPindvDrop2_docentfilt_SNPS_biall
# 146 of 149 (at 50% cutoff) 137854 sites

# Filtering sites to minor allele frequency of at least 0.05
vcftools --vcf mis07_mac3_Q20_DP8_lowDPindvDrop2_docentfilt_SNPS_biall.recode.vcf --max-missing 0.8 --maf 0.05 --recode --recode-INFO-all --out mis08_mac3_Q20_DP8_lowDPindvDrop2_docentfilt_SNPS_biall_maf05
# After filtering, kept 45103 out of a possible 137854 Sites

# Only one snp per rad tag (seq are 150 bp long, this thins to select only SNPps that are at least 150 bp apart):
vcftools --vcf mis08_mac3_Q20_DP8_lowDPindvDrop2_docentfilt_SNPS_biall_maf05.recode.vcf --thin 150 --recode --recode-INFO-all --out Allsamps_maf05_final_snps
# kept 12761 out of a possible 45103 Sites
