#2.1 ldsc regression(#bin!bash)
cd ldsc-2.0.1
python ./munge_sumstats.py \
--p P_BOLT_LMM \
--N 400000 \
--snp SNP \
--a1 ALLELE1 \
--a2 ALLELE0 \
--sumstats /Users/caoshu/WISC/620/final_project/data/Chronotype.txt.gz \
--merge-alleles /Users/caoshu/WISC/620/ldsc_inputs/w_hm3.snplist \
--out /Users/caoshu/WISC/620/final_project/data/munged_Chronotype

python ./ldsc.py \
--h2 /Users/caoshu/WISC/620/final_project/data/munged_Chronotype.sumstats.gz \
--ref-ld-chr /Users/caoshu/WISC/620/ldsc_inputs/for_h2/eur_w_ld_chr/ \
--w-ld-chr /Users/caoshu/WISC/620/ldsc_inputs/for_enrichment/weights/weights.hm3_noMHC. \
--out /Users/caoshu/WISC/620/final_project/result/ldsc_Chronotype 

#2.2 Partition trait heritability by tissue annotations 
python ./ldsc.py \
--h2 /Users/caoshu/WISC/620/final_project/data/munged_Chronotype.sumstats.gz \
--ref-ld-chr /Users/caoshu/WISC/620/ldsc_inputs/for_enrichment/Baseline/baseline.,/Users/caoshu/WISC/620/ldsc_inputs/for_enrichment/GenoSkylinePlus/GSplus_Tier3_1KGphase3. \
--w-ld-chr /Users/caoshu/WISC/620/ldsc_inputs/for_enrichment/weights/weights.hm3_noMHC. \
--overlap-annot \
--frqfile-chr /Users/caoshu/WISC/620/ldsc_inputs/for_enrichment/genotype/1000G.EUR.QC. \
--out /Users/caoshu/WISC/620/final_project/result/ldsc_enrich_Chronotype

