if (!requireNamespace("remotes")) install.packages("remotes")
remotes::install_github("MRCIEU/TwoSampleMR")

install.packages("MendelianRandomization")

library(data.table)
library(dplyr)
library(TwoSampleMR)
library(MendelianRandomization)


setwd("/Users/caoshu/WISC/620/final_project/mr")

#==============================
# Exposure: Chronotype (GWAS summary stats)
#==============================
Chronotype <- fread("/Users/caoshu/WISC/620/final_project/data/Chronotype.txt.gz")
Chronotype = Chronotype %>%
  mutate(
    beta_exp = BETA, 
    p = P_BOLT_LMM,
    A1 = ALLELE1,
    A2 = ALLELE0,
    se_exp = SE,
  ) %>%
  select(SNP, CHR, BP, A1, A2, beta_exp, se_exp, p)

# Strong instruments (P < 5e-8)
Chronotype_strong <- Chronotype %>% filter(p < 5e-8)

# Format for clumping
Chronotype_formatted <- Chronotype_strong %>%
  rename(
    pval.exposure = p,
    beta.exposure = beta_exp,
    se.exposure = se_exp,
    effect_allele.exposure = A1,
    other_allele.exposure = A2
  ) %>%
  mutate(
    chr_name = CHR,
    chrom_start = BP,
    id.exposure = "Chronotype",
    samplesize.exposure = 449734  #总的样本量
  )

#==============================
# LD Clumping (locally)
#==============================
Chronotype_clumped <- clump_data(
Chronotype_formatted,
  clump_kb = 10000,
  clump_r2 = 0.001,
  bfile = "/Users/caoshu/WISC/620/final_project/mr/data/bfile/1kg_hm3_QCed_noM",
  plink_bin = "/Users/caoshu/WISC/620/plink_mac_20241022/plink" #plink路径
)
Chronotype_clumped$exposure <- "Chronotype"  #只要对Chronotype进行clumping

fwrite(
  Chronotype_clumped,
  file = "/Users/caoshu/WISC/620/final_project/mr/Chronotype_clumped.txt",
  sep = "\t",
)
#==============================
# Outcome: EA
#==============================
ea <- fread("/Users/caoshu/WISC/620/Data2/EA.gz") 
ea = ea %>%
  mutate(
    #Beta = log(OR) ,
    outcome = "SCZ",
    id.outcome = "SCZ" 
  ) %>%
  rename(
    beta.outcome = Beta,
    se.outcome = SE,
    effect_allele.outcome = A1,
    other_allele.outcome = A2
  ) %>%
  select(SNP, beta.outcome, se.outcome, effect_allele.outcome, other_allele.outcome, outcome, id.outcome)


#==============================
# Harmonise exposure and outcome
#==============================
harmonised <- harmonise_data(
  exposure_dat = Chronotype_clumped,
  outcome_dat = ea
   action = 2
)

# MRInput 
mr_obj <- mr_input(
  bx = harmonised$beta.exposure,
  bxse = harmonised$se.exposure,
  by = harmonised$beta.outcome,
  byse = harmonised$se.outcome,
  snps = harmonised$SNP
)

# IVW 
res_ivw <- mr_ivw(mr_obj)

# run Egger regression
res_egger <- mr_egger(mr_obj)

# IVW result
sink("/Users/caoshu/WISC/620/final_project/mr/result/ivw_result_Chronotype_to_EA_strongIV.txt")  #改
cat("Chronotype ➜ EA - IVW Result\n\n")  #改
cat("Number of SNPs:   ", res_ivw@SNPs, "\n")
cat("Estimate:         ", res_ivw@Estimate, "\n")
cat("Std. Error:       ", res_ivw@StdError, "\n")
cat("95% CI:           (", res_ivw@CILower, ", ", res_ivw@CIUpper, ")\n")
cat("P-value:          ", res_ivw@Pvalue, "\n")
cat("Heterogeneity Q:  ", res_ivw@Heter.Stat[1], "\n")
cat("Heterogeneity p:  ", res_ivw@Heter.Stat[2], "\n")
cat("F-statistic:      ", res_ivw@Fstat, "\n")
sink()

# Egger result
sink("/Users/caoshu/WISC/620/final_project/mr/result/egger_result_Chronotype_to_EA_strongIV.txt")
cat("Chronotype ➜ EA - Egger Regression Result\n\n")
cat("Number of SNPs:   ", res_egger@SNPs, "\n")
cat("Causal Estimate:  ", res_egger@Estimate, "\n")
cat("Std. Error:       ", res_egger@StdError.Est, "\n")
cat("P-value:          ", res_egger@Pvalue.Est, "\n")
cat("Intercept (pleiotropy): ", res_egger@Intercept, "\n")
cat("Intercept SE:     ", res_egger@StdError.Int, "\n")
cat("Intercept P-value:", res_egger@Pvalue.Int, "\n")
sink()


# IVW scatter plots
top_df <- harmonised %>%
  arrange(se.exposure) %>%
  slice(1:500)  # Top 500 SNPs

top_mr_obj <- mr_input(
  bx = top_df$beta.exposure,
  bxse = top_df$se.exposure,
  by = top_df$beta.outcome,
  byse = top_df$se.outcome,
  snps = top_df$SNP
)

p <- mr_plot(top_mr_obj, error = TRUE, line = "ivw", interactive = FALSE)
png("/Users/caoshu/WISC/620/final_project/mr/result/ivw_plot_top500_Chronotype_to_EA_labeled.png", width = 800, height = 600)
p + 
  ggplot2::labs(
    x = "SNP effect on Chronotype", 
    y = "SNP effect on EA",
    title = "IVW Scatter Plot "
  ) +
  ggplot2::theme_minimal()
dev.off()