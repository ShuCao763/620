library(data.table)
sumstats <- fread("/Users/caoshu/WISC/620/final_project/Chronotype.txt.gz")

#1.1manhattan plot
library(qqman)
library(ggplot2)
gwas <- as.data.frame(sumstats, header=T)
colnames(gwas)[colnames(gwas) == "P_BOLT_LMM"] <- "P"

save_path <- "/Users/caoshu/WISC/620/final_project/"
png(file=paste0(save_path, "manhattan.png"), units="in", width=15, height=8, res=300)
manhattan(gwas, chr="CHR", bp="BP", p="P", snp="SNP", suggestiveline=F, genomewideline=-log10(5e-08))
dev.off()

#1.2qqplot
png(file=paste0(save_path, "qqplot.png"), units="in", width=8, height=8, res=300)
qq(gwas$P)
dev.off()

