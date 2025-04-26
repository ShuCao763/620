library(biomaRt)
library(data.table)
library(ggplot2)
library(ggrepel)
library(qqman) 


## ==================False Discovery Rate (FDR)

twas <- fread("/Users/caoshu/WISC/620/final_project/twas/data/twas_brain_Frontal_Cortex.csv")
stopifnot("pvalue" %in% colnames(twas))

# Perform FDR correction using the Benjamini-Hochberg method
# This controls the expected proportion of false positives among declared significant results
twas$fdr <- p.adjust(twas$pvalue, method = "fdr")

# Filter significant genes with FDR < 0.05
sig_genes <- subset(twas, fdr < 0.05)

cat("Number of significant genes (FDR < 0.05):", nrow(sig_genes), "\n")

# Save the full TWAS result with FDR column
write.csv(twas, "/Users/caoshu/WISC/620/final_project/twas/processed_data/twas_brain_Frontal_Cortex_withFDR.csv", row.names = FALSE)

# Save the subset of significant genes (FDR < 0.05)
write.csv(sig_genes, "/Users/caoshu/WISC/620/final_project/twas/processed_data/twas_brain_Frontal_Cortex_sig_gene.csv", row.names = FALSE)


#===========================Query Ensembl using biomaRt================

twas <- fread("/Users/caoshu/WISC/620/final_project/twas/data/twas_brain_Frontal_Cortex.csv")
genes <- twas$gene_name
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

positions <- getBM(
  attributes = c("hgnc_symbol", "chromosome_name", "start_position"),
  filters = "hgnc_symbol",
  values = genes,
  mart = mart
)
# Clean and merge
setDT(positions)
setnames(positions, c("gene_name", "CHR", "BP"))

merged <- merge(twas, positions, by = "gene_name", all.x = TRUE)
merged <- merged[!is.na(CHR) & !is.na(BP)]
merged <- merged[grepl("^\\d+$", CHR)]
merged$CHR <- factor(merged$CHR, levels = as.character(1:22))
merged <- merged[order(as.numeric(as.character(merged$CHR)), BP)]
merged$pos_index <- 1:nrow(merged)

# Annotate significance threshold 
merged$logp <- -log10(merged$pvalue)
sig_thresh <- 0.05 / nrow(merged)
merged$significant <- merged$pvalue < sig_thresh

# Prepare chromosome axis breaks 
axisdf <- merged[, .(center = mean(pos_index)), by = CHR]


# ======================== Manhattan plot ===

p_manhattan <- ggplot(merged, aes(x = pos_index, y = logp, color = CHR)) +
  geom_point(size = 1.2) +
  geom_hline(yintercept = -log10(sig_thresh), color = "red", linetype = "dashed") +
  geom_text_repel(
    data = merged[significant == TRUE],
    aes(label = gene_name),
    size = 2.8,
    max.overlaps = 50
  ) +
  scale_x_continuous(label = axisdf$CHR, breaks = axisdf$center) +
  labs(title = "S-PrediXcan Manhattan Plot (Frontal Cortex)",
       x = "Chromosome", y = "-log10(P-value)") +
  theme_classic() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 0, vjust = 0.5))

# Save plot
ggsave("/Users/caoshu/WISC/620/final_project/twas/result/twas_manhattan_FrontalCortex.png",
       p_manhattan, width = 10, height = 5)


##  === QQ Plot ===

twas$pvalue[twas$pvalue == 0] <- 1e-300
# 计算 expected 和 observed
expected <- -log10((1:length(twas$pvalue)) / (length(twas$pvalue) + 1))
observed <- -log10(sort(twas$pvalue))
max_val <- max(c(expected, observed)) * 1.05  # 稍微放大一点防止挤在边缘

png("/Users/caoshu/WISC/620/final_project/twas/result/twas_qq_FrontalCortex.png", 
    width = 3000, height = 3000, res = 300)

plot(expected, observed,
     xlab = "Expected -log10(P)", 
     ylab = "Observed -log10(P)",
     xlim = c(0, max_val), 
     ylim = c(0, max_val),
     asp = 1,   # 关键：xy单位一致
     pch = 20, cex = 0.6,
     main = "S-PrediXcan QQ plot (Frontal Cortex)")

abline(0, 1, col = "red", lwd = 2)
dev.off()


##==========Volcano Plot（
twas <- fread("/Users/caoshu/WISC/620/final_project/twas/processed_data/twas_brain_Frontal_Cortex_withFDR.csv")
twas$logp <- -log10(twas$pvalue)
twas$significant <- twas$fdr < 0.05

# Only highly significant genes are annotated
twas$label <- ifelse(twas$fdr < 0.0001 & abs(twas$zscore) > 4, twas$gene_name, NA)

p <- ggplot(twas, aes(x = zscore, y = logp)) +
  geom_point(aes(color = significant), alpha = 0.6, size = 1.5) +
  geom_vline(xintercept = 0, linetype = "dotted", color = "black") +
  scale_color_manual(values = c("grey", "skyblue")) +
  geom_text_repel(aes(label = label), size = 3, max.overlaps = 50) +
  labs(title = "S-PrediXcan Volcano Plot (FrontalCortex)",
       x = "Z-score", y = "-log10(P-value)") +
  theme_minimal()

png("/Users/caoshu/WISC/620/final_project/twas/result/twas_volcano_FrontalCortex.png",width = 3000, height = 2000, res = 300)
print(p)
dev.off()



## PS: Top genes, Ranking and screening of the top 20 significant genes
twas <- fread("/Users/quincyqi/Documents/25SP/stat620/project/TWAS/results/ADHD_twas_BA9_withFDR.csv")
top_genes <- twas[fdr < 0.05][order(fdr)][1:20, .(gene_name, zscore, effect_size, fdr)]
fwrite(top_genes, "/Users/quincyqi/Documents/25SP/stat620/project/TWAS/results/Top20_ADHD_genes.csv")

