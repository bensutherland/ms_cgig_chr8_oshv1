# Plot multiple fastman results together
# B. Sutherland (VIU)
# 2024-10-22
library(fastman)

## Set working directory
current.path <- dirname(rstudioapi::getSourceEditorContext()$path)
current.path <- gsub(pattern = "\\/01_scripts", replacement = "", x = current.path) # take main directory
setwd(current.path)
rm(current.path)


# Set filenames
denovo_panel_data <- "./../../CHR8_amp/ms_cgig_chr8_oshv1/03_results_denovo_v.0.4/output/post_Manhattan_plotting.RData"
ai2_data          <- "../impute_workflow_v.0.8_AF/07_GWAS/ai2_imputed_survival_state_pheno/post-plot.RData"
fi3_data          <- "../impute_workflow_v.0.8_AF/07_GWAS/fi3_imputed_survival_state_pheno/post-plot.RData"

highlight_snps.FN <- "../impute_workflow_v.0.8_AF/07_GWAS/denovo_snp_ids.txt"

# Output plot
output.FN <- "03_results/multipanel_Manhattan_plot.pdf"

#### 01. Load data ####
## Load denovo data
load(denovo_panel_data)

# Save gemma output data
gemma_gwas_denovo_panel <- gemma_gwas
colnames(gemma_gwas_denovo_panel)[grep(pattern = "pos.true", x = colnames(gemma_gwas_denovo_panel))] <- "pos"
gemma_gwas_denovo_panel <- gemma_gwas_denovo_panel[grep(pattern = "^Chr", x = gemma_gwas_denovo_panel$chr), ] # only keep loci on chr
rm(gemma_gwas)


## Load ai2 data
load(ai2_data)

# Save gemma output data
gemma_gwas_ai2 <- gemma_gwas
rm(gemma_gwas)

## Load fi3 data
load(fi3_data)

# Save gemma output data
gemma_gwas_fi3 <- gemma_gwas
rm(gemma_gwas)

# Read in the SNPs to highlight (single column with marker name, e.g., 'chr__position')
highlight_snps.FN <- "../impute_workflow_v.0.8_AF/07_GWAS/denovo_snp_ids.txt" # need to reassign, as this was written over by loading previous data
highlight_snps.df <- read.delim(file = highlight_snps.FN, header = F, sep = "\t")
head(highlight_snps.df)

# Limit highlight vector 
# Create highlight SNP vector, making sure that they are actually in the gemma output, and accounting
highlight.vec <- highlight_snps.df$V1
length(highlight.vec)

highlight_ai2.vec <- highlight.vec[highlight.vec %in% gemma_gwas_ai2$rs]
length(highlight_ai2.vec)

highlight_fi3.vec <- highlight.vec[highlight.vec %in% gemma_gwas_fi3$rs]
length(highlight_fi3.vec)


#### 02. Plot individual Manhattan plots ####
# denovo panel
pdf(file = "03_results/Manhattan_plot_denovo_panel.pdf", width = 9.5, 4)
par(mar = c(5,6,2,2) +0.1, mgp = c(3,1,0))
fastman(m = gemma_gwas_denovo_panel
        , chr = "chr", bp = "pos", p = "p_wald", snp = "rs"
        , genomewideline = -log10(0.05/nrow(gemma_gwas_denovo_panel))
        , suggestiveline = NULL
        , cex = 0.7, cex.lab = 1, cex.axis = 1
        #, ylim = c(0,10)
        , col = "Set2"
        , annotateHighlight = F
        , annotationAngle = 55
        , annotationCol = "black"
        , maxP = plot_maxP
        #, annotatePval = -log10(0.05/nrow(gemma_gwas_denovo_panel))
        #, annotateTop = T
        
        
)
dev.off()

# ai2
pdf(file = "03_results/Manhattan_plot_ai2.pdf", width = 9.5, 4)
par(mar = c(5,6,2,2) +0.1, mgp = c(3,1,0))

fastman(m = gemma_gwas_ai2
        , chr = "chr", bp = "pos", p = "p_wald", snp = "rs"
        , genomewideline = -log10(0.05/nrow(gemma_gwas_ai2))
        , suggestiveline = NULL
        , cex = 0.7, cex.lab = 1, cex.axis = 1
        #, ylim = c(0,10)
        , col = "Set2"
        , annotateHighlight = F
        , annotationAngle = 55
        , annotationCol = "black"
        , highlight = highlight_ai2.vec
        , maxP = plot_maxP
        # , annotatePval = 0.05/nrow(gemma_gwas_ai2)
        # , annotateTop = T
        
)
dev.off()


# fi3
pdf(file = "03_results/Manhattan_plot_fi3.pdf", width = 9.5, 5)
par(mar = c(5,6,2,2) +0.1, mgp = c(3,1,0))

fastman(m = gemma_gwas_fi3
        , chr = "chr", bp = "pos", p = "p_wald", snp = "rs"
        , genomewideline = -log10(0.05/nrow(gemma_gwas_fi3))
        , suggestiveline = NULL
        , cex = 0.7, cex.lab = 1, cex.axis = 1
        #, ylim = c(0,10)
        , col = "Set2"
        , annotateHighlight = F
        , annotationAngle = 55
        , annotationCol = "black"
        , highlight = highlight_fi3.vec
        , maxP = plot_maxP
        
)

dev.off()



#### 02. Plot multi-panel ####
pdf(file = output.FN, height = 7, width =  6)

par(mfrow = c(3,1), mar = c(5,6,2,2) +0.1, mgp = c(3,1,0))

# denovo panel
par(mar = c(5,6,2,2) +0.1, mgp = c(3,1,0))
fastman(m = gemma_gwas_denovo_panel
        , chr = "chr", bp = "pos", p = "p_wald", snp = "rs"
        , genomewideline = -log10(0.05/nrow(gemma_gwas_denovo_panel))
        , suggestiveline = NULL
        , cex = 0.7, cex.lab = 1, cex.axis = 1
        #, ylim = c(0,10)
        , col = "Set2"
        , annotateHighlight = F
        , annotationAngle = 55
        , annotationCol = "black"
        , maxP = plot_maxP
        #, annotatePval = -log10(0.05/nrow(gemma_gwas_denovo_panel))
        #, annotateTop = T
        
        
)

mtext(text = "A", side = 2, line = 3
      , at = (max(-log10(gemma_gwas_denovo_panel$p_wald)) + ( 0.1  * max(-log10(gemma_gwas_denovo_panel$p_wald))))
      , las = 1)


# ai2
par(mar = c(5,6,2,2) +0.1, mgp = c(3,1,0))
fastman(m = gemma_gwas_ai2
        , chr = "chr", bp = "pos", p = "p_wald", snp = "rs"
        , genomewideline = -log10(0.05/nrow(gemma_gwas_ai2))
        , suggestiveline = NULL
        , cex = 0.7, cex.lab = 1, cex.axis = 1
        #, ylim = c(0,10)
        , col = "Set2"
        , annotateHighlight = F
        , annotationAngle = 55
        , annotationCol = "black"
        , highlight = highlight_ai2.vec
        , maxP = plot_maxP
        # , annotatePval = 0.05/nrow(gemma_gwas_ai2)
        # , annotateTop = T
        
)

mtext(text = "B", side = 2, line = 3
      , at = (max(-log10(gemma_gwas_ai2$p_wald)) + ( 0.1  * max(-log10(gemma_gwas_ai2$p_wald))))
      , las = 1)


# fi3
par(mar = c(5,6,2,2) +0.1, mgp = c(3,1,0))
fastman(m = gemma_gwas_fi3
        , chr = "chr", bp = "pos", p = "p_wald", snp = "rs"
        , genomewideline = -log10(0.05/nrow(gemma_gwas_fi3))
        , suggestiveline = NULL
        , cex = 0.7, cex.lab = 1, cex.axis = 1
        #, ylim = c(0,10)
        , col = "Set2"
        , annotateHighlight = F
        , annotationAngle = 55
        , annotationCol = "black"
        , highlight = highlight_fi3.vec
        , maxP = plot_maxP
        
)

mtext(text = "C", side = 2, line = 3
      , at = (max(-log10(gemma_gwas_fi3$p_wald)) + ( 0.1  * max(-log10(gemma_gwas_fi3$p_wald))))
      , las = 1)


dev.off()

