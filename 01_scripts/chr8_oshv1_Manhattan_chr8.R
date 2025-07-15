# Plot GEMMA results for de novo panel, AI2-imputed, and FI3-imputed
# B. Sutherland (VIU)
# 2025-07-15

# Load libraries
library("tidyr")

## Set working directory
current.path <- dirname(rstudioapi::getSourceEditorContext()$path)
current.path <- gsub(pattern = "\\/01_scripts", replacement = "", x = current.path) # take main directory
setwd(current.path)
rm(current.path)


# Set filenames
denovo_panel_results.FN <- "03_results/additional_file_s6_GEMMA_gwas_output_denovo_panel_survival.txt"
ai2_results.FN          <- "03_results/additional_file_s6_GEMMA_gwas_output_ai2_imputed_survival.txt"
fi2_results.FN          <- "03_results/additional_file_s6_GEMMA_gwas_output_fi3_imputed_survival.txt"

# User set variables
chr8_marker.pos <- 9719736


#### 01. Read in and prepare data ####
## Panel
denovo_results.df <- read.delim(file = denovo_panel_results.FN, header = T, sep = "\t")
head(denovo_results.df)
dim(denovo_results.df)

# Prepare chr and pos 
denovo_results.df <- separate(data = denovo_results.df, col = "rs", into = c("true.chr", "true.pos"), sep = "__", remove = F)
denovo_results.df$true.pos <- as.numeric(denovo_results.df$true.pos)
head(denovo_results.df)

# Prepare cutoff
denovo_cutoff.p <- 0.05 / nrow(denovo_results.df)


## AI2
ai2_results.df <- read.delim(file = ai2_results.FN, header = T, sep = "\t")
head(ai2_results.df)
tail(ai2_results.df) # note there are some NA rows
dim(ai2_results.df)

ai2_results.df <- ai2_results.df[!is.na(ai2_results.df$p_wald), ] 

dim(ai2_results.df)
table(ai2_results.df$panel_snp)

# Prepare chr and pos 
ai2_results.df <- separate(data = ai2_results.df, col = "rs", into = c("true.chr", "true.pos"), sep = "__", remove = F)
ai2_results.df$true.pos <- as.numeric(ai2_results.df$true.pos)
head(ai2_results.df)

# Prepare cutoff
ai2_cutoff.p <- 0.05 / nrow(ai2_results.df)


## FI3
fi3_results.df <- read.delim(file = fi2_results.FN, header = T, sep = "\t")
head(fi3_results.df)
tail(fi3_results.df)
dim(fi3_results.df)
table(fi3_results.df$panel_snp)

# Remove any NAs if present
fi3_results.df <- fi3_results.df[!is.na(fi3_results.df$p_wald), ] 
dim(fi3_results.df)

# Prepare chr and pos 
fi3_results.df <- separate(data = fi3_results.df, col = "rs", into = c("true.chr", "true.pos"), sep = "__", remove = F)
fi3_results.df$true.pos <- as.numeric(fi3_results.df$true.pos)
head(fi3_results.df)

# Prepare cutoff
fi3_cutoff.p <- 0.05 / nrow(fi3_results.df)


#### 02. Plot preparations ####
## denovo
# Keep only target chr
denovo_plot.df <- denovo_results.df[denovo_results.df$true.chr=="NC_047568.1",]
dim(denovo_plot.df) # 177 loci

# Max position
denovo_max_pos <- max(denovo_plot.df$true.pos)

# Keep only target chr, Ai2
ai2_plot.df <- ai2_results.df[ai2_results.df$true.chr=="NC_047568.1", ]
dim(ai2_plot.df) # 6192

# Max position
ai2_max_pos <- max(ai2_plot.df$true.pos)

# Keep only target chr, Fi3
fi3_plot.df <- fi3_results.df[fi3_results.df$true.chr=="NC_047568.1", ]
dim(fi3_plot.df)

# Max position
fi3_max_pos <- max(fi3_plot.df$true.pos)

all_max_pos <- max(c(denovo_max_pos, ai2_max_pos, fi3_max_pos))


#### 03. Plotting ####

pdf(file = "03_results/multipanel_Manhattan_plot_chr_NC_047568.1.pdf", width = 6, height = 9)
par(mfrow = c(3,1))

# denovo panel
plot(x = (denovo_plot.df$true.pos / 1000000), y = -log10(denovo_plot.df$p_wald)
     , pch = 16
     , xlim = c(0, all_max_pos/1000000)
     , xlab = "NC_047568.1 (Mbp)"
     , ylab = "-log10(p_wald)"
     , ylim = c(0, max(-log10(denovo_plot.df$p_wald)) + 1)
     , las = 1
     , col = "black"
)

abline(h = -log10(denovo_cutoff.p), lty = 2)
abline(v = chr8_marker.pos/1000000, lty = 3)

#abline(v = c(5, 15, 25, 35, 45, 55), lty = 4)


# AI2
plot(x = (ai2_plot.df$true.pos / 1000000), y = -log10(ai2_plot.df$p_wald)
     , pch = 16
     , xlim = c(0, all_max_pos/1000000)
     , xlab = "NC_047568.1 (Mbp)"
     , ylab = "-log10(p_wald)"
     , ylim = c(0, max(-log10(ai2_plot.df$p_wald)) + 1)
     , las = 1
     , col = "grey60"
       
)

points(x = (ai2_plot.df$true.pos[ai2_plot.df$panel_snp==TRUE] / 1000000)
       , y = -log10(ai2_plot.df$p_wald[ai2_plot.df$panel_snp==TRUE])
       , pch = 16 
       , col = "black"
       )

abline(h = -log10(ai2_cutoff.p), lty = 2)
abline(v = chr8_marker.pos/1000000, lty = 3)


# FI3
plot(x = (fi3_plot.df$true.pos / 1000000), y = -log10(fi3_plot.df$p_wald)
     , pch = 16
     , xlim = c(0, all_max_pos/1000000)
     , xlab = "NC_047568.1 (Mbp)"
     , ylab = "-log10(p_wald)"
     , ylim = c(0, max(-log10(fi3_plot.df$p_wald)) + 1)
     , las = 1
     , col = "grey60"
     
)


points(x = (fi3_plot.df$true.pos[fi3_plot.df$panel_snp==TRUE] / 1000000)
       , y = -log10(fi3_plot.df$p_wald[fi3_plot.df$panel_snp==TRUE])
       , pch = 16 
       , col = "black"
)


abline(h = -log10(fi3_cutoff.p), lty = 2)
abline(v = chr8_marker.pos/1000000, lty = 3)

dev.off()
