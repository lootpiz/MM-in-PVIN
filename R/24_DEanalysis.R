##### ##### ##### ##### ##### ##### ##### ##### ##### #####
# Parameter settings
##### ##### ##### ##### ##### ##### ##### ##### ##### #####
setwd(file.path(getwd(), "R"))
source("00_directorySettings.R") # load baseDir, xeniumDir
source("00_settings.R")

dataDir <- file.path(baseDir, "data")
outDir <- file.path(baseDir, "results", "Xenium")
dir.create(file.path(outDir), showWarnings = FALSE)

##### ##### ##### ##### ##### ##### ##### ##### ##### #####
library(SingleCellExperiment)
library(presto)
library(ggplot2)
library(ggrepel)

roiL <- readRDS(file.path(outDir, "03_ROIs.RDS"))

for (section in names(sectionCols)) {
        roi <- roiL[[section]]
        pvalues <- wilcoxauc(roi, 'condition')
        pvalues <- pvalues[which(pvalues$group == 1),] # TG vs. WT

        df <- as.data.frame(pvalues)
        df$direction <- "none"
        df$direction[intersect(which(df$logFC < (wilcoxonDegFcThre * -1)), which(df$padj < wilcoxonDegPvalThre))] <- "negative"
        df$direction[intersect(which(df$logFC > wilcoxonDegFcThre), which(df$padj < wilcoxonDegPvalThre))] <- "positive"
        df$direction <- factor(df$direction, levels=c("positive", "none", "negative"))
        write.table(df, file.path(outDir, paste0("Wilcoxon_", section, ".txt")), row.names=F, col.names=T, quote=F, sep="\t")

        df$label <- df$feature
        df$label[which(df$direction == "none")] <- NA
        if (length(which(df$pval == 0)) > 0) {
                df$pval[which(df$pval == 0)] <- min(df$pval[df$pval > 0]) # pseudo value for p=0
        }

        volcanoPlot <- ggplot(data = df, aes(x = logFC, y = -log10(padj), col = direction, label = label)) +
                        geom_point(cex = 3) +
                        theme_minimal() +
                        geom_text_repel(max.overlaps = 10) +
                        geom_vline(xintercept = c(-wilcoxonDegFcThre, wilcoxonDegFcThre), col = "grey40", lty = 2) +
                        geom_hline(yintercept = -log10(wilcoxonDegPvalThre), col = "grey40", lty = 2) +
                        scale_color_manual(
                                name = "DEG", values = c("#AE3C32", "#F9FBCB", "#46639C"),
                                labels = c("Up-reg. in TG", "Not significant", "Down-reg. in TG")
                        ) +
                        labs(
                                title = paste0("TG vs. WT in ", section),
                                x = "Log2 fold-change", 
                                y = "-Log10(Adjusted P-value)"
                        ) +
                        theme(
                                axis.line = element_line(colour = "black"),
                                panel.grid.major = element_blank(),
                                panel.grid.minor = element_blank(),
                                panel.border = element_blank(),
                                panel.background = element_blank(),
                                text = element_text(size = 12)
                        )

        pdf(file.path(outDir, paste0("Wilcoxon_", section, ".pdf")), width=8, height=6)
        print(volcanoPlot)
        dev.off()
}

q("no")
