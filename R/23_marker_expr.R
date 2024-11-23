##### ##### ##### ##### ##### ##### ##### ##### ##### #####
# Parameter settings
##### ##### ##### ##### ##### ##### ##### ##### ##### #####
setwd(file.path(getwd(), "R"))
source("00_directorySettings.R") # load baseDir <- "/WORKING_DIR"

outDir <- file.path(baseDir, "results")

section <- "RSC"

##### ##### ##### ##### ##### ##### ##### ##### ##### #####
library(stringr)
library(UpSetR)
library(ggplot2)

allExpr <- readRDS(file.path(outDir, "Xenium", "04_ExprDf.RDS"))
subExpr <- allExpr[which(allExpr$section == section),]
expr <- subExpr[,c(1:(ncol(subExpr)-3))]

kMeansLabels <- readRDS(file.path(outDir, "Xenium", "kMeans.RDS"))
kMeansLabels <- kMeansLabels$cluster

geneListFiles <- list.files(file.path(outDir, "GeneLists"), pattern = "txt", full.names = T)

genesL <- lapply(seq_along(geneListFiles), function(idx) {
        genes <- read.table(geneListFiles[idx], header = F, stringsAsFactors = F)
        return(genes$V1)
})
names(genesL) <- sapply(str_split(basename(geneListFiles), "\\."), "[[", 1)

dodge <- position_dodge(width = 0.3)
if (all(names(kMeansLabels) == rownames(expr))) { # check the order of cell IDs
        for (listName in c("PV_M1", "PV_M2", "PV_M3", "NeuN_M1", "NeuN_M2", "NeuN_M3", "NeuN_M4", "UP_in_PV")) {
                markers <- intersect(genesL[["Xenium_panel"]], genesL[[listName]])
                for (gene in markers) {
                        df <- data.frame(
                                Cell = names(kMeansLabels), # cell ID
                                Cluster = kMeansLabels, # kMeans cluster
                                Expression = as.numeric(expr[,which(colnames(expr) == gene)])
                        )
                        df$Cluster <- factor(df$Cluster)

                        violin <- ggplot(data = df, aes(x = Cluster, y = asinh(Expression), fill = Cluster)) +
                                geom_violin(position = dodge, size = 0) +
                                geom_boxplot(width = 0.1, position = dodge, fill="white") +
                                scale_fill_manual(values = kmeansCols) +
                                labs(
                                        title = listName,
                                        subtitle = gene,
                                        x = "kMeans", 
                                        y = "Expression"
                                ) +
                                theme_bw() +
                                theme(
                                        axis.line = element_line(colour = "black"),
                                        panel.grid.major = element_blank(),
                                        panel.grid.minor = element_blank(),
                                        panel.border = element_blank(),
                                        panel.background = element_blank(),
                                        axis.text.x = element_text(angle = 0, vjust = 0, hjust = 0.5),
                                        legend.position = "none", 
                                        text = element_text(size = 12)
                                )
                        
                        pdf(file.path(outDir, "Xenium", paste0("Expression_", listName, "_", gene, ".pdf")), width = 7, height = 7)
                        print(violin)
                        dev.off()
                }
        }
}

q("no")
