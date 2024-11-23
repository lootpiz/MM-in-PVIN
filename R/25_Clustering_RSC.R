##### ##### ##### ##### ##### ##### ##### ##### ##### #####
# Parameter settings
##### ##### ##### ##### ##### ##### ##### ##### ##### #####
setwd(file.path(getwd(), "R"))
source("00_directorySettings.R") # load baseDir, xeniumDir
source("00_settings.R")

dataDir <- file.path(baseDir, "data")
outDir <- file.path(baseDir, "results", "Xenium")
dir.create(file.path(outDir, "Cluster"), showWarnings = TRUE)

section <- "RSC"
k <- 7 # 3 PV metagenes + 4 NeuN metagenes

##### ##### ##### ##### ##### ##### ##### ##### ##### #####
library(SingleCellExperiment)
library(arrow)

allExpr <- readRDS(file.path(outDir, "04_ExprDf.RDS"))
subExpr <- allExpr[which(allExpr$section == section),]
expr <- subExpr[,c(1:(ncol(subExpr)-3))]
annot <- subExpr[,c((ncol(subExpr)-2):ncol(subExpr))]

if (file.exists(file.path(outDir, "kMeans.RDS"))) {
        kmRes <- readRDS(file.path(outDir, "kMeans.RDS"))
} else {
        kmRes <- kmeans(expr, k, nstart = 25)
        saveRDS(kmRes, file.path(outDir, "kMeans.RDS"))
}

cls <- kmRes$cluster
if (all(names(cls) == rownames(cls))) {
        for (region in names(modelCols[["Xenium"]])) {
                subCls <- cls[which(annot$region == region)]
                outDf <- data.frame(
                        cell_id = names(subCls),
                        group = subCls
                )
                write.csv(outDf, file.path(outDir, "Cluster", paste0(region, "_RSC.csv")), row.names=F)

                
                slices <- table(subCls)
                pct <- round(slices/sum(slices)*100)
                labels <- paste0(slices, "\n(", pct, "%)")
                pdf(file.path(outDir, "Cluster", paste0(region, "_RSC.pdf")), width=5, height=5)
                pie(slices, labels = labels, col = kmeansCols, clockwise = TRUE)
                dev.off()

                pct <- round(slices/sum(slices)*100, 4)
                labels <- paste0(slices, "\n(", pct, "%)")
                names(labels) <- names(slices)
                write.csv(labels, file.path(outDir, "Cluster", paste0(region, "_RSC_summary.csv")), row.names=T)
        }
}

q("no")
