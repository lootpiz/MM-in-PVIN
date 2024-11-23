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
library(HDF5Array)
library(SingleCellExperiment)
library(scMerge)
library(stringr)

outputFolders <- list.files(path = xeniumDir, pattern = "output")
featuresAnnot <- read.delim(file.path(xeniumDir, outputFolders[1], "cell_feature_matrix", "features.tsv.gz"), header = FALSE, stringsAsFactors = FALSE)

arExprL <- readRDS(file.path(xeniumDir, "Results", "01_analysisReady_expression.RDS"))

selectCellFiles <- list.files(
        path = file.path(xeniumDir, "ROIs"), 
        pattern = "cellsstats", 
        full.names = TRUE
)

cellIdsL <- lapply(seq_along(selectCellFiles), function(idx) {
        selectCellFile <- selectCellFiles[[idx]]
        selectCell <- read.csv(selectCellFile, header = T, comment.char = "#")
        return(selectCell$Cell.ID)
})
names(cellIdsL) <- sapply(str_split(basename(selectCellFiles), "__cells"), "[[", 1)

geneAnnot <- featuresAnnot[,2]
names(geneAnnot) <- featuresAnnot[,1]

subExprL <- sapply(seq_along(cellIdsL), function(idx) {
        selectionName <- names(cellIdsL)[idx]
        regionName <- unlist(str_split(selectionName, "__"))[1]
        
        cellIds <- cellIdsL[[idx]]
        expr <- arExprL[[regionName]]
        rownames(expr) <- geneAnnot[rownames(expr)]
        subExpr <- expr[,which(colnames(expr) %in% cellIds)]
        
        geneSymbols <- as.vector(rownames(subExpr))
        if (startsWith(regionName, "T")) {
                condi = 1 # TG
        } else {
                condi = 0 # WT
        }
        
        meta <- data.frame(
                individual = regionName,
                cell_id = colnames(subExpr),
                condition = condi
        )
        rownames(meta) <- colnames(subExpr)
        meta$condition <- factor(meta$condition, levels=c(0, 1))
        obj <- SingleCellExperiment(
                assays = list(counts = as.matrix(subExpr)), 
                colData = meta, 
                rowData = geneSymbols
        )
        return(obj)
})
names(subExprL) <- names(cellIdsL)
saveRDS(subExprL, file.path(outDir, "02_ROIs_inSCE.RDS"))

roiL <- list()
for (section in names(sectionCols)) {
        roi <- sce_cbind(subExprL[which(str_detect(names(subExprL), section))], method = "intersect", exprs = c("counts"), colData_names = TRUE, cut_off_batch = 0)
        roiL[[section]] <- roi
        # saveRDS(roi, file.path(outDir, paste0("03_", section, ".RDS")))
}
names(roiL) <- names(sectionCols)
saveRDS(roiL, file.path(outDir, "03_ROIs.RDS")) # presto input

allExprL <- lapply(seq_along(cellIdsL), function(idx) {
        selectionName <- names(cellIdsL)[idx]
        regionName <- unlist(str_split(selectionName, "__"))[1]
        sectionName <- unlist(str_split(selectionName, "__"))[2]
        
        cellIds <- cellIdsL[[idx]]
        expr <- arExprL[[regionName]]
        rownames(expr) <- geneAnnot[rownames(expr)]
        subExpr <- expr[,which(colnames(expr) %in% cellIds)]
        
        geneSymbols <- as.vector(rownames(subExpr))
        
        if (startsWith(regionName, "T")) {
                group = "TG"
        } else {
                group = "WT"
        }

        expr <- data.frame(
                t(subExpr),
                group = group,
                region = regionName,
                section = sectionName
        )
        return(expr)
})
allExpr <- Reduce(rbind, allExprL)
saveRDS(allExpr, file.path(outDir, "04_ExprDf.RDS")) # UMAP, tSNE input

q("no")
