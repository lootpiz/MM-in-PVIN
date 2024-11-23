##### ##### ##### ##### ##### ##### ##### ##### ##### #####
# Parameter settings
##### ##### ##### ##### ##### ##### ##### ##### ##### #####
setwd(file.path(getwd(), "R"))
source("00_directorySettings.R") # load baseDir <- "/WORKING_DIR"; xeniumDir <- "/XENIUM_DIR"
source("00_settings.R")

dataDir <- file.path(baseDir, "data")
outDir <- file.path(baseDir, "results", "Xenium")
dir.create(file.path(outDir, "InSituType"), showWarnings = FALSE)

##### ##### ##### ##### ##### ##### ##### ##### ##### #####
library(SingleCellExperiment)
library(InSituType)
library(stringr)
library(jpeg)
library(ggplot2)
library(grid)
library(arrow)

deconL <- readRDS(file.path(outDir, "InSituType.RDS"))

dummy <- sapply(seq_along(deconL), function(idx) {
        roiName <- names(deconL)[idx]
        region <- unlist(str_split(roiName, "__"))[1]
        selection <- unlist(str_split(roiName, "__"))[2]

        regionFolder <- file.path(xeniumDir, paste0("output-XETG00116__0010348__", region, "__20240427__002555"))
        
        # Load ROI info
        selectionFile <- file.path(xeniumDir, "ROIs", paste0(roiName, "__coordinates.csv"))
        selectRegion <- read.table(selectionFile, header=T, sep=",")
        xMin <- min(selectRegion$X)
        xMax <- max(selectRegion$X)
        yMin <- min(selectRegion$Y)
        yMax <- max(selectRegion$Y)

        # Load ROI img
        imgFile <- file.path(xeniumDir, "ROIs", paste0(roiName, "_Sel.jpg"))
        img <- readJPEG(imgFile) # fetch image dimension

        # Load QC passed cell-/gene-id
        subExprL <- readRDS(file.path(xeniumDir, "Results", "02_ROIs_inSCE.RDS"))
        subExpr <- subExprL[[roiName]]
        expr <- assays(subExpr)$counts
        geneIds <- rownames(expr)
        cellIds <- colnames(expr)

        # Load Cell Info
        cells <- arrow::read_parquet(file.path(regionFolder, "cells.parquet"), as_data_frame = TRUE)
        subCells <- cells[which(cells$cell_id %in% cellIds),]
        cellDf <- data.frame(
                CellID = subCells$cell_id,
                X = subCells$x_centroid,
                Y = subCells$y_centroid
        )
        
        # Load Tx info
        tx <- arrow::read_parquet(file.path(regionFolder, "transcripts.parquet"), as_data_frame = TRUE)

        # Cluster info
        sup <- deconL[[roiName]]
        clusterInfo <- sup$clust

        cellDf$Cluster <- clusterInfo[cellDf$CellID]
        cellDf$Segment <- sapply(str_split(cellDf$Cluster, "_"), "[[", 1)
        cellDf$Segment <- factor(cellDf$Segment, levels=segments)
        cellDf$Metagene <- sapply(str_split(cellDf$Cluster, "_"), "[[", 2)
        cellDf$Metagene <- factor(cellDf$Metagene)
        
        # Plot
        segRoi <- ggplot(cellDf, aes(x = X, y = Y, color = Segment, shape = Segment)) + 
                geom_point(size = 3) +
                scale_y_reverse() +
                theme_bw() +
                labs(
                        caption = roiName
                ) +
                theme(
                        axis.line = element_blank(),
                        panel.grid.major = element_blank(),
                        panel.grid.minor = element_blank(),
                        panel.border = element_blank(),
                        panel.background = element_blank(),
                        axis.text.x = element_blank(),
                        axis.title.x = element_blank(),
                        axis.text.y = element_blank(),
                        axis.title.y = element_blank(),
                        axis.ticks = element_blank(),
                        # legend.position = "none",
                        text = element_text(size = 20)
                )
        png(file.path(outDir, "InSituType", paste0(roiName, ".png")), width = dim(img)[2], height = dim(img)[1])
        print(segRoi)
        dev.off()

        for (cluster in sort(unique(cellDf$Cluster))) {
                subCellDf <- cellDf[which(cellDf$Cluster == cluster),]
                
                metaRoi <- ggplot(cellDf, aes(x = X, y = Y)) + 
                        geom_point(shape = 16, color = "grey50") +
                        geom_point(data = subCellDf, aes(x = X, y = Y), size = 2, shape = 15, color = "red") +
                        scale_y_reverse() +
                        theme_bw() +
                        labs(
                                caption = paste0(roiName, "\n", cluster)
                        ) +
                        theme(
                                axis.line = element_blank(),
                                panel.grid.major = element_blank(),
                                panel.grid.minor = element_blank(),
                                panel.border = element_blank(),
                                panel.background = element_blank(),
                                axis.text.x = element_blank(),
                                axis.title.x = element_blank(),
                                axis.text.y = element_blank(),
                                axis.title.y = element_blank(),
                                axis.ticks = element_blank(),
                                legend.position = "none",
                                text = element_text(size = 20)
                        )
                png(file.path(outDir, "InSituType", paste0(roiName, "_", cluster, ".png")), width = dim(img)[2], height = dim(img)[1])
                print(metaRoi)
                dev.off()
        }

        for (gene in genes) {
                subTx <- tx[which(tx$feature_name == gene),] # subset a given gene
                if (nrow(subTx) > 0) {
                        dir.create(file.path(outDir, gene), showWarnings = FALSE)

                        subTx <- subTx[which(subTx$cell_id %in% cellIds),] # subset cells in ROI

                        geneDf <- data.frame(
                                Type = subTx$feature_name,
                                X = subTx$x_location,
                                Y = subTx$y_location
                        )

                        txRoi <- ggplot(cellDf, aes(x = X, y = Y)) + 
                                geom_point(shape = 16, color = "grey50") +
                                geom_point(data = subCellDf, aes(x = X, y = Y), size = 2, shape = 15, color = "red") +
                                geom_point(data = geneDf, aes(x = X, y = Y), size = 3, shape = 4, color = "darkgreen") +
                                scale_y_reverse() +
                                theme_bw() +
                                labs(
                                        caption = paste0(roiName, "\n", cluster, "\n", gene)
                                ) +
                                theme(
                                        axis.line = element_blank(),
                                        panel.grid.major = element_blank(),
                                        panel.grid.minor = element_blank(),
                                        panel.border = element_blank(),
                                        panel.background = element_blank(),
                                        axis.text.x = element_blank(),
                                        axis.title.x = element_blank(),
                                        axis.text.y = element_blank(),
                                        axis.title.y = element_blank(),
                                        axis.ticks = element_blank(),
                                        legend.position = "none",
                                        text = element_text(size = 20)
                                )

                        png(file.path(outDir, gene, paste0(roiName, "_", cluster, "_", gene, ".png")), width = dim(img)[2], height = dim(img)[1])
                        print(txRoi)
                        dev.off()
                } else {
                        dir.create(file.path(outDir, paste0(gene, "_noHit")), showWarnings = FALSE)
                }
        }

        return(roiName)
})

q("no")
