##### ##### ##### ##### ##### ##### ##### ##### ##### #####
# Parameter settings
##### ##### ##### ##### ##### ##### ##### ##### ##### #####
setwd(file.path(getwd(), "R"))
source("00_directorySettings.R") # load baseDir <- "/WORKING_DIR"
source("00_settings.R")

dataDir <- file.path(baseDir, "data")
nmfDir <- file.path(baseDir, "data", "NMF")
outDir <- file.path(baseDir, "results", "Xenium")
dir.create(file.path(outDir), showWarnings = FALSE)

selectRanks <- c(3, 4) # PV, NeuN
names(selectRanks) <- segments

##### ##### ##### ##### ##### ##### ##### ##### ##### #####
library(SingleCellExperiment)
library(InSituType)
library(tidyverse)
library(reshape2)
library(ggplot2)
library(stringr)
library(NMF)

MatrixL <- lapply(seq_along(segments), function(idx) {
        segment <- segments[idx]
        NMFres <- readRDS(file.path(nmfDir, paste0(segment, ".RDS")))
        selectRank <- selectRanks[segment]

        w <- basis(get(as.character(selectRank), NMFres$fit))
        colSum <- apply(w, 2, sum)
        wNorm <- w/colSum
        colnames(wNorm) <- paste0(segment, "_M", c(1:selectRank))

        return(wNorm)
})
names(MatrixL) <- segments

buff <- merge(MatrixL[[1]], MatrixL[[2]], by="row.names", all=T)
rownames(buff) <- buff[,1]
wMatrix <- buff[,c(2:ncol(buff))]
wMatrix[is.na(wMatrix)] <- 0 # InSituType 2.0 does not take NA anymore

xeniumL <- readRDS(file.path(dataDir, "Xenium.Expr.RDS")) # Expression profiles
xeniumROIL <- readRDS(file.path(dataDir, "Xenium.ROI.RDS")) # 02_ROIs_inSCE.RDS, subset

deconL <- lapply(seq_along(xeniumROIL), function(idx) {
        selectionName <- names(xeniumROIL)[idx]
        regionName <- unlist(str_split(selectionName, "__"))[1]
        
        expr <- xeniumL[[regionName]]
        roiObj <- xeniumROIL[[idx]]
        roi <- assays(roiObj)$counts

        subExpr <- expr[, which(colnames(expr) %in% colnames(roi))]
        negProbesIdx <- which(str_detect(rownames(subExpr), "^NegControlProbe"))
        negCounts <- subExpr[negProbesIdx,]

        sup <- insitutypeML(x = t(roi),
                neg = Matrix::rowMeans(t(negCounts)),
                reference_profiles = as.matrix(wMatrix)
        )

        outDf <- data.frame(
                cell_id = names(sup$clust),
                group = as.matrix(sup$clust)
        )
        write.csv(outDf, file.path(outDir, "Cluster", paste0(regionName, "_InSituType.csv")), row.names=F)

        # make the flightpath plot
        fp <- flightpath_plot(flightpath_result = NULL, insitutype_result = sup, col = metagenesCols[sup$clust], showclusterconfidence = FALSE)
        pdf(file.path(outDir, paste0("Flightpath_", selectionName, ".pdf")), width=4, height=4)
        print(fp)
        dev.off()

        return(sup)
})
names(deconL) <- names(xeniumROIL)
saveRDS(deconL, file.path(outDir, "InSituType.RDS"))

q("no")
