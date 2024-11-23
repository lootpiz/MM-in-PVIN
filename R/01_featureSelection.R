##### ##### ##### ##### ##### ##### ##### ##### ##### #####
# Parameter settings
##### ##### ##### ##### ##### ##### ##### ##### ##### #####
setwd(file.path(getwd(), "R"))
source("00_directorySettings.R") # load baseDir <- "/WORKING_DIR"

dataDir <- file.path(baseDir, "data")

howManyFeatures <- 2000
coefVariationThreshold <- 0.01

##### ##### ##### ##### ##### ##### ##### ##### ##### #####
library(NanoStringNCTools)
library(GeoMxWorkflows)
library(GeomxTools)
library(stringr)

library(MetaMx) # private

geomx <- readRDS(file.path(dataDir, "GeoMx.RDS"))
annot <- pData(geomx)
expr <- assayData(geomx)$vst
colnames(expr) <- annot$Library

pvIdx <- which(annot$Segment == "PV")
neIdx <- which(annot$Segment == "NeuN")

exprL <- list(PV = expr[,pvIdx], NeuN = expr[,neIdx])

featuresL <- lapply(seq_along(exprL), function(idx) {
        segmentName <- names(exprL)[idx]
        return(getMVGs(exprL[[idx]], coefVar = coefVariationThreshold, no = howManyFeatures))
})
names(featuresL) <- names(exprL)
saveRDS(featuresL, file.path(dataDir, "FeaturesL.RDS"))

dummy <- sapply(seq_along(exprL), function(idx) {
        segmentName <- names(exprL)[idx]
        exprs <- exprL[[idx]]
        features <- featuresL[[idx]]
        outMat <- exprs[which(rownames(exprs) %in% features),]
        saveRDS(outMat, file.path(dataDir, paste0("GeoMx.", segmentName, ".RDS")))
        return(dim(outMat))
})

q("no")
