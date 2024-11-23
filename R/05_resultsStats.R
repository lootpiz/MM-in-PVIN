##### ##### ##### ##### ##### ##### ##### ##### ##### #####
# Parameter settings
##### ##### ##### ##### ##### ##### ##### ##### ##### #####
setwd(file.path(getwd(), "R"))
source("00_directorySettings.R") # load baseDir <- "/WORKING_DIR"
source("00_settings.R")

dataDir <- file.path(baseDir, "data")
outDir <- file.path(baseDir, "results", "NMF")

maxRank <- 10

##### ##### ##### ##### ##### ##### ##### ##### ##### #####
library(SummarizedExperiment)
library(NanoStringNCTools)
library(GeoMxWorkflows)
library(RColorBrewer)
library(GeomxTools)
library(circlize)
library(stringr)
library(NMF)

dir.create(file.path(outDir), showWarnings = FALSE)

obj <- readRDS(file.path(dataDir, "GeoMx.RDS"))
annotAll <- pData(obj)

rankCol <- c(
        "#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#A6761D", 
        "#B3CDE3", "#E6AB02", "#386CB0", "#FED9A6", "#FA6BA5"
)
names(rankCol) <- c(1:maxRank)

for (segment in segments) {
        NMFres <- readRDS(file.path(dataDir, "NMF", paste0(segment, ".RDS")))
        NMFres_rand <- readRDS(file.path(dataDir, "NMF", paste0(segment, "_rand.RDS")))

        # clustering results/stats
        pdf(file = file.path(outDir, paste0(segment, "_ranksurvey.pdf")), width = 10, height = 7)
        plot(NMFres, NMFres_rand)
        dev.off()
        
        annot <- annotAll[which(annotAll$Segment == segment),]

        # Heat map column annotation
        phenotype <- data.frame(
                Group   = factor(annot$Group,   levels=names(groupCols)),
                Model   = factor(annot$Model,   levels=names(modelCols[["GeoMx"]])),
                Segment = factor(annot$Segment, levels=segments)
        )
        rownames(phenotype) <- annot$Library

        annColors <- list(
                Group     = groupCols,
                Model     = modelCols[["GeoMx"]],
                Segment   = segmentCols,
                basis     = rankCol, 
                consensus = rankCol
        )
        
        # consensusmap
        pdf(file = file.path(outDir, paste0(segment, "_consensusmatrix.pdf")), width = 10, height = 10, onefile = TRUE)
        for (k in c(2:maxRank)) {
                consensusmap(
                        get(as.character(k), NMFres$fit),
                        tracks = c("basis", "consensus", "silhouette"),
                        annCol = phenotype,
                        annColors = annColors,
                        Rowv = TRUE,
                        Colv = TRUE,
                        col = colorRampPalette(c("#F8FCC6", "#BC312A"))(100),
                        main = paste0("rank = ", k)
                )
        }
        dev.off()
        
        # metagenes
        pdf(file = file.path(outDir, paste0(segment, "_H.pdf")), width = 10, height = 8)
        for (k in c(2:maxRank)) {
                coefmap(get(as.character(k), NMFres$fit),
                        Colv = TRUE,
                        annCol = phenotype,
                        annColors = annColors,
                        col = colorRampPalette(c("#F8FCC6", "#BC312A"))(100),
                        main = paste0("rank = ", k)
                )
        }
        dev.off()

        # basis
        pdf(file = file.path(outDir, paste0(segment, "_W.pdf")), width = 10, height = 40)
        for (k in c(2:maxRank)) {
                basismap(get(as.character(k), NMFres$fit),
                        Rowv = TRUE,
                        col = colorRampPalette(c("#F8FCC6", "#BC312A"))(100),
                        main = paste0("rank = ", k)
                )
        }
        dev.off()
}

q("no")
