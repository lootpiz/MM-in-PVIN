##### ##### ##### ##### ##### ##### ##### ##### ##### #####
# Parameter settings
##### ##### ##### ##### ##### ##### ##### ##### ##### #####
setwd(file.path(getwd(), "R"))
source("00_directorySettings.R") # load baseDir <- "/WORKING_DIR"
source("00_settings.R")

dataDir <- file.path(baseDir, "data")
nmfDir <- file.path(baseDir, "data", "NMF")
outDir <- file.path(baseDir, "results", "SOTK")

selectRanks <- c(3, 4) # PV, NeuN
names(selectRanks) <- segments

howMany <- 100 # Top

##### ##### ##### ##### ##### ##### ##### ##### ##### #####
library(NanoStringNCTools)
library(GeoMxWorkflows)
library(GeomxTools)
library(stringr)
library(UpSetR)
library(NMF)

for (segment in segments) {
        selectRank <- selectRanks[segment]

        # Cor (Expression, cNMF Usage) > 0.2 genes based on consensus (H)
        expr <- readRDS(file.path(dataDir, paste0("GeoMx.", segment, ".RDS")))
        NMFres <- readRDS(file.path(nmfDir, paste0(segment, ".RDS")))
        
        H <- coef(get(as.character(selectRank), NMFres$fit))
        rownames(H) <- paste0("Metagene_", c(1:selectRank))

        corGenesL <- list()
        for (kdx in c(1:selectRank)) {
                consensusValues <- as.matrix(H[kdx,])
                colnames(consensusValues) <- c("consensus")

                for (geneIdx in c(1:nrow(expr))) {
                        geneName <- rownames(expr)[geneIdx]
                        geneExpr <- as.matrix(expr[geneIdx,])
                        colnames(geneExpr) <- c("Expression")

                        dat <- merge(consensusValues, geneExpr, by="row.names", all=F)
                        rho <- cor(dat[,2], dat[,3], method="spearman")

                        if (rho > 0.2) {
                                corGenesL[[paste0("Metagene_", kdx)]] <- c(corGenesL[[paste0("Metagene_", kdx)]], geneName)
                        }
                }
        }

        # Top 'howMany' genes based on basis (W)
        W <- basis(get(as.character(selectRank), NMFres$fit))
        colnames(W) <- paste0("Metagene_", c(1:selectRank))

        topGenesL <- list()
        for (kdx in c(1:selectRank)) {
                basisValues <- W[,kdx]
                basisValues <- basisValues[order(basisValues, decreasing=T)]
                topGenesL[[paste0("Metagene_", kdx)]] <- names(basisValues)[1:howMany]
        }

        highContributingGenesL <- list()
        for (kdx in c(1:selectRank)) {
                setA <- corGenesL[[paste0("Metagene_", kdx)]]
                setB <- topGenesL[[paste0("Metagene_", kdx)]]
                HCGs <- intersect(setA, setB)

                highContributingGenesL[[paste0("Metagene_", kdx)]] <- HCGs
                write.table(sort(HCGs), file.path(outDir, segment, paste0("13_HCGs_Metagene_", kdx, "_", length(HCGs), ".txt")), row.names=F, col.names=F, quote=F)
        }

        uPlot <- upset(fromList(highContributingGenesL),
                order.by = c("freq", "degree"),
                nsets = selectRank,
                decreasing = c(T, T),
                number.angles = 0,
                point.size = 16,
                text.scale = c(4, 4, 4, 2, 4, 4),
                line.size = 2,
                mainbar.y.label = "Intersection",
                sets.x.label = "Number of genes",
                keep.order = TRUE
        )
        pdf(file.path(outDir, segment, paste0("14_HCGs_", length(unique(unlist(highContributingGenesL))), ".pdf")), width = 16, height = 12)
        print(uPlot)
        dev.off()
}

q("no")
