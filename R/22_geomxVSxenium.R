##### ##### ##### ##### ##### ##### ##### ##### ##### #####
# Parameter settings
##### ##### ##### ##### ##### ##### ##### ##### ##### #####
setwd(file.path(getwd(), "R"))
source("00_directorySettings.R") # load baseDir, xeniumDir
source("00_settings.R")

dataDir <- file.path(baseDir, "data")
outDir <- file.path(baseDir, "results")
dir.create(file.path(outDir, "Comparison"), showWarnings = FALSE)

section <- "RSC"
corMet <- "spearman"

##### ##### ##### ##### ##### ##### ##### ##### ##### #####
library(NanoStringNCTools)
library(GeoMxWorkflows)
library(stringr)
library(dplyr)
library(ggplot2)
library(ggrepel)

xenium <- readRDS(file.path(outDir, "Xenium", "04_ExprDf.RDS"))
xenium <- xenium[which(xenium$section == section),]
xExpr <- xenium[,c(1:(ncol(xenium)-3))]
xLabels <- xenium[,c((ncol(xenium)-2):ncol(xenium))]

geomx <- readRDS(file.path(dataDir, "GeoMx.RDS"))
gExpr <- assayData(geomx)$exprs # raw counts
gLabels <- pData(geomx)
colnames(gExpr) <- gLabels$Library

sampleExcludeIdx <- which(str_detect(colnames(gExpr), "TG1F")) # not on Xenium slide
gExpr <- gExpr[,-sampleExcludeIdx]
gLabels <- gLabels[-sampleExcludeIdx,]

genesInCommon <- intersect(colnames(xExpr), rownames(gExpr)) # 105 genes
xExpr <- xExpr[,which(colnames(xExpr) %in% genesInCommon)]
gExpr <- gExpr[which(rownames(gExpr) %in% genesInCommon),]
gExpr <- gExpr[colnames(xExpr),] # reorder

# across samples
coef <- c()
for (sample in names(modelCols[["Xenium"]])) {
        xCellIdx <- which(xLabels$region == sample)
        xSub <- xExpr[xCellIdx,]
        xTotal <- apply(xSub, 2, sum)

        gAoiIdx <- which(gLabels$Model == sample)
        gSub <- gExpr[,gAoiIdx]
        gTotal <- apply(gSub, 1, sum)

        if (all(names(xTotal) == names(gTotal))) {
                rho <- cor(xTotal, gTotal, method = corMet)
                coef <- c(coef, rho)
        }
}
names(coef) <- names(modelCols[["Xenium"]])
allCoef <- coef
allCoef

coef <- c()
segment <- "PV"
for (sample in names(modelCols[["Xenium"]])) {
        xCellIdx <- which(xLabels$region == sample)
        xSub <- xExpr[xCellIdx,]
        xTotal <- apply(xSub, 2, sum)

        gAoiIdx <- intersect(which(gLabels$Model == sample), which(gLabels$Segment == segment))
        gSub <- gExpr[,gAoiIdx]
        if (length(gAoiIdx) > 1) {
                gTotal <- apply(gSub, 1, sum)
        } else {
                gTotal <- gSub
        }

        if (all(names(xTotal) == names(gTotal))) {
                rho <- cor(xTotal, gTotal, method = corMet)
                coef <- c(coef, rho)
        }
}
names(coef) <- names(modelCols[["Xenium"]])
pvCoef <- coef
pvCoef

coef <- c()
segment <- "NeuN"
for (sample in names(modelCols[["Xenium"]])) {
        xCellIdx <- which(xLabels$region == sample)
        xSub <- xExpr[xCellIdx,]
        xTotal <- apply(xSub, 2, sum)

        gAoiIdx <- intersect(which(gLabels$Model == sample), which(gLabels$Segment == segment))
        gSub <- gExpr[,gAoiIdx]
        if (length(gAoiIdx) > 1) {
                gTotal <- apply(gSub, 1, sum)
        } else {
                gTotal <- gSub
        }

        if (all(names(xTotal) == names(gTotal))) {
                rho <- cor(xTotal, gTotal, method = corMet)
                coef <- c(coef, rho)
        }
}
names(coef) <- names(modelCols[["Xenium"]])
neunCoef <- coef
neunCoef

pdf(file.path(outDir, "Stats", "05_expr_corr.pdf"), width=5, height=5)
plot(pvCoef, neunCoef, pch=20, cex=2, col=modelCols[["Xenium"]], xlim=c(0, 0.65), ylim=c(0, 0.65), 
        xlab="Spearman(GeoMx-PV, Xenium)", ylab="Spearman(GeoMx-NeuN, Xenium)")
text(pvCoef, neunCoef, labels=names(pvCoef), pos=1, col=modelCols[["Xenium"]])
abline(coef = c(0,1), col="grey50", lty=2)
dev.off()

outDf <- data.frame(
        Sample = names(allCoef),
        All = as.numeric(allCoef),
        PV = as.numeric(pvCoef),
        NeuN = as.numeric(neunCoef)
)
write.table(outDf, file.path(outDir, "Stats", "05_expr_corr.txt"), row.names=F, col.names=T, quote=F, sep="\t")

# across genes
coef <- c()
for (gene in genesInCommon) {
        xGeneIdx <- which(colnames(xExpr) == gene)
        xDf <- data.frame(
                Sample = xLabels$region,
                Expression = xExpr[,xGeneIdx]
        )
        xDf$Sample <- factor(xDf$Sample, levels=names(modelCols[["Xenium"]]))
        xStats <- xDf %>% group_by(Sample) %>%  summarise(Total = sum(Expression))

        gGeneIdx <- which(rownames(gExpr) == gene)
        gDf <- data.frame(
                Sample = gLabels$Model,
                Expression = as.numeric(gExpr[gGeneIdx,])
        )
        gStats <- gDf %>% group_by(Sample) %>%  summarise(Total = sum(Expression))

        if (all(names(gDf$Sample) == names(xDf$Sample))) {
                rho <- cor(xStats$Total, gStats$Total, method = corMet)
                coef <- c(coef, rho)
        }
}
names(coef) <- genesInCommon
sampleCoef <- coef
sampleCoef

coef <- c()
segment <- "PV"
for (gene in genesInCommon) {
        xGeneIdx <- which(colnames(xExpr) == gene)
        xDf <- data.frame(
                Sample = xLabels$region,
                Expression = xExpr[,xGeneIdx]
        )
        xDf$Sample <- factor(xDf$Sample, levels=names(modelCols[["Xenium"]]))
        xStats <- xDf %>% group_by(Sample) %>%  summarise(Total = sum(Expression))

        gGeneIdx <- which(rownames(gExpr) == gene)
        gDf <- data.frame(
                Sample = gLabels$Model,
                Segment = gLabels$Segment,
                Expression = as.numeric(gExpr[gGeneIdx,])
        )
        gDf <- gDf[which(gDf$Segment == segment),]
        gStats <- gDf %>% group_by(Sample) %>%  summarise(Total = sum(Expression))

        if (all(names(gDf$Sample) == names(xDf$Sample))) {
                rho <- cor(xStats$Total, gStats$Total, method = corMet)
                coef <- c(coef, rho)
        }
}
names(coef) <- genesInCommon
samplePvCoef <- coef
samplePvCoef

coef <- c()
segment <- "NeuN"
for (gene in genesInCommon) {
        xGeneIdx <- which(colnames(xExpr) == gene)
        xDf <- data.frame(
                Sample = xLabels$region,
                Expression = xExpr[,xGeneIdx]
        )
        xDf$Sample <- factor(xDf$Sample, levels=names(modelCols[["Xenium"]]))
        xStats <- xDf %>% group_by(Sample) %>%  summarise(Total = sum(Expression))

        gGeneIdx <- which(rownames(gExpr) == gene)
        gDf <- data.frame(
                Sample = gLabels$Model,
                Segment = gLabels$Segment,
                Expression = as.numeric(gExpr[gGeneIdx,])
        )
        gDf <- gDf[which(gDf$Segment == segment),]
        gStats <- gDf %>% group_by(Sample) %>%  summarise(Total = sum(Expression))

        if (all(names(gDf$Sample) == names(xDf$Sample))) {
                rho <- cor(xStats$Total, gStats$Total, method = corMet)
                coef <- c(coef, rho)
        }
}
names(coef) <- genesInCommon
sampleNeuNCoef <- coef
sampleNeuNCoef

pdf(file.path(outDir, "Stats", "06_sample_cor.pdf"), width=5, height=5)
plot(samplePvCoef, sampleNeuNCoef, pch=20, cex=2, col="grey40", xlim=c(-0.8, 1), ylim=c(-0.8, 1.1), 
        xlab="Spearman(GeoMx-PV, Xenium)", ylab="Spearman(GeoMx-NeuN, Xenium)")
text(samplePvCoef, sampleNeuNCoef, labels=genesInCommon, pos=1, col="grey40")
abline(coef = c(0,1), col="grey50", lty=2)
dev.off()

outDf <- data.frame(
        Gene = names(sampleCoef),
        All = as.numeric(sampleCoef),
        PV = as.numeric(samplePvCoef),
        NeuN = as.numeric(sampleNeuNCoef)
)

write.table(outDf, file.path(outDir, "Stats", "06_sample_cor.txt"), row.names=F, col.names=T, quote=F, sep="\t")

pdf(file.path(outDir, "Stats", "06_sample_cor_all.pdf"), width=5, height=5)
hist(sampleCoef, prob = TRUE, main = "", xlab = "Spearman rho", xlim = c(-1, 1), breaks = seq(-1, 1, 0.1))
box()
dev.off()

q("no")
