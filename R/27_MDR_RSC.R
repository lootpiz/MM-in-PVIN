##### ##### ##### ##### ##### ##### ##### ##### ##### #####
# Parameter settings
##### ##### ##### ##### ##### ##### ##### ##### ##### #####
setwd(file.path(getwd(), "R"))
source("00_directorySettings.R") # load baseDir, xeniumDir
source("00_settings.R")

dataDir <- file.path(baseDir, "data")
outDir <- file.path(baseDir, "results", "Xenium")
dir.create(file.path(outDir), showWarnings = FALSE)

section <- "RSC"

# UMAP parameters
custom.config <- umap::umap.defaults
custom.config$random_state <- 123567 # seed
custom.config$min_dist <- 0.1

##### ##### ##### ##### ##### ##### ##### ##### ##### #####
library(umap)
library(Rtsne)

kMeansLabels <- readRDS(file.path(outDir, "kMeans.RDS"))
kMeansLabels <- kMeansLabels$cluster

allExpr <- readRDS(file.path(outDir, "04_ExprDf.RDS"))
subExpr <- allExpr[which(allExpr$section == section),]
expr <- subExpr[,c(1:(ncol(subExpr)-3))]
labels <- subExpr[,c((ncol(subExpr)-2):ncol(subExpr))]
labels$kmeans <- kMeansLabels[rownames(labels)]

# Uniform Manifold Approximation and Projection (UMAP)
umapObj <- umap(expr, config = custom.config)

layout <- umapObj$layout
xylim <- range(layout)
xylim <- xylim + ((xylim[2]-xylim[1])*0.1)*c(-0.5, 0.5)

# Group
pdf(file.path(outDir, "UMAP_Group_RSC.pdf"), width=10, height=10)
plot(xylim, xylim, type="n", axes = F, frame = F, xlab="", ylab="")
points(layout[,1], layout[,2], col = groupCols[labels$group], cex = 0.6, pch = 20)
legend("topleft", legend = names(groupCols), fill = groupCols, inset = 0.03, bty = "n")
dev.off()

# Region
pdf(file.path(outDir, "UMAP_Reion_RSC.pdf"), width=10, height=10)
plot(xylim, xylim, type="n", axes = F, frame = F, xlab="", ylab="")
points(layout[,1], layout[,2], col = modelCols[["Xenium"]][labels$region], cex = 0.6, pch = 20)
legend("topleft", legend = names(modelCols[["Xenium"]]), fill = modelCols[["Xenium"]], inset = 0.03, bty = "n")
dev.off()

# kMeans
pdf(file.path(outDir, "UMAP_kMeans_RSC.pdf"), width=10, height=10)
plot(xylim, xylim, type="n", axes = F, frame = F, xlab="", ylab="")
points(layout[,1], layout[,2], col = kmeansCols[labels$kmeans], cex = 0.6, pch = 20)
legend("topleft", legend = names(kmeansCols), fill = kmeansCols, inset = 0.03, bty = "n")
dev.off()

# t-distributed Stochastic Neighbor Embedding (t-SNE)
tsneObj <- Rtsne(expr)

layout <- tsneObj$Y
xylim <- range(layout)
xylim <- xylim + ((xylim[2]-xylim[1])*0.1)*c(-0.5, 0.5)

# Group
pdf(file.path(outDir, "tSNE_Group_RSC.pdf"), width=10, height=10)
plot(xylim, xylim, type="n", axes = F, frame = F, xlab="", ylab="")
points(layout[,1], layout[,2], col = groupCols[labels$group], cex = 0.6, pch = 20)
legend("topleft", legend = names(groupCols), fill = groupCols, inset = 0.03, bty = "n")
dev.off()

# Region
pdf(file.path(outDir, "tSNE_Region_RSC.pdf"), width=10, height=10)
plot(xylim, xylim, type="n", axes = F, frame = F, xlab="", ylab="")
points(layout[,1], layout[,2], col =  modelCols[["Xenium"]][labels$region], cex = 0.6, pch = 20)
legend("topleft", legend = names(modelCols[["Xenium"]]), fill = modelCols[["Xenium"]], inset = 0.03, bty = "n")
dev.off()

# kMeans
pdf(file.path(outDir, "tSNE_kMeans_RSC.pdf"), width=10, height=10)
plot(xylim, xylim, type="n", axes = F, frame = F, xlab="", ylab="")
points(layout[,1], layout[,2], col = kmeansCols[labels$kmeans], cex = 0.6, pch = 20)
legend("topleft", legend = names(kmeansCols), fill = kmeansCols, inset = 0.03, bty = "n")
dev.off()

q("no")
