##### ##### ##### ##### ##### ##### ##### ##### ##### #####
# Parameter settings
##### ##### ##### ##### ##### ##### ##### ##### ##### #####
setwd(file.path(getwd(), "R"))
source("00_directorySettings.R") # load baseDir <- "/WORKING_DIR"

outDir <- file.path(baseDir, "results", "GeneLists")

##### ##### ##### ##### ##### ##### ##### ##### ##### #####
library(stringr)
library(UpSetR)

geneListFiles <- list.files(file.path(outDir), pattern = "txt", full.names = T)

genesL <- lapply(seq_along(geneListFiles), function(idx) {
        genes <- read.table(geneListFiles[idx], header = F, stringsAsFactors = F)
        return(genes$V1)
})
names(genesL) <- sapply(str_split(basename(geneListFiles), "\\."), "[[", 1)

uPlot <- upset(fromList(genesL[c("DOWN_in_PV", "UP_in_PV", "Xenium_panel")]),
        order.by = c("freq", "degree"),
        decreasing = c(T, T),
        number.angles = 0,
        point.size = 20,
        text.scale = c(4, 4, 4, 2, 4, 4),
        line.size = 2,
        mainbar.y.label = "Intersection",
        sets.x.label = "Number of genes",
        keep.order = TRUE
)

pdf(file.path(outDir, "01_DEGs__Panel.pdf"), width = 18, height = 12)
print(uPlot)
dev.off()

uPlot <- upset(fromList(genesL[c("PV_M1", "PV_M2", "PV_M3", "Xenium_panel")]),
        order.by = c("freq", "degree"),
        decreasing = c(T, T),
        number.angles = 0,
        point.size = 20,
        text.scale = c(4, 4, 4, 2, 4, 4),
        line.size = 2,
        mainbar.y.label = "Intersection",
        sets.x.label = "Number of genes",
        keep.order = TRUE
)

pdf(file.path(outDir, "02_PV_MAG__Panel.pdf"), width = 18, height = 12)
print(uPlot)
dev.off()

uPlot <- upset(fromList(genesL[c("NeuN_M1", "NeuN_M2", "NeuN_M3", "NeuN_M4", "Xenium_panel")]),
        order.by = c("freq", "degree"),
        decreasing = c(T, T),
        number.angles = 0,
        point.size = 20,
        text.scale = c(4, 4, 4, 2, 4, 4),
        line.size = 2,
        mainbar.y.label = "Intersection",
        sets.x.label = "Number of genes",
        keep.order = TRUE
)

pdf(file.path(outDir, "03_NeuN_MAG__Panel"), width = 18, height = 12)
print(uPlot)
dev.off()

uPlot <- upset(fromList(genesL[c("DOWN_in_PV", "UP_in_PV", "PV_M1", "PV_M2", "PV_M3", "Xenium_panel")]),
        order.by = c("freq", "degree"),
        decreasing = c(T, T),
        number.angles = 0,
        point.size = 20,
        text.scale = c(4, 4, 4, 2, 4, 4),
        line.size = 2,
        nsets = 6,
        mainbar.y.label = "Intersection",
        sets.x.label = "Number of genes",
        keep.order = TRUE
)

pdf(file.path(outDir, "04_PV_MAG__DEGs__Panel"), width = 18, height = 12)
print(uPlot)
dev.off()

uPlot <- upset(fromList(genesL[c("DOWN_in_PV", "UP_in_PV", "NeuN_M1", "NeuN_M2", "NeuN_M3", "NeuN_M4", "Xenium_panel")]),
        order.by = c("freq", "degree"),
        decreasing = c(T, T),
        number.angles = 0,
        point.size = 20,
        text.scale = c(4, 4, 4, 2, 4, 4),
        line.size = 2,
        nsets = 7,
        mainbar.y.label = "Intersection",
        sets.x.label = "Number of genes",
        keep.order = TRUE
)

pdf(file.path(outDir, "05_NeuN_MAG__DEGs__Panel"), width = 18, height = 12)
print(uPlot)
dev.off()

uPlot <- upset(fromList(genesL[c("PV_features", "NeuN_features", "Xenium_panel")]),
        order.by = c("freq", "degree"),
        decreasing = c(T, T),
        number.angles = 0,
        point.size = 20,
        text.scale = c(4, 4, 4, 2, 4, 4),
        line.size = 2,
        mainbar.y.label = "Intersection",
        sets.x.label = "Number of genes",
        keep.order = TRUE
)

pdf(file.path(outDir, "06_Features__Panel.pdf"), width = 18, height = 12)
print(uPlot)
dev.off()

uPlot <- upset(fromList(genesL[c("DOWN_in_TG", "UP_in_TG", "PV_M1", "PV_M2", "PV_M3")]),
        order.by = c("freq", "degree"),
        decreasing = c(T, T),
        number.angles = 0,
        point.size = 20,
        text.scale = c(4, 4, 4, 2, 4, 4),
        line.size = 2,
        nsets = 5,
        mainbar.y.label = "Intersection",
        sets.x.label = "Number of genes",
        keep.order = TRUE
)

pdf(file.path(outDir, "07_DEGs__Wilcoxon.pdf"), width = 18, height = 12)
print(uPlot)
dev.off()

q("no")
