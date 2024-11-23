##### ##### ##### ##### ##### ##### ##### ##### ##### #####
# Parameter settings
##### ##### ##### ##### ##### ##### ##### ##### ##### #####
setwd(file.path(getwd(), "R"))
source("00_directorySettings.R") # load baseDir <- "/WORKING_DIR"

dataDir <- file.path(baseDir, "data")
outDir <- file.path(baseDir, "results")

##### ##### ##### ##### ##### ##### ##### ##### ##### #####
library(stringr)
library(UpSetR)

featureL <- readRDS(file.path(dataDir, "FeaturesL.RDS"))

uPlot <- upset(fromList(featureL),
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

dir.create(file.path(outDir, "Stats"), showWarnings = FALSE)

pdf(file.path(outDir, "Stats", "01_Feature_comparison.pdf"), width = 18, height = 12)
print(uPlot)
dev.off()

q("no")
