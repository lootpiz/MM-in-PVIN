##### ##### ##### ##### ##### ##### ##### ##### ##### #####
# Parameter settings
##### ##### ##### ##### ##### ##### ##### ##### ##### #####
setwd(file.path(getwd(), "R"))
source("00_directorySettings.R") # load baseDir <- "/WORKING_DIR"
source("00_settings.R")

dataDir <- file.path(baseDir, "data")
nmfDir <- file.path(baseDir, "data", "NMF")
outDir <- file.path(baseDir, "results")

##### ##### ##### ##### ##### ##### ##### ##### ##### #####
library(SOTK) # private

dir.create(file.path(outDir, "SOTK"), showWarnings = FALSE)

## Load NMF/cNMF results/objects
for (segment in segments) {
        dir.create(file.path(outDir, "SOTK", segment), showWarnings = FALSE)
         
        nmfObj <- readRDS(file.path(nmfDir, paste0(segment, ".RDS")))

        dataL <- list(nmfObj); names(dataL) <- segment
        if (segment == "PV") {
                rankL <- list(c(2:6)); names(rankL) <- segment # half of the population
        } else if (segment == "NeuN") {
                rankL <- list(c(2:7)); names(rankL) <- segment
        }

        ## Concat W matrics and calculate pairwise correlations
        soSet <- SOSet(NMFobjL = dataL, NMFrankL = rankL, dataCol = segmentCols, corMet = corMet)
        plotCorrDensity(soSet, filename = file.path(outDir, "SOTK", segment, "01_Stats_correlation_density.pdf"), width = 6, height = 6)

        ## Generate a correlation network and apply community search
        soObj <- SOTK(soSet, coefThre = corrCoefThre, seed = seed, niter = niter, commWeight = commWeight, cohortWeight = cohortWeight)
        saveRDS(soObj, file.path(dataDir, paste0("sotkObj.", segment, ".RDS")))

        plotNetwork(soObj, annot = "cohort", edgeAlpha = 0.8, weighted = FALSE, filename = file.path(outDir, "SOTK", segment, "02_Network_Unweighted.pdf"), vertexSize = 10, vertexLabelCex = 2) # Correlation network at the Metagene-level
        plotNetwork(soObj, label = TRUE, annot = "cohort", edgeAlpha = 0.8, weighted = FALSE, filename = file.path(outDir, "SOTK", segment, "03_Network_Unweighted_lbl.pdf"), vertexSize = 10, vertexLabelCex = 2) 
        plotNetwork(soObj, annot = "community", edgeAlpha = 0.8, weighted = FALSE, filename = file.path(outDir, "SOTK", segment, "04_Network_Unweighted_Comm.pdf"), vertexSize = 10, vertexLabelCex = 2) # Correlation network at the Metagene-level
        plotNetwork(soObj, label = TRUE, annot = "community", edgeAlpha = 0.8, weighted = FALSE, filename = file.path(outDir, "SOTK", segment, "05_Network_Unweighted_Comm_lbl.pdf"), vertexSize = 10, vertexLabelCex = 2)
        plotNetwork(soObj, annot = "community", edgeAlpha = 0.8, weighted = TRUE, filename = file.path(outDir, "SOTK", segment, "06_Network_Community.pdf"), vertexSize = 10, vertexLabelCex = 2) # Community info
        plotNetwork(soObj, label = TRUE, annot = "community", edgeAlpha = 0.8, weighted = TRUE, filename = file.path(outDir, "SOTK", segment, "07_Network_Community_lbl.pdf"), vertexSize = 10, vertexLabelCex = 2)
        plotNetwork(soObj, annot = "cohort", edgeAlpha = 0.8, weighted = TRUE, filename = file.path(outDir, "SOTK", segment, "08_Network_Weighted.pdf"), vertexSize = 10, vertexLabelCex = 2) # Metagenes-level correlation network with weighted values
        plotNetwork(soObj, label = TRUE, annot = "cohort", edgeAlpha = 0.8, weighted = TRUE, filename = file.path(outDir, "SOTK", segment, "09_Network_Weighted_lbl.pdf"), vertexSize = 10, vertexLabelCex = 2)
}

q("no")
