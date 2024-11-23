##### ##### ##### ##### ##### ##### ##### ##### ##### #####
# Parameter settings
##### ##### ##### ##### ##### ##### ##### ##### ##### #####
setwd(file.path(getwd(), "R"))
source("00_directorySettings.R") # load baseDir <- "/WORKING_DIR"

dataDir <- file.path(baseDir, "data")
outDir <- file.path(baseDir, "NMFres")

filename <- "GeoMx.PV" # GeoMx.NeuN

# NMF parameters
rank <- 2 # seq(2, 10) # Run each rank for parallelization
method <- "brunet" # default
nr <- 1000 # iteration
seed <- 123456 # 123456 from the tutorial :p
options <- "tv2p20" # track + verbose more + 20 processors
sharedMemeoryFlag <- TRUE

##### ##### ##### ##### ##### ##### ##### ##### ##### #####
library(NMF)

expr <- readRDS(file.path(dataDir, paste0(filename, ".RDS"))) # NMF input

# Run
resultFile <- file.path(outDir, paste0(filename, "_rank_", rank, ".RDS"))
if (file.exists(resultFile)) {
    NMFres <- readRDS(resultFile)
    message("NMF object loaded.")
} else {
    nmf.options(shared.memory = sharedMemeoryFlag)
    start <- Sys.time()
    NMFres <- nmf(expr, rank = rank, method = method, nrun = nr, seed = seed, .options = options)
    end <- Sys.time()
    message(end - start)
    saveRDS(NMFres, resultFile)
    message("NMF object saved.")
}

# Permute data and recalculate to estimate if overfitting is happening
resultRandFile <- file.path(outDir, paste0(filename, "_rank_", rank, "_rand.RDS"))
if (file.exists(resultRandFile)) {
    NMFresRand <- readRDS(resultRandFile)
    message("NMF_rand object loaded.")
} else {
    exprRand <- randomize(expr)
    start <- Sys.time()
    NMFresRand <- nmf(exprRand, rank = rank, method = method, nrun = nr, seed = seed, .options = options)
    end <- Sys.time()
    message(end - start)
    saveRDS(NMFresRand, resultRandFile)
    message("NMF_rand object saved.")
}

q("no")
