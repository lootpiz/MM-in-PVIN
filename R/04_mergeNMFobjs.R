##### ##### ##### ##### ##### ##### ##### ##### ##### #####
# Parameter settings
##### ##### ##### ##### ##### ##### ##### ##### ##### #####
setwd(file.path(getwd(), "R"))
source("00_directorySettings.R") # load baseDir <- "/WORKING_DIR"
source("00_settings.R")

dataDir <- file.path(baseDir, "data", "NMF")

##### ##### ##### ##### ##### ##### ##### ##### ##### #####
library(NMF)

library(SOTK) # private

filename <- "PV" # PV, NeuN - run it for each cohort/region

res2 <-  readRDS(file.path(dataDir, filename, paste0(filename, "_rank_2.RDS")));  attr(res2, "fit")
res3 <-  readRDS(file.path(dataDir, filename, paste0(filename, "_rank_3.RDS")));  attr(res3, "fit")
res4 <-  readRDS(file.path(dataDir, filename, paste0(filename, "_rank_4.RDS")));  attr(res4, "fit")
res5 <-  readRDS(file.path(dataDir, filename, paste0(filename, "_rank_5.RDS")));  attr(res5, "fit")
res6 <-  readRDS(file.path(dataDir, filename, paste0(filename, "_rank_6.RDS")));  attr(res6, "fit")
res7 <-  readRDS(file.path(dataDir, filename, paste0(filename, "_rank_7.RDS")));  attr(res7, "fit")
res8 <-  readRDS(file.path(dataDir, filename, paste0(filename, "_rank_8.RDS")));  attr(res8, "fit")
res9 <-  readRDS(file.path(dataDir, filename, paste0(filename, "_rank_9.RDS")));  attr(res9, "fit")
res10 <- readRDS(file.path(dataDir, filename, paste0(filename, "_rank_10.RDS"))); attr(res10, "fit")

objL <- list(`2` = res2, `3` = res3, `4` = res4, `5` = res5, 
        `6` = res6, `7` = res7, `8` = res8, `9` = res9, `10` = res10
)
nmfObj <- mergeNMFObjs(objL)
saveRDS(nmfObj, file.path(dataDir, paste0(filename, ".RDS")))

res2_rand <-  readRDS(file.path(dataDir, filename, paste0(filename, "_rank_2_rand.RDS")));  attr(res2_rand, "fit")
res3_rand <-  readRDS(file.path(dataDir, filename, paste0(filename, "_rank_3_rand.RDS")));  attr(res3_rand, "fit")
res4_rand <-  readRDS(file.path(dataDir, filename, paste0(filename, "_rank_4_rand.RDS")));  attr(res4_rand, "fit")
res5_rand <-  readRDS(file.path(dataDir, filename, paste0(filename, "_rank_5_rand.RDS")));  attr(res5_rand, "fit")
res6_rand <-  readRDS(file.path(dataDir, filename, paste0(filename, "_rank_6_rand.RDS")));  attr(res6_rand, "fit")
res7_rand <-  readRDS(file.path(dataDir, filename, paste0(filename, "_rank_7_rand.RDS")));  attr(res7_rand, "fit")
res8_rand <-  readRDS(file.path(dataDir, filename, paste0(filename, "_rank_8_rand.RDS")));  attr(res8_rand, "fit")
res9_rand <-  readRDS(file.path(dataDir, filename, paste0(filename, "_rank_9_rand.RDS")));  attr(res9_rand, "fit")
res10_rand <- readRDS(file.path(dataDir, filename, paste0(filename, "_rank_10_rand.RDS"))); attr(res10_rand, "fit")

objL_rand <- list(`2` = res2_rand, `3` = res3_rand, `4` = res4_rand, `5` = res5_rand, 
        `6` = res6_rand, `7` = res7_rand, `8` = res8_rand, `9` = res9_rand, `10` = res10_rand
)
nmfObj_rand <- mergeNMFObjs(objL_rand)
saveRDS(nmfObj_rand, file.path(dataDir, paste0(filename, "_rand.RDS")))

q("no")
