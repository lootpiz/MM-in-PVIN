##### ##### ##### ##### ##### ##### ##### ##### ##### #####
# Parameter settings
##### ##### ##### ##### ##### ##### ##### ##### ##### #####
setwd(file.path(getwd(), "R"))
source("00_directorySettings.R") # load baseDir <- "/WORKING_DIR"
source("00_settings.R")

dataDir <- file.path(baseDir, "data")
nmfDir <- file.path(baseDir, "data", "NMF")
outDir <- file.path(baseDir, "results", "SOTK")

##### ##### ##### ##### ##### ##### ##### ##### ##### #####
library(stringr)
library(igraph)

library(SOTK) # private

for (segment in segments) {
        soObj <- readRDS(file.path(dataDir, paste0("sotkObj.", segment, ".RDS")))

        graph <- soObj@corNetwork
        graph <- soObj@corNetwork
        if (segment == "PV") {
                layout <- soObj@unweightedLay
        } else if (segment == "NeuN") {
                layout <- soObj@weightedLay
        }
        nodes <- igraph::V(graph)$name
        clusterInfo <- soObj@sample2metagene
        
        for (sample in names(modelCols[["GeoMx"]])) {
                vertexPie <- as.list(rep(1, length(nodes)))
                names(vertexPie) <- nodes
                vertexPieCol <- as.list(rep("white", length(nodes)))
                names(vertexPieCol) <- nodes
                vertexSize <- rep(0, length(nodes))
                names(vertexSize) <- nodes
                
                for (ndx in seq_along(nodes)) {
                        node <- nodes[ndx]
                        samples <- clusterInfo[[segment]][[node]]
                        sdx <- which(sapply(stringr::str_split(samples, "-"), "[[", 1) == sample)

                        if (length(samples) > 0 && length(sdx) > 0) {
                                vertexSize[[node]] <- 5
                                vertexPie[[node]] <- .getModelStats(unique(samples[sdx]))
                                vertexPieCol[[node]] <- modelCols[["GeoMx"]]
                        }
                }
        
                pdf(file.path(outDir, segment, paste0("10_Model_Annotation_", sample, ".pdf")), width=10, height=10)
                plot(graph,
                        layout = layout,
                        vertex.label = NA,
                        vertex.size = as.numeric(vertexSize),
                        vertex.shape = "pie",
                        vertex.pie = vertexPie,
                        vertex.pie.color = vertexPieCol,
                        edge.color = scales::alpha("grey80", 0.3)
                )
                title(main = sample)
                dev.off()
        }
}

q("no")
