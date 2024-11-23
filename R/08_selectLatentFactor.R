##### ##### ##### ##### ##### ##### ##### ##### ##### #####
# Parameter settings
##### ##### ##### ##### ##### ##### ##### ##### ##### #####
setwd(file.path(getwd(), "R"))
source("00_directorySettings.R") # load baseDir <- "/WORKING_DIR"
source("00_settings.R")

dataDir <- file.path(baseDir, "data")
outDir <- file.path(baseDir, "results", "SOTK")

##### ##### ##### ##### ##### ##### ##### ##### ##### #####
library(ComplexHeatmap)

library(SOTK) # private

for (segment in segments) {
        soObj <- readRDS(file.path(dataDir, paste0("sotkObj.", segment, ".RDS")))

        graph <- soObj@corNetwork
        if (segment == "PV") {
                layout <- soObj@unweightedLay
        } else if (segment == "NeuN") {
                layout <- soObj@weightedLay
        }
        colComm <- soObj@commCols
        nodes <- igraph::V(graph)$name
        ranks <- sort(as.numeric(unique(sapply(stringr::str_split(nodes, "\\$"), "[[", 3))))

        pdf(file.path(outDir, segment, "11_coverage.pdf"), width = 10, height = 10)
        for (rank in ranks) {
                vertexLabel <- paste(
                        sapply(stringr::str_split(nodes, "\\$"), "[[", 2),
                        sapply(stringr::str_split(nodes, "\\$"), "[[", 3),
                        sep = "/"
                )
                vertexLabel[which(as.numeric(sapply(stringr::str_split(nodes, "\\$"), "[[", 3)) != rank)] <- ""

                vertexColor <- colComm[igraph::V(graph)$community]
                liveCols <- unique(vertexColor[which(as.numeric(sapply(stringr::str_split(nodes, "\\$"), "[[", 3)) == rank)])
                if (length(liveCols) > 0) {
                        vertexColor[-which(vertexColor %in% liveCols)] <- "white"
                }
                vertexColor[which(vertexLabel != "")] <- "black"
                liveColNode <- length(which(vertexColor == "black"))

                plot(graph,
                        layout = layout,
                        vertex.label = vertexLabel,
                        vertex.color = vertexColor,
                        vertex.label.color = "white",
                        vertex.size = 8,
                        vertex.label.cex = 0.75,
                        edge.color = "grey80"
                )
                title(paste0("Rank = ", rank), cex.main = 2, col.main = "black")
                legend("bottomleft",
                        legend = c(
                                paste0(liveColNode, " metagenes"),
                                paste0(length(colComm), " communities"),
                                paste0("Ratio = ", round(length(colComm) / liveColNode, 2))
                        ),
                        bty = "n"
                )
        }
        dev.off()
}

q("no")
