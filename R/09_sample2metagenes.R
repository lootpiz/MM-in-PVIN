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

##### ##### ##### ##### ##### ##### ##### ##### ##### #####
library(RColorBrewer)
library(stringr)
library(igraph)
library(FELLA) # triangle

library(SOTK) # private

add_shape("triangle", clip = shapes("circle")$clip, plot = FELLA:::mytriangle)

for (segment in segments) {
        # Load NMF objects
        nmfObj <- readRDS(file.path(nmfDir, paste0(segment, ".RDS")))

        usage <- NMF::coef(get(as.character(selectRanks[segment]), nmfObj$fit))
        rownames(usage) <- paste0("Metagene_", c(1:selectRanks[segment]))

        # Define nodes/vertex
        nodeSamples <- data.frame(
                name = colnames(usage),
                label = sapply(str_split(colnames(usage), "-"), "[[", 1),
                lblSize = rep(samLblSize, ncol(usage)),
                shape = sapply(colnames(usage), .assignShape),
                size = sapply(colnames(usage), .assignShapeSize),
                color = modelCols[["GeoMx"]][sapply(str_split(colnames(usage), "-"), "[[", 1)]
        )

        nodeMetagenes <- data.frame(
                name = rownames(usage),
                label = paste0("M", sapply(str_split(rownames(usage), "_"), "[[", 2)),
                lblSize = rep(metageneLblSize, nrow(usage)),
                shape = rep("square", nrow(usage)),
                size = rep(metageneNodeSize, nrow(usage)),
                color = "white"
        )

        nodes <- rbind(nodeSamples, nodeMetagenes)

        # Define edges
        cl <- apply(usage, 2, function(x) {
                return(which(x == max(x)))
        })

        edges <- data.frame(
                sample = names(cl),
                metagenes = paste0("Metagene_", cl),
                weight = 1
        )

        graph <- graph.empty()
        graph <- add_vertices(
                        graph, 
                        nv = nrow(nodes),
                        name = as.character(nodes[, 1]),
                        label = as.character(nodes[, 2]),
                        lblSize = as.character(nodes[, 3]),
                        shape = as.character(nodes[, 4]),
                        size = as.character(nodes[, 5]),
                        color = as.character(nodes[, 6])
        )
        graph <- add_edges(graph, t(edges[, c(1:2)]), weight = as.numeric(edges[, 3]))
        graph <- as.undirected(graph) # set.seed() loaded in 00_directorySettings.R
        lay <- layout.fruchterman.reingold(graph)

        pdf(file.path(outDir, segment, "12_Sample2metagene.pdf"), width = figSize, height = figSize)
        plot(graph,
                layout = lay,
                vertex.label = V(graph)$label,
                vertex.label.color = "black",
                vertex.label.cex = as.numeric(V(graph)$lblSize),
                vertex.color = V(graph)$color,
                vertex.shape = V(graph)$shape,
                vertex.size = as.numeric(V(graph)$size),
                edge.weight = E(graph)$weight,
                edge.width = E(graph)$weight
        )
        title(main = segment)
        legend("topleft",
                legend = c("WT", "TG"),
                pch = c(19, 17),
                col = "grey70",
                bty = "n",
                title = "Group",
                pt.cex = 1.5
        )

        legend("topright",
                legend = names(modelCols[["GeoMx"]]),
                col = modelCols[["GeoMx"]],
                pch = 15,
                bty = "n",
                title = "Mice",
                pt.cex = 1.5
        )
        legend("bottomleft",
                legend = c(paste0("NMF rank = ", selectRanks[segment]), paste0(str_to_title(corMet), " >", corrCoefThre)),
                bty = "n"
        )
        legend("bottomright", legend = c("M, Metagene"), bty = "n")
        dev.off()
}

q("no")
