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
library(SummarizedExperiment)
library(RColorBrewer)
library(tidyverse)
library(gridExtra)
library(ggridges)
library(stringr)
library(ggplot2)
library(NMF)

consensusL <- list()
for (segment in segments) {
        nmfRes <- readRDS(file.path(nmfDir, paste0(segment, ".RDS")))
        matH <- coef(get(as.character(selectRanks[segment]), nmfRes$fit))
        rownames(matH) <- paste0("M", c(1:selectRanks[segment]))
        consensusL[[segment]] <- matH
}
sampleIds <- names(modelCols[["GeoMx"]])

plotL <- list(); consensusRange <- c()
for (sampleId in sampleIds) { 
        message(sampleId)
        consensusGeoMean <- sapply(seq_along(consensusL), function(idx) {
                consensus <- consensusL[[idx]]
                sampleIdx <- which(str_detect(colnames(consensus), sampleId))
                
                if (length(sampleIdx) == 0) {
                        usage <- rep(0, selectRanks[idx])
                } else if (length(sampleIdx) == 1) {
                        usage <- consensus[,sampleIdx]
                } else {
                        usage <- apply(consensus[,sampleIdx], 1, .geoMean)
                }

                return(usage)
        })

        df <- data.frame(
                Sample = sampleId,
                Time = c(1:sum(selectRanks)),
                Segment = rep(segments, selectRanks),
                Metagene = names(unlist(consensusGeoMean)),
                Consensus = unlist(consensusGeoMean)
        )
        df$Segment <- factor(df$Segment, levels=segments)
        consensusRange <- c(consensusRange, df$Consensus)
        
        emptyBar <- 1 # gaps/spaces among segments
        toAdd <- data.frame( matrix(NA, emptyBar * nlevels(df$Segment), ncol(df)) )
        colnames(toAdd) <- colnames(df)
        toAdd$Segment <- rep(levels(df$Segment), each=emptyBar)
        df <- rbind(df, toAdd)
        df <- df %>% arrange(Segment)
        df$id <- seq(1, nrow(df))

        labelData <- df
        numberOfBar <- nrow(labelData)
        angle <- 90 - 360 * (labelData$id-0.5) /numberOfBar
        labelData$hjust <- ifelse(angle < -90, 1, 0)
        labelData$angle <- ifelse(angle < -90, angle+180, angle)
        
        baseData <- df %>% 
                group_by(Segment) %>% 
                summarize(start=min(id), end=max(id) - emptyBar) %>% 
                rowwise() %>% 
                mutate(title=mean(c(start, end)))

        gridData <- baseData
        gridData$end <- gridData$end[ c( nrow(gridData), 1:nrow(gridData)-1)] + 1
        gridData$start <- gridData$start - 1
        gridData <- gridData[-1,]

        cb <- 
        
        ggplot(df, aes(x=as.factor(id), y=Consensus, fill=Segment)) +
                geom_bar(aes(x=as.factor(id), y=Consensus, fill=Segment), stat="identity", alpha=0.5) +
                
                # Add a val=100/75/50/25 lines. I do it at the beginning to make sur barplots are OVER it.
                geom_segment(data=gridData, aes(x = end, y = 0.1, xend = start, yend = 0.1), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
                geom_segment(data=gridData, aes(x = end, y = 0.2, xend = start, yend = 0.2), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
                geom_segment(data=gridData, aes(x = end, y = 0.3, xend = start, yend = 0.3), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
                
                # Add text showing the value of each 100/75/50/25 lines
                annotate("text", x = rep(max(df$id), 3), y = c(0.1, 0.2, 0.3), label = c(".1", ".2", ".3") , color="grey", size=3 , angle=0, fontface="bold", hjust=1) +
                
                geom_bar(aes(x=as.factor(id), y=Consensus, fill=Segment), stat="identity", alpha=0.5) +
                scale_fill_manual(values = segmentCols) + 
                ylim(-0.3,0.5) +
                theme_minimal() +
                theme(
                        legend.position = "none",
                        axis.text = element_blank(),
                        axis.title = element_blank(),
                        panel.grid = element_blank(),
                        plot.margin = unit(rep(-1,4), "cm")
                ) +
                coord_polar() +
                geom_text(data=labelData, aes(x=id, y=Consensus+0.02, label=Metagene, hjust=hjust), color="black", fontface="bold",alpha=0.6, size=5, angle=labelData$angle, inherit.aes = FALSE ) +
                
                # Add base line information
                geom_segment(data=baseData, aes(x = start, y = -0.03, xend = end, yend = -0.03), colour = "black", alpha=0.8, size=1 , inherit.aes = FALSE ) +
                geom_text(data=baseData, aes(x = title, y = -0.1, label=Segment), hjust=c(1,1), colour = "black", alpha=0.8, size=4, fontface="bold", inherit.aes = FALSE) +
                annotate(geom="text", x=0, y=0.5, label=sampleId, color="red", size=3)

        plotL[[sampleId]] <- cb
}

pdf(file.path(outDir, "ConsensusUsage_per_sample.pdf"), width = 4 * 4, height = 4 * 2)
G <- grid.arrange(grobs = plotL, ncol = 4, nrow = 2)
print(G)
dev.off()

q("no")
