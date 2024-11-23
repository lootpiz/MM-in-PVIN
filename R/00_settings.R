##### ##### ##### ##### ##### ##### ##### ##### ##### #####
# GeoMx parameters
corMet = "spearman"
corrCoefThre <- 0.5 # correlation coefficient threshold values to include
seed <- 1118
niter <- 1000 # graph layout iteration
commWeight <- 100 # graph layout weights on community info
cohortWeight <- 1 # graph layout weights on cohort info

samLblSize <- "1" # sample (ROIs) label size
metageneNodeSize <- "10"
metageneLblSize <- "1.5" # metagene label size

figSize <- 10 # width and height
vertexSize <- 5
vertexLabelCex <- 1
edgeAlpha <- 0.9

##### ##### ##### ##### ##### ##### ##### ##### ##### #####
# Xenium parameters
wilcoxonDegFcThre <- 0.5
wilcoxonDegPvalThre <- 0.05

##### ##### ##### ##### ##### ##### ##### ##### ##### #####
# Color key
groupCols <- c(
        "TG"  = "#7FFF00", # transgenic mice 
        "WT" = "#FF931E"   # wild-type
)

segmentCols <- c(
        "PV"      = "#7570B3",
        "NeuN"    = "#1B9E77"
        # "Amyloid" = "#E7298A",
        # "TN"      = "#D4AF37"
)
segments <- c("PV", "NeuN")

modelCols <- list(
        "GeoMx" = c(
                "TG1F" = "#084594",
                "TG2F" = "#2171B5",
                "TG3F" = "#4292C6",
                "TG4F" = "#6BAED6",
                "WT1F" = "#99000D",
                "WT2F" = "#CB181D",
                "WT3F" = "#EF3B2C"
        ),
        "Xenium" = c(
                "TG2F" = "#2171B5",
                "TG3F" = "#4292C6",
                "TG4F" = "#6BAED6",
                "WT1F" = "#99000D",
                "WT2F" = "#CB181D",
                "WT3F" = "#EF3B2C"
        )
)

sectionCols <- c(
        "RSC" = "#1B9E77",
        "SUB" = "#D95F02",
        "VIS" = "#7570B3",
        "ENT" = "#E7298A",
        "CA1" = "#FFBF00"
)

kmeansCols <- c("#7FC97F", "#BEAED4", "#FDC086", "#FFFF99", "#386CB0", "#F0027F", "#BF5B17")
names(kmeansCols) <- c(1:7)

metagenesCols <- c("#7570b3", "#9c98c8", "#c3c1de", "#1b9e77", "#5cb99d", "#7cc7b1", "#bde3d8")
names(metagenesCols) <- c("PV_M1", "PV_M2", "PV_M3", "NeuN_M1", "NeuN_M2", "NeuN_M3", "NeuN_M4")


##### ##### ##### ##### ##### ##### ##### ##### ##### #####
# User defined functions
.geoMean <- function(x) exp(mean(log(x)))

.listIntersection <- function(listInput, sort = TRUE) {
        listInputmat <- fromList(listInput) == 1
        listInputunique <- unique(listInputmat)
        grouplist <- list()
        for (i in 1:nrow(listInputunique)) {
                currentRow <- listInputunique[i, ]
                myelements <- which(apply(listInputmat, 1, function(x) all(x == currentRow)))
                attr(myelements, "groups") <- currentRow
                grouplist[[paste(colnames(listInputunique)[currentRow], collapse = ":")]] <- myelements
                myelements
        }
        if (sort) {
                grouplist <- grouplist[order(sapply(grouplist, function(x) length(x)), decreasing = TRUE)]
        }
        attr(grouplist, "elements") <- unique(unlist(listInput))
        return(grouplist)
}

.ttest <- function(dat, idx1, idx2){
        obj <- try(t.test(dat[idx1], dat[idx2], var.equal=FALSE), silent=TRUE)
        if(is(obj, "try-error")) {
                value <- NA
        }else{
                value <- obj$p.value
        }
        return(value)
}

.log2fc <- function(dat, idx1, idx2){
        fc <- log2(mean(dat[idx1]) / mean(dat[idx2]))
        return(fc)
}

.reformP <- function(p) {
    if(!is.na(p)) {
        if (p < 0.0001) {
            pString <- formatC(p, format = "e", digits = 2)
        } else {
            pString <- format(round(p, 4), nsmall = 4)
        }
    } else {
        pString <- "NA"
    }
    return(pString)
}

.getModelStats <- function(x) {
        TG1F <- 0; TG2F <- 0; TG3F <- 0; TG4F <- 0; 
        WT1F <- 0; WT2F <- 0; WT3F <- 0

        if (!is.null(x)) {
                for (i in x) {
                        buff <- unlist(stringr::str_split(i, "-"))
                        if (buff[1] == "TG1F") {
                                TG1F <- TG1F + 1
                        } else if (buff[1] == "TG2F") {
                                TG2F <- TG2F + 1
                        } else if (buff[1] == "TG3F") {
                                TG3F <- TG3F + 1
                        } else if (buff[1] == "TG4F") {
                                TG4F <- TG4F + 1
                        } else if (buff[1] == "WT1F") {
                                WT1F <- WT1F + 1
                        } else if (buff[1] == "WT2F") {
                                WT2F <- WT2F + 1
                        } else if (buff[1] == "WT3F") {
                                WT3F <- WT3F + 1
                        } else {
                                message(paste0("WARNING::", buff[1], " is uncategorized."))
                        }
                }
        }
        return(c(TG1F, TG2F, TG3F, TG4F,
                WT1F, WT2F, WT3F))
}

.assignShape <- function(x) {
        group <- substring(x, 1, 1)
        if (group == "T") {
                shape <- "triangle" # transgenic
        } else if (group == "W") {
                shape <- "circle" # wild-type
        }
        return(shape)
}

.assignShapeSize <- function(x) {
        group <- substring(x, 1, 1)
        if (group == "T") {
                size <- "9" # triangle
        } else if (group == "W") {
                size <- "7" # circle
        }
        return(size)
}
