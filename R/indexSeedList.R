#' Calculate gene weight
#' 
#' This function is the core function to calculate weight score for a gene.
#' 
#' @param seedListName A list of vectors which contain gene names in each region.
#' @param dbMatrix Database matrix
#' 
#' @author Jiantao Shi
#' @references
#' 
#' Jiantao Shi, Franziska Michor, Winston Hide: NetGPA, a package for network-based gene prioritization analysis.
#' 
#' @return
#' A list of vectors which contain gene indexes in each region.
#' 
#' @examples
#' # load package and data
#' library("NetGPA")
#' data("Example_NetGPA")
#' 
#' # load data base
#' data(text_2006_12_NetGPA)
#' 
#' # Compare FGF pathway and a random gene set
#' exampleSeedFGFR <- Example_NetGPA$Cancer_GeneSet[["SIGNALING_BY_FGFR"]]
#' 
#' pGenes <- intersect(exampleSeedFGFR, colnames(text_2006_12_NetGPA))
#' rGenes <- sample(colnames(text_2006_12_NetGPA), length(pGenes)) 
#' 
#' res    <- NetGPA(as.list(pGenes), pGenes, text_2006_12_NetGPA, progressBar = FALSE)
#' pTable <- res$queryTable
#' pTable$group <- "Pathway"
#' 
#' res    <- NetGPA(as.list(rGenes), rGenes, text_2006_12_NetGPA, progressBar = FALSE)
#' rTable <- res$queryTable
#' rTable$group <- "Random"
#' 
#' mT <- rbind(pTable, rTable)
#' 
#' boxplot(-log10(queryP)~group, data = mT, ylab = "-log10(p-value)")
#'
#' @export

indexSeedList <- function(seedListName, dbMatrix){

    PG   <- colnames(dbMatrix)
    refIndex <- 1:length(PG)
    names(refIndex) <- PG

    nameL <- sapply(seedListName, function(x) intersect(x, PG))
    nameL <- nameL[sapply(nameL, length) > 0]
    seedList <- sapply(nameL, function(x) refIndex[x])

    cat(length(seedListName), "regions loaded successfully.\n")
    cat(length(seedList),     "regions could be found in database.\n")

    return(seedList)
}
