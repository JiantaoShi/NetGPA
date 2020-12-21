#' NetGPA
#' 
#' The workflow to score a group of query genes
#' 
#' @param seedListName A list which contains the seed gene set.
#' @param queryName A vector of gene names.
#' @param dbMatrix Network database in a matrix format.
#' @param Pfcutoff Percentage of top genes considered, should be less than 0.1.
#' @param progressBar Whether to show process status, logical.
#' 
#' @importFrom utils setTxtProgressBar
#' @importFrom utils txtProgressBar
#' 
#' @author Jiantao Shi
#' @references
#' 
#' Jiantao Shi NetGPA, a package for network-based gene prioritization analysis.
#' 
#' @return
#' A list of two tables
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

NetGPA <- function(seedListName, queryName, dbMatrix, Pfcutoff = 0.1, progressBar = TRUE){

    PG   <- colnames(dbMatrix)
    Nref <- length(PG)

    # check and substract network database
    Pcut <- nrow(dbMatrix)/Nref
    if(Pfcutoff > Pcut)
        stop("The maximum Pfcutoff allowed is", Pcut, "\n")

    Norder <- 1:Nref
    names(Norder) <- PG

    Ncut     <- as.integer(Nref*Pfcutoff)
    dbMatrix <- dbMatrix[1:Ncut, , drop = FALSE]

    # index seed
    seedList <- indexSeedList(seedListName, dbMatrix)
    if(length(seedList) < 10)
        stop("Number of regions are too small, try a larger set of seed genes.\n")

    # check query genes
    queryVector <- intersect(queryName, PG)
    cat(length(queryVector), "genes found in database.\n")

    if(length(queryVector) == 0)
        stop("None of the query of genes are present in database.\n")
    
    # update queryVector to include seed genes
    queryVector <- c(queryVector, unlist(seedListName))
    queryVector <- intersect(queryVector, PG)
    
    Ng <- length(queryVector)

    # assess each gene
    pv <- rep(NA, Ng)

    if(progressBar)
        pb <- txtProgressBar(min = 0, max = Ng, style = 3)

    for(i in 1:Ng){

        g <- queryVector[i]

        queryIndex <- Norder[g]
        topIndex   <- dbMatrix[, g]

        gw    <- geneWeight(queryIndex, Nref, topIndex, seedList, Pfcutoff)
        pv[i] <- geneScoring(gw[1], gw[2], Pfcutoff)

        if(progressBar)
            setTxtProgressBar(pb, i)
    }
    cat("\n")

    # summary results
    names(pv) <- queryVector
    rList <- summaryResults(seedListName, queryName, pv)

    return(rList)
}
