#' Summary of prioritization results
#' 
#' Post-processing of prioritization results
#' 
#' @param seedListName A list which contains the seed gene set.
#' @param queryName A vector of gene names.
#' @param pvalue A vector of p-values.
#' 
#' @importFrom stats p.adjust
#' 
#' @author Jiantao Shi
#' @references
#' 
#' Jiantao Shi: NetGPA, a package for network-based gene prioritization analysis.
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

summaryResults <- function(seedListName, queryName, pvalue){

    # get best prediction for each region
    TG    <- names(pvalue)
    Nr    <- length(seedListName)
    bestP <- rep(NA, Nr)
    bestG <- rep("", Nr)
    N     <- rep(NA, Nr)

    for(i in 1:Nr){

        g <- intersect(seedListName[[i]], TG)
        N[i] <- length(g)

        if(length(g) == 0)
            next

        ming <- sort(pvalue[g])[1]
        bestG[i] <- names(ming)
        bestP[i] <- ming
    }
    Genes     <- sapply(seedListName, function(x) paste0(x, collapse = " "))
    bestFDR   <- p.adjust(bestP)
    seedTable <- data.frame(Genes, N, bestG, bestP, bestFDR)

    # p-value for query
    sharedg    <- intersect(queryName, TG)
    queryP     <- pvalue[sharedg]
    queryFDR   <- p.adjust(queryP)
    queryTable <- data.frame(queryP, queryFDR)
    rownames(queryTable) <- sharedg
    
    # return
    return(list(seedTable = seedTable, queryTable = queryTable))
}
