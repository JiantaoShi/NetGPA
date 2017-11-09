#' Calculate gene weight
#' 
#' This function is the core function to calculate weight score for a gene.
#' 
#' @param queryIndex The index of query gene.
#' @param Nref Total number of genes in the database. 
#' @param topIndex Top indexes that are associated with queryIndex gene
#' @param seedList A list of intger vectors, each of which represents indexes of seed genes.
#' @param Pfcutoff Percentage of top genes considered, should be less than 0.1.
#' 
#' @author Jiantao Shi
#' 
#' @references
#' Jiantao Shi: NetGPA, a package for network-based gene prioritization analysis.
#' 
#' @return
#' A vector with 2 variables, Nseeds and Cg.
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

geneWeight <- function(queryIndex, Nref, topIndex, seedList, Pfcutoff){

    Ncut     <- length(topIndex)
    topRank  <- 1:Ncut

    # effective seed regions
    check <- sapply(seedList, function(x) (queryIndex %in% x))
    seedL <- seedList[!check]
    Nseed <- length(seedL)

    # work on shared Rank
    sharedRank <- sapply(seedL, function(x) topRank[topIndex %in% x])

    Index <- sapply(sharedRank, length) > 0

    if(sum(Index) == 0)
        return(c(Nseed = Nseed, Cg = 0))

    sharedRank <- sharedRank[Index]

    # weight
    pvalue <- sapply(sharedRank, function(x) (1 - (1 - min(x)/Nref)^length(x)) )
    pvalue[pvalue > Pfcutoff] <- Pfcutoff
    Cg <- -sum(log(pvalue/Pfcutoff))

    return(c(Nseed = Nseed, Cg = Cg))
}
