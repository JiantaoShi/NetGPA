#' The core function to score a gene
#' 
#' This function is the core function to test whether a query gene is statistically associated with a set of seed genes.
#' 
#' @param Nseed Number of seed genes
#' @param Cg Pre-calcuated weight for this gene
#' @param Pfcutoff Percentage of top genes considered, should be less than 0.1.
#' 
#' @importFrom stats dpois
#' @importFrom stats pgamma
#' 
#' @author Jiantao Shi
#' @references
#' 
#' Jiantao Shi: NetGPA, a package for network-based gene prioritization analysis.
#' 
#' @return
#' p-value
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

geneScoring <- function(Nseed, Cg, Pfcutoff = 0.1) {

    Pg <- 0
    lambda <- Nseed*Pfcutoff

    for(i in 0:Nseed)
        Pg <- Pg + dpois(i, lambda)* pgamma(Cg, shape = i, rate = 1, lower.tail = FALSE, log.p = FALSE)

    return(Pg)
}
