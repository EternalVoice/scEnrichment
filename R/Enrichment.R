#' @title Enrichment
#' 
#' @param geneset
#' @param genes
#' @param enrichment.matrix
#' @param seurat_obj
#' @param type
#' @param do.sort
#' @param pval.thresh
#' @param pval.n
#' @param do.plot
#' @param do.print
#' 
#' @importFrom stats phyper
#' @export
#' 
Enrichment <- function(geneset, genes=NULL, db=NULL,
                       enrichment.matrix=NULL,
                       seurat_obj = NULL,
                       type = 'GO',
                       species = 'Human',
                       do.sort = FALSE,
                       pval.thresh = 1,
                       pval.n = Inf,
                       do.plot = F,
                       do.print = F) {
  if (is.null(genes) & is.null(enrichment.matrix)) {
    genes <- rownames(seurat_obj)
    em <- BuildEnrichMatrix(genes = genes, type = type,db = db, species = species)
  } else if (!is.null(genes) & is.null(enrichment.matrix)) {
    em <- BuildEnrichMatrix(genes = genes, type = type, db = db, species = species)
  } else if (is.null(genes) & !is.null(enrichment.matrix)) {
    genes <- rownames(enrichment.matrix)
  }
  
  geneset.log <- genes %in% geneset
  present <- geneset.log %*% enrichment.matrix
  pvals <- stats::phyper(present,colSums(enrichment.matrix),colSums(1 - enrichment.matrix),sum(geneset.log),lower.tail = F,log.p = F)
  names(pvals) <- colnames(enrichment.matrix)
  
  if (do.sort) {
    pvals <- sort(pvals)
  }
  pvals <- pvals[pvals <= pval.thresh]
  if (length(pvals) >= pval.n) {
    pvals <- sort(pvals)
    pvals <- pvals[1:pval.n]
  }
  if (do.print) {
    print(pvals)
  }
  if (do.plot) {
    PlotEnrichment(pvals,do.sort = do.sort,pval.thresh = pval.thresh,pval.n = pval.n,do.print = do.print)
  }
  return(pvals)
}
