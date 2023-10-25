#' @title BuildEnrichMatrix
#' @description Build enrichment matrix
#' @param genes
#' @param type 
#' @param db
#' 
#' @export
#' 
BuildEnrichMatrix <- function(genes, type='GO', db=NULL, species='Human'){
  if (is.null(db)) {
    db <- BuildMSigDB(type, species = species)
  }
  em <- sapply(db, function(x) {genes %in% x})
  rownames(em) <- genes
  em[, colSums(em) == 0] <- NA
  return(em)
}