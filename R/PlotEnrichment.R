#' @title PlotEnrichment
#' 
#' @param pvals
#' @param do.sort
#' @param pval.thresh
#' @param pval.n
#' @param do.plot
#' @param do.print
#' @param ...
#' 
#' @importFrom graphics barplot text
#' @export
#' 
PlotEnrichment <- function(pvals, 
                           do.sort = F,
                           pval.thresh = 1,
                           pval.n = Inf,
                           do.print = F) {
  pvals <- pvals[pvals < pval.thresh]
  if (do.sort) {
    pvals <- sort(pvals)
  }
  pvals <- pvals[pvals < pval.thresh]
  if (length(pvals) >= pval.n) {
    pvals <- pvals[1:pval.n]
  }
  if (do.print) {
    print(pvals)
  }
  z <- pvals == 0
  pvals[z] <- 10^-10
  names(pvals)[z] <- paste0(names(pvals)[z], '*')
  if (length(pvals) > 0) {
    barplot(-log10(pvals), horiz = T, xlab = '-log10(p-value)', space=0, names.arg='')
    text(0, (1:length(pvals))-0.5, labels=names(pvals), pos=4, cex=0.7)
  }
}
