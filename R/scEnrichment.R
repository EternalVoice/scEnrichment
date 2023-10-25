#' @title scEnrichment
#' @description Perform enrichment analysis on single-cell level
#' 
#' @param markers Output of `FindMarkers` or `FindAllmarkers` in Seurat package.
#' @param seurat_obj
#' @param enrich.method
#' @param type 
#' @param output
#' 
#' @import clusterProfiler
#' @import org.Hs.eg.db
#' @import org.Mm.eg.db
#' @importFrom dplyr pull
#' 
#' @export
#' 
scEnrichment <- function(markers, genes=NULL, db=NULL,
                         seurat_obj = NULL,
                         enrich.method='clusterProfiler',
                         species = 'Human',
                         type = 'GO',
                         do.sort = FALSE,
                         pval.thresh = 1,
                         pval.n = Inf,
                         do.plot = F,
                         do.print = F) {
  
  markers <- markers[order(markers$cluster, decreasing = FALSE),]
  df <- split(markers, markers$cluster)
  geneset <- lapply(df, function(x) x$gene)
  
  if (enrich.method == 'clusterProfiler') {
    if (species == 'Human') {
      orgdb <- 'org.Hs.eg.db'
      organism <- 'hsa'
    } else if (species == 'Mouse') {
      orgdb <- 'org.Mm.eg.db'
      organism <- 'mmu'
    }
    
    if (type == 'GO') {
      ego_BP_All <- lapply(geneset, function(gene){
        ego_BP <- enrichGO(
          gene = gene, OrgDb = orgdb, keyType = 'SYMBOL', ont = 'BP',
          pAdjustMethod = 'BH', pvalueCutoff = 0.01, qvalueCutoff = 0.05
        )
        ego_BP@result$Description <- substring(ego_BP@result$Description, 1, 70)
        ego_BP@result
      })
      names(ego_BP_All) <- names(geneset)
      return(ego_BP_All)
    } else if (type == 'KEGG') {
      ekegg_All <- lapply(geneset, function(gene){
        genelist <- bitr(gene, fromType = 'SYMBOL', toType = 'ENTREZID', OrgDb = orgdb)
        genelist <- dplyr::pull(genelist, ENTREZID)
        ekegg <- enrichKEGG(gene = genelist, organism = organism)
        ekegg@result
      })
      names(ekegg_All) <- names(geneset)
      return(ekegg_All)
    }
  } else if (enrich.method == 'custom') {
    em <- BuildEnrichMatrix(rownames(seurat_obj), type = type, species = species)
    enrich.out <- sapply(geneset, function(gene){
      e <- Enrichment(
        geneset = gene, genes = genes, db = db, enrichment.matrix = em, seurat_obj = seurat_obj,
        type = type, species = species, do.sort = do.sort, pval.thresh = pval.thresh,
        pval.n = pval.n, do.plot = do.plot, do.print = do.print
      )
      names(sort(e))
    })
    return(enrich.out)
  }
}
