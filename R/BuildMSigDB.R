
#' @title BuildMSigDB
#' @description Create a reference MSigDB
#' @param type GO, KEGG, HALLMARK, MOTIF, PATHWAYS, BIOCARTA, REACTOME, CELLTYPE.
#' @param species Human or Mouse
#' 
#' @export
#' 
BuildMSigDB <- function(type, species) {
  switch (type,
    GO = {
      db <- MSigDB[[species]]$C5
    },
    HALLMARK = {
      db <- MSigDB[[species]]$H
    },
    MOTIF = {
      db <- MSigDB[[species]]$C3[-grep('UNKNOWN|MIR', names(MSigDB[[species]]$C3))]
    },
    PATHWAYS = {
      db <- MSigDB[[species]]$C2[grep('BIOCARTA|REACTOME|KEGG', names(MSigDB[[species]]$C2))]
    },
    BIOCARTA = {
      db <- MSigDB[[species]]$C2[grep('BIOCARTA', names(MSigDB[[species]]$C2))]
    },
    KEGG = {
      db <- MSigDB[[species]]$C2[grep('KEGG', names(MSigDB[[species]]$C2))]
    },
    REACTOME = {
      db <- MSigDB[[species]]$C2[grep('REACTOME', names(MSigDB[[species]]$C2))]
    },
    CELLTYPE = {
      db <- MSigDB[[species]]$C8
    }
  )
}
