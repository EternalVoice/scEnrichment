## scEnrichment: single-cell functional enrichment analysis

#### Install

```r
devtools::install_github('EternalVoice/scEnrichment')
```

#### Usage

```r
library(Seurat)
data('MSigDB')

markers <- FindAllmarkers(pbmc_small, only.pos = TRUE)
enrich_out <- scEnrichment(
    markers = markers,
    seurat_obj = pbmc_small,
    enrich.method = 'clusterProfiler',
    species = 'Human',
    type = 'GO'
)
enrich_out_2 <- scEnrichment(
    markers = markers,
    seurat_obj = pbmc_small,
    enrich.method = 'custom',
    species = 'Human',
    type = 'REACTOME'
)
```
