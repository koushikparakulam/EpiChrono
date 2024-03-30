# functional_enrichment.R


load_orgdb <- function() {
  # Try loading org.Hs.eg.db for human gene annotation
  if(!requireNamespace("org.Hs.eg.db", quietly=TRUE)) {
    BiocManager::install("org.Hs.eg.db")
  }
  library(org.Hs.eg.db)
}

load_clusterProfiler <- function() {
  # clusterProfiler for GO/KEGG enrichment
  if(!requireNamespace("clusterProfiler", quietly=TRUE)) {
    BiocManager::install("clusterProfiler")
  }
  library(clusterProfiler)
}

enrich_go <- function(gene_vector, universe_vector=NULL, ont="BP", pval_cutoff=0.05) {
  load_orgdb()
  load_clusterProfiler()
  if(is.null(universe_vector)) {
    res <- clusterProfiler::enrichGO(
      gene = gene_vector,
      OrgDb = org.Hs.eg.db,
      ont = ont,
      pvalueCutoff = pval_cutoff,
      keyType = "ENTREZID"
    )
  } else {
    res <- clusterProfiler::enrichGO(
      gene = gene_vector,
      OrgDb = org.Hs.eg.db,
      ont = ont,
      pvalueCutoff = pval_cutoff,
      universe = universe_vector,
      keyType = "ENTREZID"
    )
  }
  res
}

enrich_kegg <- function(gene_vector, pval_cutoff=0.05) {
  load_orgdb()
  load_clusterProfiler()
  res <- clusterProfiler::enrichKEGG(
    gene = gene_vector,
    organism = "hsa",
    pvalueCutoff = pval_cutoff
  )
  res
}

gsea_analysis <- function(ranked_list) {
  load_orgdb()
  load_clusterProfiler()
  res <- GSEA(
    geneList = ranked_list,
    pvalueCutoff = 0.05,
    minGSSize = 10,
    maxGSSize = 500
  )
  res
}

hypergeom_test <- function(gene_set, background_genes, interesting_genes) {
  k <- length(intersect(interesting_genes, gene_set))
  m <- length(gene_set)
  n <- length(background_genes) - m
  x <- length(interesting_genes)
  pval <- phyper(k-1, m, n, x, lower.tail=FALSE)
  pval
}

perform_functional_enrichment <- function(df, gene_col="nearest_gene", background_genes=NULL) {
  
  custom_logger("Starting functional enrichment")
  
  # gather genes
  genes <- unique(na.omit(df[[gene_col]]))
  if(!is.null(background_genes)) {
    # keep only genes in background
    genes <- intersect(genes, background_genes)
  }
  
  if(length(genes) < 3) {
    custom_logger("Not enough genes for enrichment")
    return(NULL)
  }
  # run GO
  go_res <- enrich_go(gene_vector=genes, universe_vector=background_genes)
  # run KEGG
  kegg_res <- enrich_kegg(gene_vector=genes)
  # pack results
  list(go=go_res, kegg=kegg_res)
}

perform_multi_ont_enrichment <- function(genes, ont_list=c("BP","MF","CC")) {
  load_orgdb()
  load_clusterProfiler()
  out_list <- list()
  for(ont in ont_list) {
    enr <- clusterProfiler::enrichGO(
      gene=genes,
      OrgDb=org.Hs.eg.db,
      ont=ont,
      pvalueCutoff=0.05
    )
    out_list[[ont]] <- enr
  }
  out_list
}

simulate_ranked_list <- function(genes) {
  # make a numeric vector for GSEA
  set.seed(123)
  vals <- runif(length(genes), -2, 2)
  names(vals) <- genes
  sort(vals, decreasing=TRUE)
}

run_full_enrichment_pipeline <- function(df, gene_col="nearest_gene") {
  custom_logger("Running full enrichment pipeline")
  
  gvec <- unique(na.omit(df[[gene_col]]))
  if(length(gvec) < 5) {
    custom_logger("Not enough genes, skipping full pipeline")
    return(NULL)
  }
  
  # GO and KEGG
  basic_res <- perform_functional_enrichment(df, gene_col=gene_col)
  
  # GSEA with a simulated ranking
  rnk <- simulate_ranked_list(gvec)
  gsea_res <- gsea_analysis(rnk)
  
  # Multi-ontology
  multi_go <- perform_multi_ont_enrichment(gvec)
  
  list(basic=basic_res, gsea=gsea_res, multi_go=multi_go)
}

batch_enrichment <- function(df_list, gene_col="nearest_gene") {
  # apply enrichment to multiple data splits
  res_list <- list()
  for(i in seq_along(df_list)) {
    custom_logger(paste("Enrichment on subset", i))
    res_list[[i]] <- run_full_enrichment_pipeline(df_list[[i]], gene_col=gene_col)
  }
  res_list
}
