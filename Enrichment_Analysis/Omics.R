load_edgeR <- function() {
  if(!requireNamespace("edgeR", quietly=TRUE)) {
    BiocManager::install("edgeR")
  }
  library(edgeR)
}

load_DESeq2 <- function() {
  if(!requireNamespace("DESeq2", quietly=TRUE)) {
    BiocManager::install("DESeq2")
  }
  library(DESeq2)
}

merge_expression_with_bed <- function(bed_df, expr_df, gene_col_bed="nearest_gene", gene_col_expr="gene_id") {
  # merges by gene
  custom_logger("Merging expression with BED info")
  bed_expr <- dplyr::inner_join(bed_df, expr_df, by=setNames(gene_col_expr, gene_col_bed))
  bed_expr
}

perform_dge_analysis_edger <- function(expr_matrix, group_vector) {
  load_edgeR()
  d <- DGEList(counts=expr_matrix, group=group_vector)
  d <- calcNormFactors(d)
  design <- model.matrix(~group_vector)
  d <- estimateDisp(d, design)
  fit <- glmFit(d, design)
  lrt <- glmLRT(fit, coef=2)
  topTags(lrt, n=Inf)
}

perform_dge_analysis_deseq2 <- function(expr_matrix, group_vector) {
  load_DESeq2()
  colData <- data.frame(row.names=colnames(expr_matrix), condition=group_vector)
  dds <- DESeqDataSetFromMatrix(countData=expr_matrix, colData=colData, design=~condition)
  dds <- DESeq(dds)
  res <- results(dds)
  as.data.frame(res)
}

integrate_methylation <- function(bed_df, meth_df, key_col_bed="nearest_gene", key_col_meth="gene", method="avg") {
  # Merges bed info with methylation data
  custom_logger("Integrating methylation with BED data")
  df_merged <- dplyr::inner_join(bed_df, meth_df, by=setNames(key_col_meth, key_col_bed))
  if(method=="avg" && "meth_value" %in% names(df_merged)) {
    df_merged$meth_scaled <- df_merged$meth_value / max(df_merged$meth_value, na.rm=TRUE)
  }
  df_merged
}

multi_omics_integration <- function(bed_df, expr_df, meth_df) {
  # high-level integration
  custom_logger("Running multi-omics integration")
  merged_expr <- merge_expression_with_bed(bed_df, expr_df)
  merged_all <- integrate_methylation(merged_expr, meth_df)
  merged_all
}

correlate_omics_features <- function(df, features=c("score","expression","meth_value")) {
  # compute correlation matrix among selected features
  valid_feats <- features[features %in% colnames(df)]
  if(length(valid_feats) < 2) return(NULL)
  cor(df[valid_feats], use="pairwise.complete.obs")
}

evaluate_multiomics_clusters <- function(df, features=c("score","expression","meth_value"), k=3) {
  subdf <- df[features]
  subdf <- subdf[complete.cases(subdf), ]
  km <- kmeans(subdf, centers=k)
  km
}

run_integration_pipeline <- function(bed_df, expr_df, meth_df, group_vec_expr=NULL) {
  if(!is.null(group_vec_expr) && ncol(expr_df) > 2) {
    # assume expr_df is a matrix of counts
    custom_logger("Performing expression analysis with edgeR")
    edger_res <- perform_dge_analysis_edger(as.matrix(expr_df), group_vec_expr)
    custom_logger("Performing expression analysis with DESeq2")
    deseq_res <- perform_dge_analysis_deseq2(as.matrix(expr_df), group_vec_expr)
  } else {
    edger_res <- NULL
    deseq_res <- NULL
  }
  merged <- multi_omics_integration(bed_df, expr_df, meth_df)
  cormat <- correlate_omics_features(merged)
  custom_logger("Integration pipeline complete")
  list(edgeR=edger_res, DESeq2=deseq_res, integrated=merged, cor_matrix=cormat)
}