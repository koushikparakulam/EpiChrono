
# Load all modules
source("Utilities_Preprocessing/utilities.R")
source("Utilities_Preprocessing/preprocessing.R")
source("Functional_Analysis/annotation.R")
source("Functional_Analysis/statistics.R")
source("Enrichment_Analysis/functional_enrichment.R")
source("Enrichment_Analysis/integration_omics.R")
source("Enrichment_Analysis/visualization.R")

main_pipeline <- function(
  bed_path="data/GSM466737_UCSD.H1.H3K36me3.LL160.bed",
  blacklist_path=NULL,
  genome_version="hg19",
  expr_file=NULL,
  meth_file=NULL,
  output_dir="analysis_results"
) {
  init_env()
  check_dependencies(c("rtracklayer","dplyr","GenomicRanges","GenomicFeatures","ggplot2","reshape2"))
  if(!dir.exists(output_dir)) dir.create(output_dir)
  config <- build_config(output_dir)
    
  # Preprocessing
  bed_df <- prepare_bed_data(bed_path, blacklist_path)
  bed_df <- rescale_scores(bed_df)
  bed_df <- finalize_preprocessing(bed_df)
  
  # Annotation
  txdb <- build_txdb(genome=genome_version)
  bed_df <- annotate_regions(bed_df, txdb, add_promoters=TRUE, add_gc=TRUE)
  
  # Stats
  stats <- compute_summary_stats(bed_df)
  model1 <- fit_linear_model(bed_df, "score ~ (end - start)")
  diff_res <- group_comparison(bed_df, group_col="group", metric_col="score")
  
  # Enrichment
  enrich_res <- run_full_enrichment_pipeline(bed_df, gene_col="nearest_gene")
  
  # Integration
  if(!is.null(expr_file) && file.exists(expr_file)) {
    expr_data <- read.csv(expr_file, header=TRUE)
  } else {
    # optionally simulate expression
    genes <- unique(na.omit(bed_df$nearest_gene))
    expr_data <- simulate_expr_data(genes, n_samples=5)
  }
  if(!is.null(meth_file) && file.exists(meth_file)) {
    meth_data <- read.csv(meth_file, header=TRUE)
  } else {
    # simulate
    genes <- unique(na.omit(bed_df$nearest_gene))
    meth_data <- simulate_meth_data(genes)
  }
  group_vec <- rep(c("Control","Treated"), length.out=ncol(expr_data))
  integration_res <- run_integration_pipeline(bed_df, expr_data, meth_data, group_vec)
  
  # Visualization
  plots <- comprehensive_visualization(
    bed_df,
    enrich_res=if(!is.null(enrich_res$basic)) enrich_res$basic else NULL,
    cor_mat=integration_res$cor_matrix,
    kmeans_res=NULL
  )
  if(!is.null(plots$length_hist)) {
    ggsave(file.path(output_dir,"length_hist.png"), plot=plots$length_hist)
  }
  if(!is.null(plots$score_density)) {
    ggsave(file.path(output_dir,"score_density.png"), plot=plots$score_density)
  }
  if(!is.null(plots$enrichment_plot)) {
    ggsave(file.path(output_dir,"enrichment_plot.png"), plot=plots$enrichment_plot)
  }
  if(!is.null(plots$correlation_plot)) {
    ggsave(file.path(output_dir,"correlation_plot.png"), plot=plots$correlation_plot)
  }

  # Save results
  write.csv(bed_df, file.path(output_dir,"final_bed.csv"), row.names=FALSE)
  write.csv(stats, file.path(output_dir,"summary_stats.csv"), row.names=FALSE)
  if(!is.null(diff_res)) {
    capture.output(diff_res, file=file.path(output_dir,"diff_test.txt"))
  }
  if(!is.null(enrich_res$basic$go)) {
    capture.output(enrich_res$basic$go, file=file.path(output_dir,"go_results.txt"))
  }
  if(!is.null(integration_res$cor_matrix)) {
    write.table(integration_res$cor_matrix, file.path(output_dir,"omics_correlation.txt"), sep="\t")
  }

  session_info <- capture.output(sessionInfo())
  writeLines(session_info, con=file.path(output_dir,"session_info.txt"))
}


if(interactive()) {
  main_pipeline()
}
