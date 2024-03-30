basic_length_hist <- function(df, bin_count=50) {
  ggplot(df, aes(x=end-start)) +
    geom_histogram(bins=bin_count) +
    ggtitle("Region Length Distribution")
}

score_density_plot <- function(df) {
  ggplot(df, aes(x=score)) +
    geom_density() +
    ggtitle("Score Density")
}

plot_gene_enrichment <- function(enrichment_res, top_n=10) {
  if(is.null(enrichment_res)) return(NULL)
  if(nrow(enrichment_res@result) == 0) return(NULL)
  df <- as.data.frame(enrichment_res)
  df <- df[order(df$pvalue), ]
  df_top <- head(df, top_n)
  ggplot(df_top, aes(x=reorder(Description, -pvalue), y=-log10(pvalue))) +
    geom_bar(stat="identity") +
    coord_flip() +
    xlab("Term") + ylab("-log10(P-value)") +
    ggtitle("Top Enriched Terms")
}

visualize_cor_matrix <- function(cor_mat) {
  if(is.null(cor_mat)) return(NULL)
  molten <- reshape2::melt(cor_mat)
  ggplot(molten, aes(x=Var1, y=Var2, fill=value)) +
    geom_tile() + 
    scale_fill_gradient2(mid=0) +
    geom_text(aes(label=round(value,2))) +
    ggtitle("Correlation Matrix")
}

plot_kmeans_clusters <- function(df, km_res, xcol="score", ycol="expression") {
  if(is.null(km_res)) return(NULL)
  df$cluster <- factor(km_res$cluster)
  ggplot(df, aes_string(x=xcol, y=ycol, shape="cluster")) +
    geom_point() +
    ggtitle("K-Means Clustering")
}

comprehensive_visualization <- function(bed_df, enrich_res=NULL, cor_mat=NULL, kmeans_res=NULL) {
  p1 <- basic_length_hist(bed_df)
  p2 <- score_density_plot(bed_df)
  p3 <- if(!is.null(enrich_res)) plot_gene_enrichment(enrich_res$go) else NULL
  p4 <- if(!is.null(cor_mat)) visualize_cor_matrix(cor_mat) else NULL
  p5 <- if(!is.null(kmeans_res)) plot_kmeans_clusters(bed_df, kmeans_res) else NULL
  
  list(
    length_hist=p1,
    score_density=p2,
    enrichment_plot=p3,
    correlation_plot=p4,
    kmeans_plot=p5
  )
}