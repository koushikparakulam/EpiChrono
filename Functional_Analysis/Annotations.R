
build_txdb <- function(genome="hg19") {
  custom_logger(paste("Building TxDb for", genome))
  txdb <- tryCatch(
    {
      GenomicFeatures::makeTxDbFromUCSC(genome=genome, tablename="refGene")
    },
    error=function(e) {
      custom_logger("TxDb creation failed")
      stop(e)
    }
  )
  txdb
}

annotate_nearest_gene <- function(df, txdb) {
  gr <- GenomicRanges::makeGRangesFromDataFrame(df, keep.extra.columns=TRUE)
  gene_gr <- GenomicFeatures::genes(txdb)
  hits <- GenomicRanges::nearest(gr, gene_gr)
  df$nearest_gene <- NA
  df$nearest_gene_dist <- NA
  idx <- !is.na(hits)
  df$nearest_gene[idx] <- gene_gr$gene_id[hits[idx]]
  df$nearest_gene_dist[idx] <- GenomicRanges::distance(gr[idx], gene_gr[hits[idx]])
  df
}

annotate_promoters <- function(df, txdb, upstream=2000, downstream=200) {
  tss <- GenomicFeatures::transcripts(txdb)
  prom <- GenomicRanges::promoters(tss, upstream=upstream, downstream=downstream)
  bed_gr <- GenomicRanges::makeGRangesFromDataFrame(df, keep.extra.columns=TRUE)
  overlap <- IRanges::findOverlaps(bed_gr, prom)
  df$in_promoter <- FALSE
  df$in_promoter[S4Vectors::queryHits(overlap)] <- TRUE
  df
}

annotate_gc_content <- function(df, ref_genome="BSgenome.Hsapiens.UCSC.hg19") {
  if(!requireNamespace(ref_genome, quietly=TRUE)) return(df)
  custom_logger("Annotating GC content")
  library(ref_genome, character.only=TRUE)
  bed_gr <- GenomicRanges::makeGRangesFromDataFrame(df, keep.extra.columns=TRUE)
  seqs <- BSgenome::getSeq(get(ref_genome), bed_gr)
  gc_vals <- letterFrequency(seqs, letters=c("G","C"))
  gc_pct <- rowSums(gc_vals) / width(seqs)
  df$gc_content <- gc_pct
  df
}

annotate_regions <- function(df, txdb, add_promoters=TRUE, add_gc=FALSE) {
  df2 <- annotate_nearest_gene(df, txdb)
  if(add_promoters) {
    df2 <- annotate_promoters(df2, txdb)
  }
  if(add_gc) {
    df2 <- annotate_gc_content(df2)
  }
  df2
}