
#Extracting BED file data
read_bed_file <- function(bed_path) {
  custom_logger(paste("Reading BED:", bed_path))
  gr <- tryCatch(
    {
      rtracklayer::import(bed_path, format="BED")
    },
    error=function(e) {
      custom_logger("Error reading BED")
      stop(e)
    }
  )
  gr
}

# Instantiating Dataframe
convert_gr_to_df <- function(gr) {
  df <- as.data.frame(gr)
  df$seqnames <- as.character(df$seqnames)
  if(!"score" %in% names(df)) df$score <- NA_real_
  df
}


apply_qc_filters <- function(df, min_width=50, max_width=1e6) {
  df$width <- df$end - df$start
  df <- df[df$width >= min_width & df$width <= max_width, ]
  df
}

fix_seqlevels <- function(df, drop_unusual=TRUE) {
  df <- df[grepl("^chr[0-9XY]+$", df$seqnames), ]
  if(drop_unusual) {
    df <- df[!df$seqnames %in% c("chrM", "chrUn"), ]
  }
  df
}

remove_blacklisted_regions <- function(df, blacklist_gr) {
  if(!inherits(blacklist_gr, "GRanges")) return(df)
  df_gr <- GenomicRanges::makeGRangesFromDataFrame(df, keep.extra.columns=TRUE)
  hits <- IRanges::overlapsAny(df_gr, blacklist_gr)
  df_clean <- df[!hits, ]
  df_clean
}

prepare_bed_data <- function(bed_path, blacklist_path=NULL) {
  custom_logger("Starting preprocessing")
  raw_gr <- read_bed_file(bed_path)
  bed_df <- convert_gr_to_df(raw_gr)
  bed_df <- apply_qc_filters(bed_df)
  bed_df <- fix_seqlevels(bed_df)

  if(!is.null(blacklist_path) && file.exists(blacklist_path)) {
    custom_logger("Removing blacklisted regions")
    bl_gr <- tryCatch(
      {
        rtracklayer::import(blacklist_path, format="BED")
      },
      error=function(e) { NULL }
    )
    if(!is.null(bl_gr)) {
      bed_df <- remove_blacklisted_regions(bed_df, bl_gr)
    }
  }
  
  custom_logger("Preprocessing complete")
  bed_df
}

rescale_scores <- function(df, col="score") {
  vals <- df[[col]]
  rng <- range(vals, na.rm=TRUE)
  if(diff(rng)==0) {
    df[[col]] <- 0
  } else {
    df[[col]] <- (vals - rng[1]) / (rng[2] - rng[1])
  }
  df
}

finalize_preprocessing <- function(df) {
  df <- df[order(df$seqnames, df$start), ]
  rownames(df) <- NULL
  df
}
