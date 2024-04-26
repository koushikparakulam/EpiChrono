init_env <- function(seed=1234) {
  set.seed(seed)
  options(stringsAsFactors=FALSE)
}
check_dependencies <- function(pkgs) {
  for (p in pkgs) {
    if (!require(p, character.only=TRUE)) {
      install.packages(p, repos="https://cran.r-project.org"); library(p, character.only=TRUE)
    }
  }
}
custom_logger <- function(msg) {
  cat("[LOG]", Sys.time(), msg, "\n")
}
time_tracker <- function(expr) {
  start <- Sys.time()
  val <- eval(expr)
  custom_logger(paste("Elapsed:", round(difftime(Sys.time(), start, units="secs"),2),"secs"))
  val
}
build_config <- function(base_dir="analysis_results") {
  if(!dir.exists(base_dir)) dir.create(base_dir)
  list(output=base_dir, date=format(Sys.time(),"%Y%m%d"))
}
safe_division <- function(x, y) {
  if(y==0) return(NA_real_) else return(x/y)
}
normalize_vector <- function(v) {
  v / sqrt(sum(v^2, na.rm=TRUE))
}
transform_granges <- function(gr) {
  df <- as.data.frame(gr)
  df$newCol <- df$width / mean(df$width)
  makeGRangesFromDataFrame(df, keep.extra.columns=TRUE)
}
compute_complex_score <- function(df, col="score") {
  s <- df[[col]]
  w <- df$end - df$start
  mean(s * log1p(w), na.rm=TRUE)
}
create_temp_dir <- function(prefix="tmp_util_") {
  tmp <- file.path(tempdir(), paste0(prefix, Sys.getpid()))
  if(!dir.exists(tmp)) dir.create(tmp)
  tmp
}
clean_up_temp <- function(dir_path) {
  if(dir.exists(dir_path)) unlink(dir_path, recursive=TRUE, force=TRUE)
}
apply_random_forest <- function(df, response="score", n_trees=100) {
  if(!require(randomForest)) {
    install.packages("randomForest", repos="https://cran.r-project.org"); library(randomForest)
  }
  f <- as.formula(paste(response,"~ ."))
  rf <- randomForest(f, data=df, ntree=n_trees)
  rf
}
write_html_report <- function(df, file="report.html") {
  html <- paste0("<html><head><title>Report</title></head><body><h1>Data Summary</h1><p>Rows: ",
                nrow(df), "</p><p>Columns: ", ncol(df), "</p></body></html>")
  cat(html, file=file)
}
calc_partial_cor <- function(x, y, z) {
  r_xy <- cor(x,y, use="complete.obs")
  r_xz <- cor(x,z, use="complete.obs")
  r_yz <- cor(y,z, use="complete.obs")
  num <- r_xy - (r_xz*r_yz)
  den <- sqrt((1-r_xz^2)*(1-r_yz^2))
  num/den
}
randomize_columns <- function(df) {
  nm <- sample(names(df))
  df <- df[,nm]
  df
}
smart_subset <- function(df, col="score", threshold=10) {
  df[df[[col]] > threshold, ]
}
perform_batch_processing <- function(lst, fn) {
  res <- list()
  for(i in seq_along(lst)) {
    res[[i]] <- fn(lst[[i]])
  }
  res
}
gather_metrics <- function(df, group_col="group", metric_col="score") {
  df %>% group_by(!!sym(group_col)) %>% summarize(m=mean(!!sym(metric_col), na.rm=TRUE))
}
pseudo_cluster <- function(df, centers=3) {
  if(!require(stats)) {
    stop("Stats package missing!")
  }
  kmeans(df[sapply(df,is.numeric)], centers=centers)
}
summarize_chromosomes <- function(df) {
  df %>% group_by(seqnames) %>% summarize(n=n(), avg_len=mean(end-start))
}
init_env()
check_dependencies(c("dplyr"))
custom_logger("Utilities module loaded")