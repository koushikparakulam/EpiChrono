
compute_summary_stats <- function(df, length_col="width", score_col="score") {
  df[[length_col]] <- df$end - df$start
  data.frame(
    n_regions=nrow(df),
    mean_len=mean(df[[length_col]], na.rm=TRUE),
    median_len=median(df[[length_col]], na.rm=TRUE),
    sd_len=sd(df[[length_col]], na.rm=TRUE),
    mean_score=mean(df[[score_col]], na.rm=TRUE),
    sd_score=sd(df[[score_col]], na.rm=TRUE)
  )
}

fit_linear_model <- function(df, formula_str="score ~ width") {
  f <- as.formula(formula_str)
  mod <- lm(f, data=df)
  mod
}

group_comparison <- function(df, group_col="group", metric_col="score") {
  if(!group_col %in% names(df)) {
    df[[group_col]] <- sample(c("A","B"), nrow(df), replace=TRUE)
  }
  t.test(df[[metric_col]] ~ df[[group_col]])
}

fit_mixed_model <- function(df, formula_str="score ~ width + (1|seqnames)") {
  if(!requireNamespace("lme4", quietly=TRUE)) {
    install.packages("lme4"); library(lme4)
  }
  f <- as.formula(formula_str)
  m <- lmer(f, data=df)
  m
}

calc_diff_enrichment <- function(df, col="score", label_col="group") {
  grp_stats <- df %>%
    dplyr::group_by(!!rlang::sym(label_col)) %>%
    dplyr::summarise(avg=mean(!!rlang::sym(col),na.rm=TRUE))
  custom_logger("Diff Enrichment computed")
  grp_stats
}

basic_cor_analysis <- function(df, xcol="score", ycol="gc_content") {
  if(!all(c(xcol, ycol) %in% names(df))) return(NA)
  cor(df[[xcol]], df[[ycol]], use="complete.obs")
}
