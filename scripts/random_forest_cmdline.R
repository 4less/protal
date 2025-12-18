# Random forest trainer modeled after random_forest.Rmd (caret + retrain)
# Usage:
# Rscript scripts/random_forest_cmdline.R \
#   --truth-file /path/to/all.truth_annotated \
#   --output-prefix /path/to/output/random_forest_caret \
#   --ntree 512 \
#   --maxnodes 0 \
#   --test-fraction 0.2
# Optional: --seed 1234

suppressPackageStartupMessages({
  library(randomForest)
  library(pmml)
  library(caret)
  library(dplyr)
  library(doParallel)
})



sensitivity <- function(ct) {
  tp <- ct[2, 2]; fp <- ct[2, 1]; fn <- ct[1, 2]
  tp / (tp + fn)
}

precision <- function(ct) {
  tp <- ct[2, 2]; fp <- ct[2, 1]
  tp / (tp + fp)
}

f1 <- function(ct) {
  tp <- ct[2, 2]; fp <- ct[2, 1]; fn <- ct[1, 2]
  (2 * tp) / (2 * tp + fp + fn)
}

evaluate_split <- function(rf, test_data, train_data) {
  pred_test <- predict(rf, newdata = test_data)
  pred_train <- predict(rf, newdata = train_data)
  ct_test <- table(pred_test, test_data$truth)
  ct_train <- table(pred_train, train_data$truth)
  list(
    test = list(ct = ct_test, sensitivity = sensitivity(ct_test), precision = precision(ct_test), f1 = f1(ct_test)),
    train = list(ct = ct_train, sensitivity = sensitivity(ct_train), precision = precision(ct_train), f1 = f1(ct_train))
  )
}

fill_in_allele_info <- function(df) {
  df$total_a <- df$A0 + df$A1 + df$A2 + df$A3 + df$A4
  df$total_af <- df$AF0 + df$AF1 + df$AF2 + df$AF3 + df$AF4
  df$RA0 <- df$A0 / df$total_a
  df$RA1 <- df$A1 / df$total_a
  df$RA2 <- df$A2 / df$total_a
  df$RA3 <- df$A3 / df$total_a
  df$RA4 <- df$A4 / df$total_a
  df$RAF0 <- df$AF0 / df$total_af
  df$RAF1 <- df$AF1 / df$total_af
  df$RAF2 <- df$AF2 / df$total_af
  df$RAF3 <- df$AF3 / df$total_af
  df$RAF4 <- df$AF4 / df$total_af
  df$RA0[is.nan(df$RA0)] <- 0
  df$RA1[is.nan(df$RA1)] <- 0
  df$RA2[is.nan(df$RA2)] <- 0
  df$RA3[is.nan(df$RA3)] <- 0
  df$RA4[is.nan(df$RA4)] <- 0
  df$RAF0[is.nan(df$RAF0)] <- 0
  df$RAF1[is.nan(df$RAF1)] <- 0
  df$RAF2[is.nan(df$RAF2)] <- 0
  df$RAF3[is.nan(df$RAF3)] <- 0
  df$RAF4[is.nan(df$RAF4)] <- 0
  df
}

parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  opts <- list()
  i <- 1
  while (i <= length(args)) {
    key <- args[[i]]
    if (startsWith(key, "--")) {
      if (i == length(args)) stop(paste("Missing value for", key))
      val <- args[[i + 1]]
      if (key == "--truth-file") opts$truth_file <- val
      else if (key == "--ntree") opts$ntree <- as.integer(val)
      else if (key == "--maxnodes") opts$maxnodes <- as.integer(val)
      else if (key == "--output-prefix") opts$output_prefix <- val
      else if (key == "--seed") opts$seed <- as.integer(val)
      else if (key == "--test-fraction") opts$test_fraction <- as.numeric(val)
      else if (key == "--threads") opts$threads <- as.numeric(val)
      else stop(paste("Unknown option", key))
      i <- i + 2
    } else {
      stop(paste("Unexpected argument", key))
    }
  }
  if (is.null(opts$truth_file) || is.null(opts$output_prefix)) {
    stop("Required: --truth-file <path> --output-prefix <path>")
  }
  if (is.null(opts$ntree)) opts$ntree <- 256L
  if (is.null(opts$maxnodes)) opts$maxnodes <- 128L
  if (is.null(opts$test_fraction)) opts$test_fraction <- 0.2
  if (is.null(opts$threads)) opts$threads <- 4
  opts
}

save_forest <- function(model, file) {
  pmod <- pmml(model, file)
  save_pmml(pmod, file)
}

load_truth_data <- function(path) {
  df <- tryCatch(
    read.delim(path, header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE),
    error = function(e) stop("Failed to read truth file: ", e$message)
  )
  # df <- df[, seq_len(length(base_cols))]
  # names(df) <- base_cols
  df$dataset <- "dataset"
  df$truth_raw <- as.logical(df$truth)
  df$prediction <- as.logical(df$prediction)
  df$truth <- as.factor(df$truth_raw)
  inf_cols <- c("lu_gene_rate", "lsu_gene_rate", "lu_gene_rate2", "lsu_gene_rate2", "lu_gene_rate3", "lsu_gene_rate3",
                "lu_rate", "lsu_rate", "su_rate", "lsu_per_read", "lu_per_read")
  for (col in inf_cols) {
    bad <- is.infinite(df[[col]])
    if (any(bad)) df[[col]][bad] <- 0
  }
  #df <- fill_in_allele_info(df)
  df
}

main <- function() {
  opts <- parse_args()
  if (!is.null(opts$seed)) set.seed(opts$seed)

  print("Load data")
  data <- load_truth_data(opts$truth_file)

  # Feature set from random_forest.Rmd (columns_with_unique2 without truth)
  feature_cols <- c(
    "present_genes", "unique_hits", "mean_ani", "expected_gene_presence",
    "expected_gene_presence_ratio", "uniqueness", "mean_mapq", "variance1", "variance2", "A0", "A1", "A2",
    "A3", "A4", "AF0", "AF1", "AF2", "AF3", "AF4", "stddev", "hittable", "lu", "lu_genes", "lsu",
    "lsu_genes", "su_genome", "lu_genome", "lsu_genome", "total_genome", "su_rate", "lu_rate", "lsu_rate",
    "lu_gene_rate", "lsu_gene_rate", "lu_gene_rate2", "lsu_gene_rate2", "lu_gene_rate3", "lsu_gene_rate3",
    "lsu_per_read", "lu_per_read"
  )

  feature_cols <- colnames(data) %>% 
    setdiff(c("total_hits", "truth", "truth_raw", "taxon", "prediction", "dataset")) %>%
    setdiff(c(""))

  # caret CV to choose mtry
  trControl <- trainControl(
    method = "cv",
    number = 5,
    search = "grid",
    savePredictions = TRUE,
    allowParallel = TRUE
  )
  tuneGrid <- expand.grid(mtry = seq(6, 25, by = 1))

  caret_tmp <- "tmp_caret_fit.rds"
  if (file.exists(caret_tmp)) {
    cat("Load caret_tmp")
    caret_fit <- readRDS(caret_tmp)
  } else {
    cat("train caret on columns\n")
    cat(feature_cols)

    print(colnames(data))
    print("setdiff columns------ feature_cols - colnames(data)")
    print(setdiff(feature_cols, colnames(data)))
    print("------")

    print("Feature cols")
    print(feature_cols)
    print("------")

    n_cores <- min(parallel::detectCores() - 1, opts$threads)
    cl <- makeCluster(n_cores)
    registerDoParallel(cl)

    cat(paste0("Execute in parallel with ", n_cores, " cores"))

    caret_fit <- caret::train(
      x = data[, feature_cols],
      y = data[["truth"]],
      method = "rf",
      trControl = trControl,
      ntree = opts$ntree,
      tuneGrid = tuneGrid,
      allowParallel = TRUE
    )
    stopCluster(cl)
    registerDoSEQ()
    saveRDS(caret_fit, caret_tmp)
    print("caret_tmp saved")
  }




  model <- caret_fit$finalModel
  chosen_mtry <- model$mtry
  cat(paste("Mtry:", chosen_mtry))

  vi_caret <- varImp(model)
  vi_ordered <- vi_caret %>% 
    arrange(desc(Overall))

  chosen_cols <- rownames(vi_ordered)[seq_len(chosen_mtry)]

  # Retrain on chosen columns with randomForest
  train_df <- data[, c("truth", chosen_cols)]

  # split into train/test
  n <- nrow(train_df)
  test_size <- max(1, floor(opts$test_fraction * n))
  idx <- sample.int(n = n, size = test_size, replace = FALSE)
  test_data <- train_df[idx, , drop = FALSE]
  train_data <- train_df[-idx, , drop = FALSE]

  print("Train forest")
  rf <- randomForest(
    formula = truth ~ .,
    data = train_data,
    ntree = opts$ntree,
    maxnodes = if (opts$maxnodes > 0) opts$maxnodes else NULL,
    importance = TRUE
  )

  prefix <- opts$output_prefix
  pmml_path <- paste0(prefix, ".xml")
  varimp_path <- paste0(prefix, ".varimp.tsv")
  rds_path <- paste0(prefix, ".rds")
  varimp_png <- paste0(prefix, ".varimp.png")

  save_forest(rf, pmml_path)

  vi <- importance(rf, type = 1, scale = TRUE)
  vi_df <- data.frame(feature = rownames(vi), importance = vi[, 1], stringsAsFactors = FALSE, row.names = NULL)
  vi_df <- vi_df[order(vi_df$importance, decreasing = TRUE), ]
  write.table(vi_df, file = varimp_path, sep = "\t", quote = FALSE, row.names = FALSE)

  png(filename = varimp_png, width = 1200, height = 800)
  barplot(
    vi_df$importance,
    names.arg = vi_df$feature,
    las = 2,
    main = "Variable Importance",
    cex.names = 0.8
  )
  dev.off()

  saveRDS(rf, file = rds_path)

  evals <- evaluate_split(rf, test_data, train_data)
  message(sprintf("Split: train=%d test=%d (test fraction=%.3f)", nrow(train_data), nrow(test_data), opts$test_fraction))
  message(sprintf("Train   sens=%.4f prec=%.4f f1=%.4f", evals$train$sensitivity, evals$train$precision, evals$train$f1))
  message(sprintf("Test    sens=%.4f prec=%.4f f1=%.4f", evals$test$sensitivity, evals$test$precision, evals$test$f1))
  message("Chosen mtry (caret): ", chosen_mtry)
  message("Saved PMML: ", pmml_path)
  message("Saved varimp: ", varimp_path)
  message("Saved varimp plot: ", varimp_png)
  message("Saved model RDS: ", rds_path)
}

main()
