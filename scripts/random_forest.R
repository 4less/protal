library(dplyr)
library(randomForest) 
library(ranger)
library(Metrics)
library(caret)
library(rfviz)
library(pmml)
library(tuneRanger)
# install.packages("tuneRanger")

save_forest <- function(model, file) {
  pmod <- pmml(model, file)
  save_pmml(pmod, file)
}

sensitivity <- function(ct) {
  tp <- ct[2, 2]
  tn <- ct[1, 1]
  fp <- ct[2, 1]
  fn <- ct[1, 2]
  
  return(tp / (tp+fn))
}

precision <- function(ct) {
  tp <- ct[2, 2]
  tn <- ct[1, 1]
  fp <- ct[2, 1]
  fn <- ct[1, 2]
  
  return(tp / (tp+fp))
}

evaluate <- function(rf, test.data, training.data, col.names=NA) {
  predictions.test <- predict(rf,  newdata=test.data)
  predictions.train <- predict(rf,  newdata=training.data)
  ct.test <- table(predictions.test, test.data$truth)
  ct.train <- table(predictions.train, training.data$truth)
  
  if (is.na(col.names[1])) {
    col.names = c("Test:  ", "Training:  ")
  }
  
  print(paste(col.names[2], "Sensitivity: ", sensitivity(ct.train), " Precision: ", precision(ct.train), sep=''))
  print(paste(col.names[1], "Sensitivity: ", sensitivity(ct.test), " Precision: ", precision(ct.test), sep=''))
  
  test.data$prediction <- predictions.test
  training.data$prediction <- predictions.train
  
  return(rbind(test.data, training.data))
}

evaluate_ranger <- function(rf, test.data, training.data) {
  predictions.test <- predict(rf,  test.data)$predictions
  predictions.train <- predict(rf,  training.data)$predictions
  ct.test <- table(predictions.test, test.data$truth)
  ct.train <- table(predictions.train, training.data$truth)
  
  print(paste("TRAINING: Sensitivity: ", sensitivity(ct.train), " Precision: ", precision(ct.train), sep=''))
  print(paste("TEST:     Sensitivity: ", sensitivity(ct.test), " Precision: ", precision(ct.test), sep=''))
}

evaluate_tuned <- function(rf, test.data, training.data) {
  predictions.test <- predict(rf,  newdata=test.data)
  predictions.train <- predict(rf,  newdata=training.data)
  ct.test <- table(predictions.test$data$response, test.data$truth)
  ct.train <- table(predictions.train$data$response, training.data$truth)
  
  print(paste("TRAINING: Sensitivity: ", sensitivity(ct.train), " Precision: ", precision(ct.train), sep=''))
  print(paste("TEST:     Sensitivity: ", sensitivity(ct.test), " Precision: ", precision(ct.test), sep=''))
}
set.seed(1234)






data_path <- "/usr/users/QIB_fr017/fritsche/CLionProjects/protal/scripts/data/truths.csv"
data_path2 <- "/usr/users/QIB_fr017/fritsche/CLionProjects/protal/scripts/data/truth2.tsv"
data_path3 <- "/usr/users/QIB_fr017/fritsche/CLionProjects/protal/scripts/data/with_alleles/all.tsv"
data_path4 <- "/usr/users/QIB_fr017/fritsche/CLionProjects/protal/scripts/data/with_mask/all.tsv"
data_path5 <- "/usr/users/QIB_fr017/fritsche/CLionProjects/protal/scripts/data/forest_data/unique_truths.tsv"

data <- read.csv(data_path, header=FALSE, sep='\t')
colnames(data) <- c("environment", "sample", "truth", "prediction", "taxon", "present_genes", "total_hits", "unique_hits", "mean_ani", "expected_gene_presence", "expected_gene_presence_ratio", "uniqueness", "mean_mapq", "vcov")

paste("Training/Test  ", nrow(training_data3), "/", nrow(test_data3), sep='')
data2 <- read.csv(data_path2, header=FALSE, sep='\t')
colnames(data2) <- c("truth", "prediction", "taxon", "present_genes", "total_hits", "unique_hits", "mean_ani", "expected_gene_presence", "expected_gene_presence_ratio", "uniqueness", "mean_mapq", "vcov")

data3 <- read.csv(data_path3, header=FALSE, sep='\t')
colnames(data3) <- c("truth", "prediction", "taxon", "present_genes", "total_hits", "unique_hits", "mean_ani", "expected_gene_presence", "expected_gene_presence_ratio", "uniqueness", "mean_mapq", "variance1", "variance2", "A0", "A1", "A2" , "A3", "A4", "AF0", "AF1", "AF2", "AF3", "AF4","stddev")

data4 <- read.csv(data_path4, header=FALSE, sep='\t')
colnames(data4) <- c("truth", "prediction", "taxon", "present_genes", "total_hits", "unique_hits", "mean_ani", "expected_gene_presence", "expected_gene_presence_ratio", "uniqueness", "mean_mapq", "variance1", "variance2", "A0", "A1", "A2" , "A3", "A4", "AF0", "AF1", "AF2", "AF3", "AF4","stddev", "hittable")

data5 <- read.csv(data_path5, header=FALSE, sep='\t')
colnames(data5) <- c("truth", "prediction", "taxon", "present_genes", "total_hits", "unique_hits", "mean_ani", "expected_gene_presence", "expected_gene_presence_ratio", "uniqueness", "mean_mapq", "variance1", "variance2", "A0", "A1", "A2" , "A3", "A4", "AF0", "AF1", "AF2", "AF3", "AF4","stddev", "hittable", "lu", "lu_genes", "lsu", "lsu_genes", "su_genome", "lu_genome", "lsu_genome", "total_genome", "lu_rate", "lsu_rate", "lu_gene_rate", "lsu_gene_rate", "lu_gene_rate2", "lsu_gene_rate2")


data$truth_raw <- as.logical(data$truth)
data$prediction <- as.logical(data$prediction)
data$truth <- as.factor(data$truth_raw)

data2$truth_raw <- as.logical(data2$truth)
data2$prediction <- as.logical(data2$prediction)
data2$truth <- as.factor(data2$truth_raw)

data3$truth_raw <- as.logical(data3$truth)
data3$prediction <- as.logical(data3$prediction)
data3$truth <- as.factor(data3$truth_raw)

data4$truth_raw <- as.logical(data4$truth)
data4$prediction <- as.logical(data4$prediction)
data4$truth <- as.factor(data4$truth_raw)

data5$truth_raw <- as.logical(data5$truth)
data5$prediction <- as.logical(data5$prediction)
data5$truth <- as.factor(data5$truth_raw)

combined <- rbind(data[, colnames(data2)], data2)
 # Set Seed so that same sample can be reproduced in future also
# Now Selecting 75% of data as sample from total 'n' rows of the data  

columns <- c("truth", "present_genes", "total_hits", "unique_hits", "mean_ani", "expected_gene_presence", "mean_mapq")
predictors <- c("present_genes", "total_hits", "unique_hits", "mean_ani", "expected_gene_presence", "mean_mapq")

# columns2 <- c("truth", "present_genes", "total_hits", "unique_hits", "mean_ani", "expected_gene_presence", "mean_mapq", "variance1", "variance2", "stddev")
# predictors2 <- c("present_genes", "total_hits", "unique_hits", "mean_ani", "expected_gene_presence", "mean_mapq", "variance1", "variance2", "stddev")
# columns2 <- c("truth", "present_genes", "total_hits", "unique_hits", "mean_ani", "expected_gene_presence", "variance1", "stddev", "A1", "A2", "AF0")
columns2 <- c("truth", "present_genes", "total_hits", "unique_hits", "mean_ani", "expected_gene_presence", "variance1", "variance2", "stddev", "A0", "A1", "A2", "A3", "A4", "AF0", "AF1", "AF2", "AF3", "AF4")
columns2 <- c("truth", "present_genes", "expected_gene_presence", "total_hits", "unique_hits", "A2", "AF1", "AF0", "A1", "mean_ani", "variance1", "stddev")
columns3 <- c("truth", "present_genes", "expected_gene_presence", "total_hits", "unique_hits", "A2", "AF1", "AF0", "A1", "mean_ani", "variance1", "stddev", "hittable")

predictors2 <- c("present_genes", "total_hits", "unique_hits", "mean_ani", "expected_gene_presence", "variance1", "stddev")

sample <- sample.int(n = nrow(data), size = floor(.8*nrow(data)), replace = F)
training_data <- data[sample, columns]
test_data  <- data[-sample, columns]
paste("Training/Test  ", nrow(training_data), "/", nrow(test_data), sep='')

sample2 <- sample.int(n = nrow(combined), size = floor(.8*nrow(combined)), replace = F)
training_data2 <- combined[sample2, columns]
test_data2 <- combined[-sample2, columns]
paste("Training/Test  ", nrow(training_data2), "/", nrow(test_data2), sep='')

sample3 <- sample.int(n = nrow(data3), size = floor(.8*nrow(data3)), replace = F)
training_data3 <- data3[sample3, columns2]
test_data3 <- data3[-sample3, columns2]
paste("Training/Test  ", nrow(training_data3), "/", nrow(test_data3), sep='')


# Ranger model with parameter tuning
train.task <- makeClassifTask(data = training_data, target = "truth")
estimateTimeTuneRanger(train.task)
res = tuneRanger(train.task, measure = list(multiclass.brier), num.trees = 500,
                 num.threads = 8, iters = 150, save.file.path = NULL)
varImp(res$model)#(res$model, sort = TRUE , n.var , main = "Variable importance" )
str(res$model)
evaluate_tuned(res$model, test_data, training_data)

ranger.fr <- ranger(truth ~ present_genes + total_hits + unique_hits + mean_ani + expected_gene_presence + mean_mapq, 
                    data = training_data, 
                    importance = 'permutation',
                    scale.permutation.importance = TRUE,
                    mtry = res$recommended.pars$mtry,
                    min.node.size = res$recommended.pars$min.node.size,
                    sample.fraction = res$recommended.pars$sample.fraction,
                    multiclass.brier = res$recommended.pars$multiclass.brier)
ranger.fr$variable.importance
evaluate_ranger(ranger.fr, test_data, training_data)

##################################################################################################


# Random forest ##################################################################################
rf <- randomForest(formula = truth ~ ., 
                   data = training_data,
                   ntree = 1000, maxnodes=15)
evaluate(rf, test_data, training_data)

rf.combined <- randomForest(formula = truth ~ ., 
                   data = training_data2,
                   ntree = 1000, maxnodes=30)
evaluate(rf.combined, test_data2, training_data2)

#################################################################################################################################
# Divide into training and test
sample3 <- sample.int(n = nrow(data3), size = floor(.8*nrow(data3)), replace = F)
training_data3 <- data3[sample3, columns2]
test_data3 <- data3[-sample3, columns2]
paste("Training/Test  ", nrow(training_data3), "/", nrow(test_data3), sep='')

# Test
rf.var <- randomForest(formula = truth ~ ., 
                            data = training_data3,
                            ntree = 256, maxnodes=128)
evaluate(rf.var, test_data3, training_data3)

#######################

fill_in_allele_info <- function(data4) {
  data4$total_a <- data4$A0 + data4$A1 + data4$A2 + data4$A3 + data4$A4
  data4$total_af <- data4$AF0 + data4$AF1 + data4$AF2 + data4$AF3 + data4$AF4
  data4$RA0 <- data4$A0 / data4$total_a
  data4$RA1 <- data4$A1 / data4$total_a
  data4$RA2 <- data4$A2 / data4$total_a
  data4$RA3 <- data4$A3 / data4$total_a
  data4$RA4 <- data4$A4 / data4$total_a
  data4$RAF0 <- data4$AF0 / data4$total_af
  data4$RAF1 <- data4$AF1 / data4$total_af
  data4$RAF2 <- data4$AF2 / data4$total_af
  data4$RAF3 <- data4$AF3 / data4$total_af
  data4$RAF4 <- data4$AF4 / data4$total_af
  data4$RA0[is.nan(data4$RA0)] <- 0
  data4$RA1[is.nan(data4$RA1)] <- 0
  data4$RA2[is.nan(data4$RA2)] <- 0
  data4$RA3[is.nan(data4$RA3)] <- 0
  data4$RA4[is.nan(data4$RA4)] <- 0
  data4$RAF0[is.nan(data4$RAF0)] <- 0
  data4$RAF1[is.nan(data4$RAF1)] <- 0
  data4$RAF2[is.nan(data4$RAF2)] <- 0
  data4$RAF3[is.nan(data4$RAF3)] <- 0
  data4$RAF4[is.nan(data4$RAF4)] <- 0
  return(data4)
}
data4 <- fill_in_allele_info(data4)
data5 <- fill_in_allele_info(data5)

data5$lsu_per_read <- data5$lsu / data5$total_hits
data5$lu_per_read <- data5$lu / data5$total_hits
data5[is.infinite(data5$lsu_gene_rate2),"lsu_gene_rate2"] <- 0
data5[is.infinite(data5$lu_rate),"lu_rate"] <- 0
data5[is.infinite(data5$lsu_rate),"lsu_rate"] <- 0

split.tt <- function(data, columns) {
  sample <- sample.int(n = nrow(data), size = floor(.8*nrow(data)), replace = F)
  result <- list()
  result$training <- data[sample, columns]
  result$test <- data[-sample, columns]
  paste("Training/Test  ", nrow(result$training), "/", nrow(result$test), sep='')
  return(result)
}

columns3 <- c("truth", "present_genes", "expected_gene_presence", "total_hits", "unique_hits", "A2", "AF1", "AF0", "A1", "mean_ani", "variance1", "stddev")
columns3 <- c("truth", "present_genes", "expected_gene_presence", "total_hits", "unique_hits", "RA3", "RAF0", "RAF1", "mean_ani", "variance1", "stddev")#, "RA2")
columns3 <- c("truth", "present_genes", "total_hits", "unique_hits", "mean_ani", "expected_gene_presence", "variance1", "variance2", "stddev", "A0", "A1", "A2", "A3", "A4", "AF0", "AF1", "AF2", "AF3", "AF4", "RA0", "RA1", "RA2", "RA3", "RA4", "RAF0", "RAF1", "RAF2", "RAF3", "RAF4")
columns3 <- c("truth", "present_genes", "expected_gene_presence", "total_hits", "unique_hits", "RA3", "RAF0", "RAF1", "mean_ani", "variance1", "stddev")#, "RA2")

columns5 <- c("truth", "present_genes", "expected_gene_presence", "total_hits", "unique_hits", "RA3", "RAF0", "RAF1", "mean_ani", "variance1", "stddev", "lsu_per_read", "lu_per_read", "lu_genome", "lsu_genome", "lu", "lsu", "lsu_genes", "lu_genes")#, "RA2")

split_train_forest <- function(data, columns, ntree, maxnodes) {
  result <- list()
  data.tt <- split.tt(data, columns)
  result$data <- data.tt
  result$rf <- randomForest(formula = truth ~ ., 
                         data = data.tt$training,
                         ntree = ntree, maxnodes=maxnodes)
  result$eval <- evaluate(result$rf , data.tt$test, data.tt$training)
  result$varImportance <- varImpPlot(result$rf , sort = TRUE , length(columns), main = "Variable importance" )
  return(result)
}

data5 %>% filter(truth == "TRUE" | prediction == "TRUE") %>%
  ggplot(aes(x = truth, y=total_hits, fill=truth)) +
  geom_boxplot()
data5 %>% filter(truth == "TRUE" | prediction == "TRUE") %>%
  ggplot(aes(x = truth, y=lsu_per_read, fill=truth)) +
  geom_boxplot()
data5 %>% filter(truth == "TRUE" | prediction == "TRUE") %>%
  ggplot(aes(x = truth, y=lu_per_read, fill=truth)) +
  geom_boxplot()

# split into training and test
data.tt <- split.tt(data4, columns3)
data.tt <- split.tt(data5, columns5)


rf.result <- split_train_forest(data4, columns3, 256, 64)

columns5 <- c("truth", "present_genes", "expected_gene_presence", "lu", "lsu_genes", "lu_genes", "lu_gene_rate", "lsu_gene_rate", "lu_gene_rate2")#, "lsu_gene_rate2")#, "RA2")

rf.unique.result <- split_train_forest(data5, columns5, 32, 32)
varImpPlot(rf.unique.result$rf, sort = TRUE , length(columns5), main = "Variable importance" )

columns6 <- c("truth", "present_genes", "total_hits", "unique_hits", "mean_ani", "expected_gene_presence", "expected_gene_presence_ratio", "uniqueness", "mean_mapq", "variance1", "variance2", "A0", "A1", "A2" , "A3", "A4", "AF0", "AF1", "AF2", "AF3", "AF4","stddev", "hittable", "lu", "lu_genes", "lsu", "lsu_genes", "su_genome", "lu_genome", "lsu_genome", "total_genome", "lu_rate", "lsu_rate", "lu_gene_rate", "lsu_gene_rate", "lu_gene_rate2", "lsu_per_read", "lu_per_read", "RA0", "RA1", "RA2", "RA3", "RA4", "RAF0", "RAF1", "RAF2", "RAF3", "RAF4")
rf.unique.result <- split_train_forest(data5, columns6, 512, 128)
varImpPlot(rf.unique.result$rf, sort = TRUE , length(columns6), main = "Variable importance" )

#"total_hits", "unique_hits","A0", "A1", "A2" , "A3", "A4", , "AF1", "AF2", "AF3", "AF4""RA0", "RA1", "RA2", "RA3", "RA4", "RAF0", "RAF1", "RAF2", "RAF3", "RAF4"
columns6 <- c("truth", "present_genes", "mean_ani", "expected_gene_presence", "expected_gene_presence_ratio", "uniqueness", "mean_mapq", "variance1", "variance2", "AF0", "RAF0", "stddev", "hittable", "lu", "lu_genes", "lsu", "lsu_genes", "su_genome", "lu_genome", "lsu_genome", "total_genome", "lu_rate", "lsu_rate", "lu_gene_rate", "lsu_gene_rate", "lu_gene_rate2", "lsu_per_read", "lu_per_read")
rf.unique.result <- split_train_forest(data5, columns6, 512, 128)
varImpPlot(rf.unique.result$rf, sort = TRUE , length(columns6), main = "Variable importance" )

#"total_hits", "unique_hits","A0", "A1", "A2" , "A3", "A4", , "AF1", "AF2", "AF3", "AF4""RA0", "RA1", "RA2", "RA3", "RA4", "RAF0", "RAF1", "RAF2", "RAF3", "RAF4"
columns6 <- c("truth", "present_genes", "mean_ani", "expected_gene_presence", "expected_gene_presence_ratio", "uniqueness", "mean_mapq", "variance1", "variance2", "AF0", "RAF0", "stddev", "hittable", "lu", "lu_genes", "lsu", "lsu_genes", "su_genome", "lu_genome", "lsu_genome", "total_genome", "lu_rate", "lsu_rate", "lu_gene_rate", "lsu_gene_rate", "lu_gene_rate2", "lsu_per_read", "lu_per_read")
rf.unique.result <- split_train_forest(data5, columns6, 256, 64)
varImpPlot(rf.unique.result$rf, sort = TRUE , length(columns6), main = "Variable importance" )

save_forest(rf.unique.result$rf, "/usr/users/QIB_fr017/fritsche/CLionProjects/protal/scripts/forests/random_forest_unique_256_64.xml")
write.csv(columns6,file="/usr/users/QIB_fr017/fritsche/CLionProjects/protal/scripts/forests/random_forest_unique_256_64.columns",row.names=F)

#"total_hits", "unique_hits","A0", "A1", "A2" , "A3", "A4", , "AF1", "AF2", "AF3", "AF4""RA0", "RA1", "RA2", "RA3", "RA4", "RAF0", "RAF1", "RAF2", "RAF3", "RAF4"
columns6 <- c("truth", "present_genes", "mean_ani", "expected_gene_presence", "expected_gene_presence_ratio", "uniqueness", "mean_mapq", "variance1", "variance2", "AF0", "RAF0", "stddev", "hittable", "lu", "lu_genes", "lsu", "lsu_genes", "su_genome", "lu_genome", "lsu_genome", "total_genome", "lu_rate", "lsu_rate", "lu_gene_rate", "lsu_gene_rate", "lu_gene_rate2", "lsu_per_read", "lu_per_read")
rf.unique.result <- split_train_forest(data5, columns6, 128, 32)
varImpPlot(rf.unique.result$rf, sort = TRUE , length(columns6), main = "Variable importance" )

rf.unique.result$eval$correct <- rf.unique.result$eval$truth == rf.unique.result$eval$prediction
View(rf.unique.result$eval)

# Test
rf.hit <- randomForest(formula = truth ~ .,
                       data = data.tt$training,
                       ntree = 256, maxnodes=64)
evaluate(rf.hit, data.tt$test, data.tt$training)
varImpPlot(rf.hit, sort = TRUE , length(columns3), main = "Variable importance" )

varImpPlot(rf.var, sort = TRUE , 11, main = "Variable importance" )



save_forest(rf.hit, "/usr/users/QIB_fr017/fritsche/CLionProjects/protal/scripts/forests/random_forest_unique.xml")


pmod <- pmml(rf.hit)
save_pmml(pmod, "/media/fritsche/Extreme_SSD/db/protal/protal_flex16_r214/random_forest_hit.xml")

###
varImpPlot(rf.combined, sort = TRUE , 6, main = "Variable importance" )
varImpPlot(rf.var, sort = TRUE , 11, main = "Variable importance" )


evaluate(rf.combined, test_data2, training_data2)
evaluate(rf.var, test_data3, training_data3)

cat("Small model")
evaluate(rf, test_data2, test_data, col.names=c("newtest:  ", "oldtest:  "))
cat("Larger model")
evaluate(rf.combined, test_data2, test_data, col.names=c("newtest:  ", "oldtest:  "))


cat("Small model")
evaluate(rf, data2, data, col.names=c("newtest:  ", "oldtest:  "))
cat("Larger model")
evaluate(rf.combined, data2, data, col.names=c("newtest:  ", "oldtest:  "))

save_forest(rf.combined, "/usr/users/QIB_fr017/fritsche/CLionProjects/protal/scripts/forests/random_forest_1000_n15_all.xml")
save_forest(rf.var, "/usr/users/QIB_fr017/fritsche/CLionProjects/protal/scripts/forests/random_forest_1000_n30_all_var_alleles.xml")

varImpPlot(rf, sort = TRUE , 6, main = "Variable importance" )

evaluate(rf, test_data, training_data)

node_sizes <- seq(10,200,10)

rf.list <- lapply(node_sizes, FUN=function(x) {randomForest(formula = truth ~ ., 
                                                                 data = training_data,
                                                                 ntree = 500, maxnodes = x)})

# Random forest ##################################################################################

rf.restricted <- randomForest(formula = truth ~ ., 
                    data = training_data,
                   ntree = 1000, maxnodes = 15)

rf.prep <- rf_prep(x=training_data[,predictors], y=training_data$truth, ntree = 1000, maxnodes = 15)


predictions.rf2 <- predict(rf.restricted,  test_data)
predictions.rf2.all <- predict(rf.restricted,  data[,columns])

predictions.rf3 <- predict(rf.prep$rf,  test_data)
predictions.rf3.all <- predict(rf.prep$rf,  data[,columns])

ct1 <- table(predictions.rf, test_data$truth)
ct2 <- table(predictions.rf2, test_data$truth)
ct3 <- table(predictions.rf3, test_data$truth)

sensitivity(ct1)
precision(ct1)
sensitivity(ct3)
precision(ct3)

sensitivity(ct2)
precision(ct2)

c1 <- table(predictions.rf2, test_data$truth)
c1.a <- table(predictions.rf2.all, data[,columns]$truth)
c2 <-table(predictions.rf3, test_data$truth)
c2.a <- table(predictions.rf3.all, data[,columns]$truth)
sensitivity(c1)
precision(c1)
sensitivity(c1.a)
precision(c1.a)
sensitivity(c2)
precision(c2)
sensitivity(c2.a)
precision(c2.a)

varImpPlot(rf.prep$rf)

x <- node_sizes
y.sens <- lapply(rf.list, FUN=function(rf) {
  predictions.rf <- predict(rf,  test_data) > 0.5
  ct <- table(predictions.rf, test_data$truth)
  return(sensitivity(ct))
})
y.prec <- lapply(rf.list, FUN=function(rf) {
  predictions.rf <- predict(rf,  test_data) > 0.5
  ct <- table(predictions.rf, test_data$truth)
  return(precision(ct))
})

plot(x, y.sens)
plot(x, y.prec)

bcrf <- rf_viz(rf.prep$rf)#, input=TRUE, imp=TRUE, cmd=TRUE)
bcrf
#> sensitivity(ct1)
#[1] 0.9665428
#> precision(ct1)
#[1] 0.9942639
#> sensitivity(ct2)
#[1] 0.9535316
#> precision(ct2)
#[1] 0.994186


#> sensitivity(ct1)
#[1] 0.9582505
#> precision(ct1)
#[1] 0.9816701
#> sensitivity(ct2)
#[1] 0.9423459
#> precision(ct2)
#[1] 0.9773196

rf.restricted
pmod <- pmml(rf.restricted)
save_pmml(pmod, "/media/fritsche/Extreme_SSD/db/protal/protal_flex16_r214/random_forest_1000_n15.xml")



rf.prep$rf
pmod <- pmml(rf.prep$rf)
save_pmml(pmod, "/media/fritsche/Extreme_SSD/db/protal/protal_flex16_r214/random_forest2.xml")