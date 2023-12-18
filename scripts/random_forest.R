library(dplyr)
library(randomForest) 
library(ranger)
library(Metrics)
library(caret)
library(rfviz)
library(pmml)
library(tuneRanger)
install.packages("tuneRanger")

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

evaluate <- function(rf, test.data, training.data) {
  predictions.test <- predict(rf,  newdata=test.data)
  predictions.train <- predict(rf,  newdata=training.data)
  ct.test <- table(predictions.test, test.data$truth)
  ct.train <- table(predictions.train, training.data$truth)
  
  print(paste("TRAINING: Sensitivity: ", sensitivity(ct.train), " Precision: ", precision(ct.train), sep=''))
  print(paste("TEST:     Sensitivity: ", sensitivity(ct.test), " Precision: ", precision(ct.test), sep=''))
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

# Random forest
rf <- randomForest(formula = truth ~ ., 
                   data = training_data,
                   ntree = 500)
evaluate(rf, test_data, training_data)
varImpPlot(rf, sort = TRUE , 6, main = "Variable importance" )






data_path <- "/usr/users/QIB_fr017/fritsche/CLionProjects/protal/scripts/data/truths.csv"
data <- read.csv(data_path, header=FALSE, sep='\t')
colnames(data) <- c("environment", "sample", "truth", "prediction", "taxon", "present_genes", "total_hits", "unique_hits", "mean_ani", "expected_gene_presence", "expected_gene_presence_ratio", "uniqueness", "mean_mapq", "vcov")

data$truth_raw <- as.logical(data$truth)
data$prediction <- as.logical(data$prediction)
data$truth <- as.factor(data$truth_raw)

 # Set Seed so that same sample can be reproduced in future also
# Now Selecting 75% of data as sample from total 'n' rows of the data  

columns <- c("truth", "present_genes", "total_hits", "unique_hits", "mean_ani", "expected_gene_presence", "mean_mapq")
predictors <- c("present_genes", "total_hits", "unique_hits", "mean_ani", "expected_gene_presence", "mean_mapq")
sample <- sample.int(n = nrow(data), size = floor(.8*nrow(data)), replace = F)
training_data <- data[sample, columns]
test_data  <- data[-sample, columns]
paste("Training/Test  ", nrow(training_data), "/", nrow(test_data), sep='')



evaluate(rf, test_data, training_data)

node_sizes <- seq(10,200,10)

rf.list <- lapply(node_sizes, FUN=function(x) {randomForest(formula = truth ~ ., 
                                                                 data = training_data,
                                                                 ntree = 500, maxnodes = x)})

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