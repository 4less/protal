library(dplyr)
library(randomForest) 
library(ranger)
library(Metrics)
library(caret)
library(rfviz)
library(pmml)



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


data_path <- "/media/fritsche/Extreme_SSD/results/CAMI/HumanToy/protal_r214/profiles/profile_forest/truths.csv"
data <- read.csv(data_path, header=FALSE, sep='\t')
colnames(data) <- c("environment", "sample", "truth", "prediction", "taxon", "present_genes", "total_hits", "unique_hits", "mean_ani", "expected_gene_presence", "expected_gene_presence_ratio", "uniqueness", "mean_mapq", "vcov")

data$truth_raw <- as.logical(data$truth)
data$prediction <- as.logical(data$prediction)
data$truth <- as.factor(data$truth_raw)

set.seed(1234) # Set Seed so that same sample can be reproduced in future also
# Now Selecting 75% of data as sample from total 'n' rows of the data  

columns <- c("truth", "present_genes", "total_hits", "unique_hits", "mean_ani", "expected_gene_presence", "mean_mapq")
predictors <- c("present_genes", "total_hits", "unique_hits", "mean_ani", "expected_gene_presence", "mean_mapq")
sample <- sample.int(n = nrow(data), size = floor(.8*nrow(data)), replace = F)
training_data <- data[sample, columns]
test_data  <- data[-sample, columns]


rf <- randomForest(formula = truth ~ ., 
                   data = training_data,
                   ntree = 500)

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