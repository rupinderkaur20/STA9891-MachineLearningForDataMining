# STA9891: Final Project on Binary Classfication
# Group 8, Members: Rupinder Kaur, Janani RaviChandran


rm(list = ls()) #delete objects
cat("\014")

library(Matrix)
library(tidyverse)
library(randomForest)
library(glmnet)
library(dplyr)
library(pROC)
library(gridExtra)
library(tibble)

############################### Data Preparation ##################################

codon = read.csv("codon_usage.csv", header = TRUE) # num rows: 13,028

## UUU and UUC are char not num
str(codon)
codon$UUU = as.numeric(codon$UUU)
codon$UUC = as.numeric(codon$UUC)

sum(is.na(codon))
codon = na.omit(codon) # num rows: 13,026

# there are 11 unique Kingdoms
nrow(distinct(codon,Kingdom))
unique(codon$Kingdom)

# Kingdoms: 'vrl' and 'pln'     
codon = filter(codon, Kingdom == 'vrl' | Kingdom == 'pln') # num rows: 5354

# How many observations for each class
# num rows of kingdom 'vrl': 2831
vrl = filter(codon, Kingdom == 'vrl')
n_rows_vrl = dim(vrl)[1]

# num rows for kingdom 'pln': 2523
pln = filter(codon, Kingdom == 'pln')
n_rows_pln = dim(pln)[1]

codon$Kingdom = replace(codon$Kingdom, codon$Kingdom == 'vrl', 0)
codon$Kingdom = replace(codon$Kingdom, codon$Kingdom == 'pln', 1)

# final dataset: with Kingdom and codons only
X = codon[, -c(1,2,3,4,5)]
Y = as.numeric(codon[,1])

n = dim(X)[1]
p = dim(X)[2] # 64

########################## Model Fitting on splitted data ##########################
k = 50 # iteration for loop

# Matrix to store AUC values
auc.df = data.frame(matrix(nrow = 0, ncol = 4))
names(auc.df) = c("loop", "model", "auc.train", "auc.test")

for (i in c(1:k)) {
  cat("Loop", i, "\n")
  loop = i
  
  # splitting data into training and test sets
  i.mix = sample(n) #shuffling
  train = i.mix[1:floor(n*0.65)]
  
  #Training predictors and response
  X.train = X[train,]
  Y.train = Y[train]
  
  # Validation predictors and response
  X.test = X[-train,]
  Y.test = Y[-train]
  
  # Implementing Random forest, logistic elastic-net, logistic lasso, and
  # logistic ridge
  
  # Random Forest
  rf.fit = randomForest(x = X.train, y = as.factor(Y.train), mtry = sqrt(p))

  # training AUC 
  rf.pred.train = predict(rf.fit, X.train, type="prob", norm.votes=TRUE)
  rf.pred.train = data.frame(rf.pred.train)
  rf.train.roc = roc(Y.train, rf.pred.train$X1)
  auc.train = pROC::auc(rf.train.roc)
  
  # Test AUC
  rf.pred.test = predict(rf.fit, X.test, type="prob", norm.votes=TRUE)
  rf.pred.test = data.frame(rf.pred.test)
  rf.test.roc = roc(Y.test, rf.pred.test$X1)
  auc.test = pROC::auc(rf.test.roc)
  
  model = "Random Forest"
  
  # Storing AUC for Random Forest
  auc.df = rbind(auc.df, cbind(loop, model, auc.train, auc.test))
  
  # Logistic Lasso, Ridge and Elastic-net
  alphas = c(1, 0.5, 0)
  for (a in alphas){
    if (a==1){ 
      print("LASSO")
      model = "Lasso"}
    
    if (a==0.5){ 
      print("ELASTIC-NET")
      model = "Elastic-Net"}
    
    if (a==0){ 
      print("RIDGE")
      model = "Ridge"}
    
    # Tuning the lambda values and recording the time
    cv.start = Sys.time()
    cv.fit = cv.glmnet(as.matrix(X.train), Y.train, family = "binomial", alpha = a, 
                       type.measure = "auc", nfolds = 10)
    cv.time = round(Sys.time() - cv.start, 1)
    
    # cross validation curves for the 45th sample
    if (i == 45){
      plot(cv.fit)
      if (model == "Lasso"){
        title(paste("Lasso, ", "CV Time: ", cv.time, units(cv.time)), line= 2.5)
      } else if (model == "Elastic-Net"){
        title(paste("Elastic Net, ", "CV Time: ", cv.time, units(cv.time)), line= 2.5)
      } else {
        title(paste("Ridge, ", "CV Time: ", cv.time, units(cv.time)), line= 2.5)
      }
    }
    
    # fitting the model
    fit = glmnet(as.matrix(X.train), Y.train, family = "binomial", alpha = a, 
                 lambda = cv.fit$lambda.min )
    beta0.hat = fit$a0
    beta.hat = as.vector(fit$beta)
    
    prob.train = exp(as.matrix(X.train) %*% beta.hat +  beta0.hat  )/(1 + exp(as.matrix(X.train) %*% beta.hat +  beta0.hat))
    prob.test = exp(as.matrix(X.test) %*% beta.hat +  beta0.hat  )/(1 + exp(as.matrix(X.test) %*% beta.hat +  beta0.hat  ))

    thrs = seq(0, 1, 0.01)
    thrs.length = length(thrs)
    
    ## To store FPR AND TPR
    FPR.train = matrix(0, thrs.length)
    TPR.train = matrix(0, thrs.length)
    FPR.test = matrix(0, thrs.length)
    TPR.test = matrix(0, thrs.length)
    
    for (j in c(1:thrs.length)){
      # calculate the FPR and TPR for train data 
      Y.hat.train = ifelse(prob.train > thrs[j], 1, 0) 
      FP.train = sum(Y.train[Y.hat.train==1] == 0) # false positives 
      TP.train = sum(Y.hat.train[Y.train==1] == 1) # true positives  
      P.train = sum(Y.train==1) # total positives in the data
      N.train = sum(Y.train==0) # total negatives in the data
      FPR.train[j] = FP.train/N.train # false positive rate 
      TPR.train[j] = TP.train/P.train # true positive rate 
      
      # calculate the FPR and TPR for test data 
      Y.hat.test = ifelse(prob.test > thrs[j], 1, 0)
      FP.test = sum(Y.test[Y.hat.test==1] == 0)
      TP.test = sum(Y.hat.test[Y.test==1] == 1) 
      P.test = sum(Y.test==1) 
      N.test = sum(Y.test==0) 
      FPR.test[j] = FP.test/N.test
      TPR.test[j] = TP.test/P.test 
    }
    #Calculating AUC
    auc.train = sum((TPR.train[1:(thrs.length-1)] + 0.5 * diff(TPR.train)) * diff(FPR.train))
    auc.test = sum((TPR.test[1:(thrs.length-1)] + 0.5 * diff(TPR.test)) * diff(FPR.test))
    
    # storing AUC
    auc.df = rbind(auc.df, cbind(loop, model, auc.train, auc.test))
  }
}

################################## Boxplot of AUCs ######################################
# training AUC and test AUC are characters
str(auc.df) 
auc.df$auc.train = abs(as.numeric(auc.df$auc.train))
auc.df$auc.test = abs(as.numeric(auc.df$auc.test)) 

p1 <- auc.df %>%
  ggplot(mapping = aes(x = model, y = auc.train, fill = model))+
  geom_boxplot() + ylim(0.97, 1) + labs(x = "Models", y = "Training AUC") + 
  theme(legend.position = "none")

p2 <- auc.df %>%
  ggplot(mapping = aes(x = model, y = auc.test, fill = model))+
  geom_boxplot() + ylim(0.97, 1) + labs(x = "Models", y = "Test AUC") +
  theme(legend.position = "none")

grid.arrange(p1, p2, nrow = 1)

###################### Fitting Models on all data #############################
# logistic lasso
lasso.cv.start = Sys.time()
lasso.cv.fit = cv.glmnet(as.matrix(X), Y, family = "binomial", alpha = 1, 
                   type.measure = "auc", nfolds = 10)
lasso.fit = glmnet(as.matrix(X), Y, family = "binomial", alpha = 1, 
             lambda = lasso.cv.fit$lambda.min )
lasso.cv.time = round(Sys.time() - lasso.cv.start, 1)

# logistic Elastic Net
elnet.cv.start = Sys.time()
elnet.cv.fit = cv.glmnet(as.matrix(X), Y, family = "binomial", alpha = 0.5, 
                         type.measure = "auc", nfolds = 10)
elnet.fit = glmnet(as.matrix(X), Y, family = "binomial", alpha = 0.5, 
                   lambda = elnet.cv.fit$lambda.min )
elnet.cv.time = round(Sys.time() - elnet.cv.start, 1)

# logistic ridge
ridge.cv.start = Sys.time()
ridge.cv.fit = cv.glmnet(as.matrix(X), Y, family = "binomial", alpha = 0, 
                         type.measure = "auc", nfolds = 10)
ridge.fit = glmnet(as.matrix(X), Y, family = "binomial", alpha = 0, 
                   lambda = ridge.cv.fit$lambda.min )
ridge.cv.time = round(Sys.time() - ridge.cv.start, 1)

# Random Forest
rf.cv.start = Sys.time()
random.forest.fit = randomForest(x = X, y = as.factor(Y), mtry = sqrt(p))
rf.cv.time = round(Sys.time() - rf.cv.start, 1)

# median(AUC) for each model
# lasso
lasso = filter(auc.df,model == 'Lasso')
lasso.median = round(median(lasso$auc.test),7)

# Elastic Net
elnet = filter(auc.df, model == 'Elastic-Net')
elnet.median = round(median(elnet$auc.test),7)

# Ridge
ridge = filter(auc.df, model == 'Ridge')
ridge.median = round(median(ridge$auc.test),7)

# Random Forest
rf = filter(auc.df, model == 'Random Forest')
rf.median = round(median(rf$auc.test),7)

# time Vs. accuracy table
Model = c("Lasso", "Elastic Net", "Ridge", "Random Forest")
AUC.Median = c(lasso.median, elnet.median, ridge.median, rf.median)
Time = c(lasso.cv.time, elnet.cv.time, ridge.cv.time, rf.cv.time)

time.vs.acc = data.frame(Model, AUC.Median, Time)

################## Bar Plots of Standardized Coefficients #####################
beta.df = data.frame(Lasso = as.numeric(lasso.fit$beta), 
                     Elastic.Net = as.numeric(elnet.fit$beta),
                     Ridge = as.numeric(ridge.fit$beta),
                     Random.Forest = random.forest.fit$importance[,1])

beta.df <- rownames_to_column(beta.df, "Codons")

# Standardizing the coefficients & arranging the data
beta.df = beta.df %>% mutate_at(c("Lasso", "Elastic.Net", "Ridge"),
                      ~(scale(.) %>% as.vector)) %>%
  arrange(desc(Elastic.Net)) %>% 
  mutate(Order = seq(1,p,by=1))

# Sorting based on Order
beta.df = beta.df %>%
  mutate(Codons = fct_reorder(Codons, Order))

# Reshaping the data
beta.df = beta.df %>% 
  gather(Model, Coefficients, Lasso:Random.Forest) %>%
  mutate(pos = Coefficients >= 0)

# Plot the coefficients and importance
ggplot(beta.df) + 
  aes(x=Codons, y=Coefficients, fill = pos) + 
  geom_col(color = "black", show.legend = FALSE) + 
  facet_grid(Model ~., scales="free") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.text.x = element_text(size = 6, angle = 90, hjust = 1)) + 
  theme(axis.text.y = element_text(size = 6)) + 
  ylab("Coefficients/ Variable Importance")

