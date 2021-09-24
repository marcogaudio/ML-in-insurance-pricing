# models
source("utilities.R")
library(tidyverse)
#### GLM
# choose a number of fold (optimal number usually 5 or 10)
k <- 5 


# re-sample the data
set.seed(123)
freMTPL2freq <- freMTPL2freq %>% 
  sample_frac(size = 1) %>% 
  mutate(fold = rep(1:k, length = n())) %>% 
  arrange(fold)


dataset <- freMTPL2freq


# data needs to be preprocessed before fitting a poisson glm with log link
# (following paper about french case)
processed_dataset <- dataset
processed_dataset$VehGas <- factor(processed_dataset$VehGas)
processed_dataset$AreaGLM <- as.integer(processed_dataset$Area)
processed_dataset$VehPowerGLM <- as.factor(pmin(processed_dataset$VehPower,9))
VehAgeGLM <- cbind(c(0:110), c(1, rep(2,10), rep(3,100)))
processed_dataset$VehAgeGLM <- as.factor(VehAgeGLM[processed_dataset$VehAge+1,2])
processed_dataset[,"VehAgeGLM"] <-relevel(processed_dataset[,"VehAgeGLM"], ref="2")
DrivAgeGLM <- cbind(c(18:100), c(rep(1,3), rep(2,5), rep(3,5), rep(4,10), rep(5,10), rep(6,20), rep(7,30)))
processed_dataset$DrivAgeGLM <- as.factor(DrivAgeGLM[processed_dataset$DrivAge-17,2])
processed_dataset[,"DrivAgeGLM"] <-relevel(processed_dataset[,"DrivAgeGLM"], ref="5")
processed_dataset$BonusMalusGLM <- as.integer(pmin(processed_dataset$BonusMalus, 150))
processed_dataset$DensityGLM <- as.numeric(log(processed_dataset$Density))
processed_dataset[,"Region"] <-relevel(processed_dataset[,"Region"], ref="Centre")

# view(processed_dataset)

# create a vector in which store the poisson deviance per fold.

in_sample_loss_per_fold_glm <- rep(0, k)
out_sample_loss_per_fold_glm <- rep(0, k)
# write the algorithm to perform k-fold cross validation

for (i in 1:k){
  
  # use fold different than i to fit the model
  # use fold i to test the model
  
  train_dataset <- processed_dataset %>% 
    filter(fold != i)
  test_dataset <- processed_dataset %>% 
    filter(fold == i)
  
  # fit the Poisson glm model
  
  train_fit <- glm(ClaimNb ~ VehPowerGLM + VehAgeGLM + DrivAgeGLM + BonusMalusGLM
                   + VehBrand + VehGas + DensityGLM + Region + AreaGLM, 
                   data=train_dataset, offset=log(Exposure), family=poisson())
  
  # get the predictions using the test dataset
  
  train_dataset$fit <- fitted(train_fit)
  test_dataset$fit <- predict(train_fit, newdata=test_dataset, type="response")
  
  # get the Poisson deviance for both the test and train dataset (out and in sample)
  
  in_sample_loss_per_fold_glm[i] <- Poisson.Deviance(train_dataset$fit, train_dataset$ClaimNb)
  out_sample_loss_per_fold_glm[i] <- Poisson.Deviance(test_dataset$fit, test_dataset$ClaimNb)
  
}

in_sample_loss_glm <- mean(in_sample_loss_per_fold_glm)
out_sample_loss_glm <- mean(out_sample_loss_per_fold_glm)


#### GAM

library(gam)
library(mgcv)
library(mgcViz)

# data needs to be preprocessed before fitting a poisson gam with log link
processed_dataset2 <- dataset
processed_dataset2$VehGas <- factor(processed_dataset2$VehGas)
# processed_dataset2$VehPowerGLM <- as.factor(pmin(processed_dataset2$VehPower,9))
processed_dataset2$AreaGAM <- as.integer(processed_dataset2$Area)
processed_dataset2$BonusMalusGAM <- as.integer(pmin(processed_dataset2$BonusMalus, 150))
processed_dataset2$DensityGAM <- as.numeric(log(processed_dataset2$Density))
processed_dataset2[,"Region"] <-relevel(processed_dataset2[,"Region"], ref="Centre")

gam_model <- gam::gam(ClaimNb ~ s(VehPower) + s(VehAge) + s(DrivAge) + s(BonusMalus)
                      + VehBrand + VehGas + s(Density) + Region + Area, 
                      data = processed_dataset2, offset=log(Exposure), family=poisson())

# plot the fitted smooth functions
viz <- getViz(gam_model)

functions_plot <- plot(viz, allTerms = TRUE)


in_sample_loss_per_fold_gam <- rep(0, k)
out_sample_loss_per_fold_gam <- rep(0, k)

for (i in 1:k){
  
  # use fold different than i to fit the model
  # use fold i to test the model
  train_dataset <- processed_dataset2 %>% 
    filter(fold != i)
  test_dataset <- processed_dataset2 %>% 
    filter(fold == i)
  
  # fit the gam model
  train_fit <- gam::gam(ClaimNb ~ gam::s(VehPower) + gam::s(VehAge) + gam::s(DrivAge) + gam::s(BonusMalus)
                        + VehBrand + VehGas + Density + Region + Area, 
                        data = train_dataset, offset=log(Exposure), family=poisson())
  
  # get the predictions using the test dataset
  train_dataset$fit <- fitted(train_fit)
  test_dataset$fit <- gam::predict.Gam(train_fit, newdata=test_dataset, type="response")
  
  # get the Poisson deviance for both the test and train dataset (out and in sample)
  
  in_sample_loss_per_fold_gam[i] <- Poisson.Deviance(train_dataset$fit, train_dataset$ClaimNb)
  out_sample_loss_per_fold_gam[i] <- Poisson.Deviance(test_dataset$fit, test_dataset$ClaimNb)
  
}

in_sample_loss_gam <- mean(in_sample_loss_per_fold_gam)
out_sample_loss_gam <- mean(out_sample_loss_per_fold_gam)


### Regression Tree

library(rpart)
library(rpart.plot)

# we need to find the value of complexity parameter that minimizes the loss function.
cp_values_grid <- seq(from = 0,
                      to = 0.0005,
                      len = 101)

# vector for storing the error for each complexity
# parameter in cp_values_grid.
error_estimates_rt <- rep(0, times = length(cp_values_grid))

# vector for storing loss for each fold.
error_estimate_per_fold_rt <- rep(0, k)

for (j in 1:length(cp_values_grid)) {
  
  current_cp_value = cp_values_grid[j]
  
  for(i in 1:k) {
    
    train_dataset = dataset %>% 
      filter(fold != i)
    test_dataset = dataset %>% 
      filter(fold == i)
    
    # fit the regression tree
    train_tree <- rpart::rpart(cbind(Exposure , ClaimNb) ~ Area + VehPower + VehAge + DrivAge
                               + BonusMalus + VehBrand + VehGas + Density + Region ,
                               data = train_dataset , method = "poisson",
                               control = rpart.control( xval = 1,
                                                        minbucket =0.01 * nrow(train_dataset),
                                                        cp =current_cp_value))
    
    # get predictions
    test_dataset$fit <- predict(train_tree, newdata=test_dataset, type="vector") * test_dataset$Exposure
    
    # get loss
    error_estimate_per_fold[i] <- Poisson.Deviance(test_dataset$fit, test_dataset$ClaimNb)
  }
  error_estimates_rt[j] <- mean(error_estimate_per_fold_rt)
  
}

# store the cp values and the corresponding errors in a tibble, then plot them and
# get the best estimate 
table2 <- tibble(
  cp_value = cp_values_grid,
  error_estimate = error_estimates_rt
)
ggplot(table2, aes(x = cp_value, y = error_estimate)) +
  geom_point() +
  labs(x = "Complexity parameter", y = "Poisson Deviance")

best_estimate_cp <- table2 %>% slice_min(error_estimate) %>% 
  select(cp_value) %>% 
  as.double()

# using the best estimate of the complexity parameter, fit the model with cross-validation
# to obtain in and out of sample loss and compare with other models.

in_sample_loss_per_fold_rt <- rep(0, k)
out_sample_loss_per_fold_rt <- rep(0, k)

for (i in 1:k){
  
  # use fold different than i to fit the model
  # use fold i to test the model
  train_dataset <- dataset %>% 
    filter(fold != i)
  test_dataset <- dataset %>% 
    filter(fold == i)
  
  # fit the regression tree
  train_tree <- rpart::rpart(cbind(Exposure , ClaimNb) ~ Area + VehPower + VehAge + DrivAge
                             + BonusMalus + VehBrand + VehGas + Density + Region ,
                             data = train_dataset , method = "poisson",
                             control = rpart.control( xval = 1,
                                                      minbucket =0.01 * nrow(train_dataset) ,
                                                      cp =best_estimate_cp))
  
  
  # get the predictions using the test dataset
  train_dataset$fit <- predict(train_tree) * train_dataset$Exposure
  test_dataset$fit <- predict(train_tree, newdata=test_dataset, type="vector") * test_dataset$Exposure
  
  # get the Poisson deviance for both the test and train dataset (out and in sample)
  
  in_sample_loss_per_fold_rt[i] <- Poisson.Deviance(train_dataset$fit, train_dataset$ClaimNb)
  out_sample_loss_per_fold_rt[i] <- Poisson.Deviance(test_dataset$fit, test_dataset$ClaimNb)
  
}

in_sample_loss_rt <- mean(in_sample_loss_per_fold_rt)
out_sample_loss_rt <- mean(out_sample_loss_per_fold_rt)

### Random Forests
library(distRforest)
# the random forests require two tuning parameters: 1) number of trees in the forests (n_tree)
# the number of variables candidate at each split (m)
# according to hastie et al. (2008), an optimal choice for m is given by the square root of p,
# that is the number of features considered to predict the claim frequency in our case (sqrt(9))

ntree_grid <- seq(from = 50, to = 500, 25)
m <- 3
# vector for storing the error for each number of tree parameter in ntree_grid.
error_estimates_rf <- rep(0, times = length(ntree_grid))

# vector for storing loss for each fold.
error_estimate_per_fold_rf <- rep(0, k)

for (j in 1:length(ntree_grid)) {
  
  ntree_value = ntree_grid[j]
  
  for(i in 1:k) {
    
    train_dataset = dataset %>% 
      filter(fold != i)
    test_dataset = dataset %>% 
      filter(fold == i)
    
    # fit the Random Forest
    train_rf <- distRforest::rforest(
      formula = cbind(Exposure , ClaimNb) ~ Area + VehPower + VehAge + DrivAge
      + BonusMalus + VehBrand + VehGas + Density + Region,
      data = train_dataset,
      method = 'poisson',
      ntrees = ntree_value, 
      ncand = m, 
      subsample = 0.75,
      parms = list(shrink = 0.25), 
      control = rpart.control(cp = 0,
                              minbucket = 0.01 * 0.75 * nrow(train_dataset),
                              xval = 0,
                              maxcompete = 0,
                              maxsurrogate = 0),
      red_mem = TRUE 
    )
    # get predictions
    test_dataset$fit <- predict.rforest(train_rf, newdata=test_dataset) * test_dataset$Exposure
    
    # get loss
    error_estimate_per_fold_rf[i] <- Poisson.Deviance(test_dataset$fit, test_dataset$ClaimNb)
  }
  error_estimates_rf[j] <- mean(error_estimate_per_fold_rf)
  
}

# plot the errors associated to the ntree param.
table3 <- tibble(
  ntree_value = ntree_grid,
  error_estimate = error_estimates_rf
)
ggplot(table3, aes(x = ntree_value, y = error_estimate)) +
  geom_point() +
  labs(x = "number of tree", y = "Poisson Deviance")

best_estimate_ntree <- table3 %>% slice_min(error_estimate) %>% 
  select(ntree_value) %>% 
  as.double()

# using the best estimate of the number of tree parameter, fit the model with cross-validation
# to obtain in and out of sample loss and compare with other models.

in_sample_loss_per_fold_rf <- rep(0, k)
out_sample_loss_per_fold_rf <- rep(0, k)

for (i in 1:k){
  
  # use fold different than i to fit the model
  # use fold i to test the model
  train_dataset <- dataset %>% 
    filter(fold != i)
  test_dataset <- dataset %>% 
    filter(fold == i)
  
  train_rf_best <- distRforest::rforest(
    formula = cbind(Exposure , ClaimNb) ~ Area + VehPower + VehAge + DrivAge
    + BonusMalus + VehBrand + VehGas + Density + Region,
    data = train_dataset,
    method = 'poisson',
    ntrees = 143, # T in Table 3
    ncand = 7, # m in Table 3
    subsample = 0.75, # delta in Table 1
    parms = list(shrink = 0.25), # gamma in Table 1
    control = rpart.control(cp = 0, # cp in Table 1
                            minbucket = 0.01 * 0.75 * nrow(train_dataset), # kappa * delta in Table 1
                            xval = 0,
                            maxcompete = 0,
                            maxsurrogate = 0),
    red_mem = TRUE # reduces the memory footprint of individual rpart trees
  )
  
  # get the predictions using the test dataset
  train_dataset$fit <- predict.rforest(train_rf_best, newdata = train_dataset) * train_dataset$Exposure
  test_dataset$fit <- predict.rforest(train_rf_best, newdata=test_dataset) * test_dataset$Exposure  
  # get the Poisson deviance for both the test and train dataset (out and in sample)
  
  in_sample_loss_per_fold_rf[i] <- Poisson.Deviance(train_dataset$fit, train_dataset$ClaimNb)
  out_sample_loss_per_fold_rf[i] <- Poisson.Deviance(test_dataset$fit, test_dataset$ClaimNb)
  
}

in_sample_loss_rf <- mean(in_sample_loss_per_fold_rf)
out_sample_loss_rf <- mean(out_sample_loss_per_fold_rf)


# GBM params tuning

library(gbm)
ntree_grid_gbm <- seq(from = 100, to = 2000,  100)
depth <- 3

error_estimates_gbm <- rep(0, times = length(ntree_grid_gbm))


error_estimate_per_fold_gbm <- rep(0, k)

for (j in 1:length(ntree_grid_gbm)) {
  
  ntree_value = ntree_grid_gbm[j]
  
  for(i in 1:k) {
    
    train_dataset = dataset %>% 
      filter(fold != i)
    test_dataset = dataset %>% 
      filter(fold == i)
    
    # fit the gbm
    train_gbm <- gbm::gbm(
      formula = ClaimNb ~ offset(log(Exposure)) + Area + VehPower + VehAge + DrivAge
      + BonusMalus + VehBrand + as.factor(VehGas) + Density + Region,
      data = train_dataset,
      distribution = 'poisson',
      n.trees = ntree_value, 
      interaction.depth =depth, 
      shrinkage = 0.05,
      bag.fraction = 0.75, 
      n.minobsinnode = 0.01 * 0.75 * nrow(train_dataset), 
      verbose = FALSE
    )
    
    # get predictions
    test_dataset$fit <- gbm::predict.gbm(train_gbm, newdata = test_dataset, type = "response") * test_dataset$Exposure
    
    # get loss
    error_estimate_per_fold_gbm[i] <- Poisson.Deviance(test_dataset$fit, test_dataset$ClaimNb)
  }
  error_estimates_gbm[j] <- mean(error_estimate_per_fold_gbm)
  
}

# get the best estimate for ntree

table5 <- tibble(
  ntree_value = ntree_grid_gbm,
  error_estimate = error_estimates_gbm
)
ggplot(table5, aes(x = ntree_value, y = error_estimate)) +
  geom_point() +
  labs(x = "number of tree", y = "Poisson Deviance")

best_estimate_ntree <- table5 %>% slice_min(error_estimate) %>% 
  select(ntree_value) %>% 
  as.double()

# use the best estimate ntree to find best estimate depth

depth_grid <- seq(3, 9, 1)
# vector for storing the error for each complexity
# parameter in cp_values_grid.
error_estimates_gbm2 <- rep(0, times = length(depth_grid))

# vector for storing loss for each fold.
error_estimate_per_fold_gbm2 <- rep(0, k)
for (j in 1:length(depth_grid)) {
  
  depth_value = depth_grid[j]
  
  for(i in 1:k) {
    
    train_dataset = dataset %>% 
      filter(fold != i)
    test_dataset = dataset %>% 
      filter(fold == i)
    
    # a <- Sys.time()
    # fit the regression tree
    train_gbm2 <- gbm::gbm(
      formula = ClaimNb ~ offset(log(Exposure)) + Area + VehPower + VehAge + DrivAge
      + BonusMalus + VehBrand + as.factor(VehGas) + Density + Region,
      data = train_dataset,
      distribution = 'poisson',
      n.trees = 1900, # T in Table 3
      interaction.depth =depth_value, # d in Table 3
      shrinkage = 0.05, # lambda in Table 1
      bag.fraction = 0.75, # delta in Table 1
      n.minobsinnode = 0.01 * 0.75 * nrow(train_dataset), # kappa * delta in Table 1
      verbose = FALSE
    )
    
    # b <- Sys.time()
    # b-a
    # get predictions
    # train_dataset$fit <- predict.rforest(train_rf2, newdata = train_dataset) * train_dataset$Exposure
    test_dataset$fit <- gbm::predict.gbm(train_gbm2, newdata = test_dataset, type = "response") * test_dataset$Exposure
    
    # get loss
    error_estimate_per_fold_gbm2[i] <- Poisson.Deviance(test_dataset$fit, test_dataset$ClaimNb)
  }
  error_estimates_gbm2[j] <- mean(error_estimate_per_fold_gbm2)
  
}


# GBM optimal parameters already computed 

n_tree_gbm <- 1900
depth <- 3

in_sample_loss_per_fold_gbm <- rep(0, k)
out_sample_loss_per_fold_gbm <- rep(0, k)

for(i in 1:k) {
  
  train_dataset = dataset %>% 
    filter(fold != i)
  test_dataset = dataset %>% 
    filter(fold == i)
  
  # fit the gbm
  train_gbm2 <- gbm::gbm(
    formula = ClaimNb ~ offset(log(Exposure)) + Area + VehPower + VehAge + DrivAge
    + BonusMalus + VehBrand + as.factor(VehGas) + Density + Region,
    data = train_dataset,
    distribution = 'poisson',
    n.trees = n_tree_gbm, 
    interaction.depth = 3, 
    shrinkage = 0.05, 
    bag.fraction = 0.75, 
    n.minobsinnode = 0.01 * 0.75 * nrow(train_dataset),
    verbose = FALSE
  )
  
  train_dataset$fit <- exp(gbm::predict.gbm(train_gbm2, newdata = train_dataset) ) * train_dataset$Exposure
  test_dataset$fit <- gbm::predict.gbm(train_gbm2, newdata = test_dataset, type = "response") * test_dataset$Exposure
  
  # get loss
  in_sample_loss_per_fold_gbm[i] <- Poisson.Deviance(train_dataset$fit, train_dataset$ClaimNb)
  out_sample_loss_per_fold_gbm[i] <- Poisson.Deviance(test_dataset$fit, test_dataset$ClaimNb)
}


in_sample_loss_gbm <- mean(in_sample_loss_per_fold_gbm) 
out_sample_loss_gbm <- mean(out_sample_loss_per_fold_gbm)



# graph the out_sample loss for the 5 fold
a <- tibble(loss = out_sample_loss_per_fold_glm, method = "GLM")
b <- tibble(loss = out_sample_loss_per_fold_gam, method = "GAM")
c <- tibble(loss = out_sample_loss_per_fold_rt, method = "RT")
d <- tibble(loss = out_sample_loss_per_fold_rf, method = "RF")
e <- tibble(loss = out_sample_loss_per_fold_gbm, method = "GBM")
osl_to_plot <- bind_rows(a,b, c, d, e) %>% 
  group_by(method) %>% 
  mutate(fold = as.character(1:n())) %>% 
  ungroup()

library(ggplot2)


ggplot(data = osl_to_plot) +
  geom_line(aes(y = loss, x = fold, color = method, group = method)) + 
  theme_classic()


