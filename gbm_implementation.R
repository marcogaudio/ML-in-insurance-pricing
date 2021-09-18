# gbm

# parameter to be tuned in this case are: number of tree/iteration, tree depth







library(gbm)
set.seed(123)
a <- Sys.time()
gbm_freq <- gbm::gbm(
  formula = ClaimNb ~ offset(log(Exposure)) + Area + VehPower + VehAge + DrivAge
  + BonusMalus + VehBrand + as.factor(VehGas) + Density + Region,
  data = train_dataset,
  distribution = 'poisson',
  n.trees = 2500, # T in Table 3
  interaction.depth = 5, # d in Table 3
  shrinkage = 0.05, # lambda in Table 1
  bag.fraction = 0.75, # delta in Table 1
  n.minobsinnode = 0.01 * 0.75 * nrow(train_dataset), # kappa * delta in Table 1
  verbose = FALSE
)

b <- Sys.time()

b-a

# test_dataset$fit <- predict.gbm(gbm_freq, newdata = test_dataset, type = "response") * test_dataset$Exposure
test_dataset$fit <- predict(gbm_freq, newdata = test_dataset, n.trees = gbm_freq$n.trees, type = "response") * test_dataset$Exposure
Poisson.Deviance(test_dataset$fit, test_dataset$ClaimNb)

# get the Poisson deviance for both the test and train dataset (out and in sample)

in_sample_loss_per_fold_rt[i] <- Poisson.Deviance(train_dataset$fit, train_dataset$ClaimNb)
out_sample_loss_per_fold_rt[i] <- Poisson.Deviance(test_dataset$fit, test_dataset$ClaimNb)







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
a <- Sys.time()
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

b <- Sys.time()
b-a







