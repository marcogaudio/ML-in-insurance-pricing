# regression tree
library(rpart)
library(rpart.plot)
k <- 5


# data already folded



# we need to find the value of complexity parameter that minimizes the loss function.
cp_values_grid <- seq(from = 0,
                      to = 0.0005,
                      len = 101)

# vector for storing the error for each complexity
# parameter in cp_values_grid.
error_estimates <- rep(0, times = length(cp_values_grid))

# vector for storing loss for each fold.
error_estimate_per_fold <- rep(0, k)

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
  error_estimates[j] <- mean(error_estimate_per_fold)
  
}


table2 <- tibble(
  cp_value = cp_values_grid,
  error_estimate = error_estimates
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
  
  # fit the Poisson glm model
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



# to have an overview of the model, a regression tree is fitted over the whole dataset.

tree_fold5 <- rpart::rpart(cbind(Exposure , ClaimNb) ~ Area + VehPower + VehAge + DrivAge
                           + BonusMalus + VehBrand + VehGas + Density + Region ,
                           data = train_dataset , method = "poisson",
                           control = rpart.control( xval = 1,
                                                    minbucket =0.01 * nrow(train_dataset),
                                                    cp =best_estimate_cp))
tree_model

rpart.plot(tree_model)
