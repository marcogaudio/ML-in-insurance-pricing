# models

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


view(processed_dataset)



# create a vector in which store the poisson deviance per fold.

in_sample_loss_per_fold_glm <- rep(0, k)
out_sample_loss_per_fold_glm <- rep(0,k)
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
  
  # y_hat <- predict(train_fit, newdata = test_dataset, type ="response")
  
  
  train_dataset$fit <- fitted(train_fit)
  test_dataset$fit <- predict(train_fit, newdata=test_dataset, type="response")
  
  
  # c(Poisson.Deviance(train_dataset$fit, train_dataset$ClaimNb), Poisson.Deviance(test_dataset$fit, test_dataset$ClaimNb))
  
  # get the Poisson deviance for both the test and train dataset (out and in sample)
  
  in_sample_loss_per_fold_glm[i] <- Poisson.Deviance(train_dataset$fit, train_dataset$ClaimNb)
  out_sample_loss_per_fold_glm[i] <- Poisson.Deviance(test_dataset$fit, test_dataset$ClaimNb)
  
}

in_sample_loss_glm <- mean(in_sample_loss_per_fold)
out_sample_loss_glm <- mean(out_sample_loss_per_fold)


#### GAM

library(gam)
library(mgcv)
library(mgcViz)

# data needs to be preprocessed before fitting a poisson glm with log link
# (following paper about french case)
processed_dataset2 <- dataset
processed_dataset2$VehGas <- factor(processed_dataset2$VehGas)
processed_dataset2$VehPowerGLM <- as.factor(pmin(processed_dataset2$VehPower,9))
processed_dataset2$AreaGAM <- as.integer(processed_dataset2$Area)
processed_dataset2$BonusMalusGAM <- as.integer(pmin(processed_dataset2$BonusMalus, 150))
processed_dataset2$DensityGAM <- as.numeric(log(processed_dataset2$Density))
processed_dataset2[,"Region"] <-relevel(processed_dataset2[,"Region"], ref="Centre")

a <- Sys.time()
gam_model <- gam::gam(ClaimNb ~ gam::s(VehPower) + gam::s(VehAge) + gam::s(DrivAge) + gam::s(BonusMalus)
                        + VehBrand + VehGas + Density + Region + Area, 
                       data = processed_dataset2, offset=log(Exposure), family=poisson())
b <- Sys.time()
b-a
#plot the fitted smooth functions
viz <- getViz(gam_model)

functions_plot <- plot(viz, allTerms = TRUE) +
  l_points() +
  l_fitLine(linetype = 1)  +
  l_ciLine(linetype = 3) +
  l_ciBar() +
  l_rug() +
  theme_grey() 

print(functions_plot, pages = 1)


in_sample_loss_per_fold_gam <- rep(0, k)
out_sample_loss_per_fold_gam <- rep(0, k)

for (i in 1:k){
  
  # use fold different than i to fit the model
  # use fold i to test the model
  train_dataset <- dataset %>% 
    filter(fold != i)
  test_dataset <- dataset %>% 
    filter(fold == i)
  
  # fit the Poisson glm model
  train_fit <- gam::gam(ClaimNb ~ VehPower + s(VehAge) + s(DrivAge) + s(BonusMalus)
                        + VehBrand + VehGas + s(Density) + Region + Area, 
                        data = train_dataset, offset=log(Exposure), family=poisson())
  
  # get the predictions using the test dataset
  train_dataset$fit <- fitted(train_fit)
  test_dataset$fit <- predict(train_fit, newdata=test_dataset, type="response")
  
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

### Random Forests







