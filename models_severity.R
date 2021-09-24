# severity modelling 

# first thing to do is to merge the two dataset freMTPL2freq and freMTPL2sev in order to
# get only the policies with claim != 0.
library(CASdatasets)
source("utilities.R")

set.seed(123)

sev_dataset <- tibble(dplyr::left_join(freMTPL2freq, freMTPL2sev, by = c("IDpol"))) %>% 
  dplyr::filter(!is.na(ClaimAmount)) %>%  as.data.frame()

# view(sev_dataset)


# resample data and assign number of fold

k <- 5

# re-sample the data
# set a seed in order to obtain the same resample each time
sev_dataset <- sev_dataset %>% 
  sample_frac(size = 1) %>% 
  mutate(fold = rep(1:k, length = n())) %>% 
  arrange(fold)

#overview of the k fold
tot_claim_fold <- sev_dataset %>%  
  group_by(ClaimNb, fold) %>% 
  summarise(n_policy = n()) 

table_sevfold <- tot_claim_fold %>% 
  group_by(fold) %>% 
  mutate(relative_nc = n_policy/sum(n_policy))

# view(table_sevfold)


# glm implementation severity

# data needs to be preprocessed before fitting a gamma glm with log link
# (following paper about french case)
processed_dataset <- sev_dataset
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


 # create a vector in which store the poisson deviance per fold.
k <- 5 

sev_out_sample_loss_per_fold_glm <- rep(0, k)

observed_loss_per_fold <- rep(0, k)
predicted_loss_per_fold_glm <- rep(0, k)
for (i in 1:k){
  
  # use fold different than i to fit the model
  # use fold i to test the model
  train_dataset <- processed_dataset %>% 
    filter(fold != i)
  test_dataset <- processed_dataset %>% 
    filter(fold == i)
  
  # fit the  glm model
  sev_train_glm <- glm(ClaimAmount ~ VehPowerGLM + DrivAgeGLM + BonusMalusGLM
                       + VehBrand + DensityGLM + Region, weights = ClaimNb,
                   data=train_dataset, family = Gamma(link = 'log'))
  
  # get the predictions using the test dataset
  # train_dataset$fit <- fitted(sev_train_glm)
  test_dataset$fit <- predict(sev_train_glm, newdata=test_dataset, type="response") * as.vector(test_dataset$ClaimNb)
  
  # get the Poisson deviance for both the test and train dataset (out and in sample)
  
  sev_out_sample_loss_per_fold_glm[i] <- Gamma.Deviance(pred = test_dataset$fit, obs = test_dataset$ClaimAmount, alpha = as.vector(test_dataset$ClaimNb))
  
  observed_loss_per_fold[i] <- mean(train_dataset$ClaimAmount)
  predicted_loss_per_fold_glm[i] <- mean(test_dataset$fit)
}

sev_out_sample_loss_glm <- mean(sev_out_sample_loss_per_fold_glm)

observed_loss <- mean(observed_loss_per_fold)
predicted_loss_glm <- mean(predicted_loss_per_fold_glm)

# fit the model over the whole dataset to get coeffs


glm_sev_model <- glm(
  ClaimAmount ~ VehPowerGLM + DrivAgeGLM + BonusMalusGLM
  + VehBrand + DensityGLM + Region, weights = ClaimNb, 
  data=processed_dataset, family = Gamma(link = 'log')
)

summary(glm_sev_model)

## gam model for severity

library(gam)
library(mgcv)
library(mgcViz)

# data needs to be preprocessed before fitting a gamma gam with log link
processed_dataset2_sev <- sev_dataset
processed_dataset2_sev$VehGas <- factor(processed_dataset2_sev$VehGas)
# processed_dataset2_sev$VehPowerGLM <- as.factor(pmin(processed_dataset2_sev$VehPower,9))
processed_dataset2_sev$Area <- as.integer(processed_dataset2_sev$Area)
processed_dataset2_sev$BonusMalus <- as.integer(pmin(processed_dataset2_sev$BonusMalus, 150))
processed_dataset2_sev$Density <- as.numeric(log(processed_dataset2_sev$Density))
processed_dataset2_sev[,"Region"] <-relevel(processed_dataset2_sev[,"Region"], ref="Centre")

gam_model_sev <- mgcv::gam(ClaimAmount ~ s(VehPower)  + s(DrivAge) + s(BonusMalus)
                      + VehBrand  + s(Density) + Region, 
                      data = processed_dataset2_sev, weights = as.vector(ClaimNb), family=Gamma(link = 'log'))

viz_sev <- getViz(gam_model_sev)

functions_plot <- plot(viz_sev, allTerms = TRUE)




sev_out_sample_loss_per_fold_gam <- rep(0, k)


predicted_loss_per_fold_gam <- rep(0, k)

for (i in 1:k){
  
  # use fold different than i to fit the model
  # use fold i to test the model
  train_dataset <- processed_dataset2_sev %>% 
    filter(fold != i)
  test_dataset <- processed_dataset2_sev %>% 
    filter(fold == i)
  
  # fit the Poisson glm model
  sev_train_gam <- gam::gam(ClaimAmount ~ gam::s(VehPower)  + gam::s(DrivAge) + gam::s(BonusMalus)
                            + VehBrand + gam::s(Density) + Region, 
                            data = train_dataset, weights = as.vector(ClaimNb), family=Gamma(link = 'log'))
  
  # get the predictions using the test dataset
  # train_dataset$fit <- fitted(sev_train_gam)
  test_dataset$fit <- gam::predict.Gam(sev_train_gam, newdata=test_dataset, type="response") * as.vector(test_dataset$ClaimNb)
  
  # get the Poisson deviance for both the test and train dataset (out and in sample)
  
  sev_out_sample_loss_per_fold_gam[i] <- Gamma.Deviance(pred = test_dataset$fit, obs = test_dataset$ClaimAmount, alpha = as.vector(test_dataset$ClaimNb))
  
  
  predicted_loss_per_fold_gam[i] <- mean(test_dataset$fit)
}

sev_out_sample_loss_gam <- mean(sev_out_sample_loss_per_fold_gam)

predicted_loss_gam <- mean(predicted_loss_per_fold_gam)


# regression tree severity

library(rpart)
library(rpart.plot)

# we need to find the value of complexity parameter that minimizes the loss function.
sev_cp_values_grid <- seq(from = 0,
                      to = 0.007,
                      len = 500)

# vector for storing the error for each complexity
# parameter in cp_values_grid.
sev_error_estimates_rt <- rep(0, times = length(sev_cp_values_grid))

# vector for storing loss for each fold.
sev_error_estimate_per_fold_rt <- rep(0, k)

for (j in 1:length(sev_cp_values_grid)) {
  
  sev_current_cp_value = sev_cp_values_grid[j]
  
  for(i in 1:k) {
    
    train_dataset = sev_dataset %>% 
      filter(fold != i)
    test_dataset = sev_dataset %>% 
      filter(fold == i)
    
    # fit the regression tree
    sev_train_tree <- distRforest::rpart(ClaimAmount ~ VehPower + DrivAge 
                                         + BonusMalus + VehBrand  + Density + Region,
                               data = train_dataset , method = "gamma", weights = ClaimNb, 
                               control = rpart.control(xval = 0,
                                                      maxcompete = 0, maxsurrogate = 0, 
                                                      minbucket = 0.02 * nrow(train_dataset),
                                                      cp = sev_current_cp_value))
    
    # get predictions
    test_dataset$fit <- predict(sev_train_tree, newdata=test_dataset, type="vector") * as.vector(test_dataset$ClaimNb)
    
    # get loss
    sev_error_estimate_per_fold_rt[i] <- Gamma.Deviance(pred = test_dataset$fit, obs = test_dataset$ClaimAmount, alpha = as.vector(test_dataset$ClaimNb))
  }
  sev_error_estimates_rt[j] <- mean(sev_error_estimate_per_fold_rt)
  
}

# collect the cp values and the related errors, plot them and get the best estimate
table2_sev <- tibble(
  cp_value = sev_cp_values_grid,
  error_estimate = sev_error_estimates_rt
)
ggplot(table2_sev, aes(x = cp_value, y = error_estimate)) +
  geom_point() +
  labs(x = "Complexity parameter", y = "Poisson Deviance")

best_estimate_cp_sev <- table2_sev %>% slice_min(error_estimate) %>% 
  select(cp_value) %>%
  pull() 

best_estimate_cp_sev <- as.double(best_estimate_cp_sev[1])

# using the best estimate of the complexity parameter, fit the model with cross-validation
# to obtain in and out of sample loss and compare with other models.

sev_out_sample_loss_per_fold_rt <- rep(0, k)

predicted_loss_per_fold_rt <- rep(0, k)
for (i in 1:k){
  
  # use fold different than i to fit the model
  # use fold i to test the model
  train_dataset <- sev_dataset %>% 
    filter(fold != i)
  test_dataset <- sev_dataset %>% 
    filter(fold == i)
  
  # fit the regression tree
  sev_train_tree <- distRforest::rpart(ClaimAmount ~ VehPower  + DrivAge 
                                       + BonusMalus + VehBrand + Density + Region,
                             data = train_dataset , method = "gamma",
                             control = rpart.control( xval = 0,
                                                      minbucket = 0.02 * nrow(train_dataset),
                                                      
                                                      maxcompete = 0,
                                                      maxsurrogate = 0,
                                                      cp =best_estimate_cp_sev))
  
  
  # get the predictions using the test dataset
  train_dataset$fit <- predict(sev_train_tree) * as.vector(train_dataset$ClaimNb)
  test_dataset$fit <-  predict(sev_train_tree, newdata=test_dataset, type="vector") * as.vector(test_dataset$ClaimNb)
  
  # get the Poisson deviance for both the test and train dataset (out and in sample)
  
  sev_out_sample_loss_per_fold_rt[i] <- Gamma.Deviance(pred = test_dataset$fit, obs = test_dataset$ClaimAmount, alpha = as.vector(test_dataset$ClaimNb))
  
  
  
  predicted_loss_per_fold_rt[i] <- mean(test_dataset$fit)
}

sev_out_sample_loss_rt <- mean(sev_out_sample_loss_per_fold_rt)

predicted_loss_rt <- mean(predicted_loss_per_fold_rt)



# random forest
library(distRforest)

ntree_grid <- seq(from = 50, to = 300, 10)
m <- sqrt(6)
# vector for storing the error for each number of tree parameter in ntree_grid.
error_estimates_rf <- rep(0, times = length(ntree_grid))

# vector for storing loss for each fold.
error_estimate_per_fold_rf <- rep(0, k)

for (j in 1:length(ntree_grid)) {
  
  ntree_value = ntree_grid[j]
  
  for(i in 1:k) {
    
    train_dataset = sev_dataset %>% 
      filter(fold != i)
    test_dataset = sev_dataset %>% 
      filter(fold == i)
    
    # fit the random forest
    train_rf <- distRforest::rforest(
      formula = ClaimAmount ~ VehPower + DrivAge
      + BonusMalus + VehBrand + Density + Region,
      data = train_dataset,
      method = 'gamma',
      ntrees = ntree_value, 
      ncand = m, 
      subsample = 0.75,
      parms = list(shrink = 0.25), 
      control = rpart.control(cp = 0,
                              minbucket = 0.02 * 0.75 * nrow(train_dataset),
                              xval = 0,
                              maxcompete = 0,
                              maxsurrogate = 0),
      red_mem = TRUE,
      keep_data = TRUE
    )
    # get predictions
    test_dataset$fit <- predict.rforest(train_rf, newdata=test_dataset) * as.vector(test_dataset$ClaimNb) 
    
    # get loss
    error_estimate_per_fold_rf[i] <- Gamma.Deviance(pred = test_dataset$fit,obs = test_dataset$ClaimAmount, alpha = as.vector(test_dataset$ClaimNb))
  }
  error_estimates_rf[j] <- mean(error_estimate_per_fold_rf)
  
}

# get the best estimate of ntree param.
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

sev_out_sample_loss_per_fold_rf <- rep(0, k)

predicted_loss_per_fold_rf <- rep(0, k)

for (i in 1:k){
  
  # use fold different than i to fit the model
  # use fold i to test the model
  train_dataset <- sev_dataset %>% 
    filter(fold != i)
  test_dataset <- sev_dataset %>% 
    filter(fold == i)
  
  train_rf_best <- distRforest::rforest(
    formula = ClaimAmount ~ VehPower + DrivAge
    + BonusMalus + VehBrand  + Density + Region,
    data = train_dataset,
    method = 'gamma',
    ntrees = best_estimate_ntree, 
    weights = ClaimNb,
    ncand = 3, 
    subsample = 0.75, 
    parms = list(shrink = 0.25), 
    control = rpart.control(cp = 0,
                            minbucket = 0.02 * 0.75 * nrow(train_dataset),
                            xval = 0,
                            maxcompete = 0,
                            maxsurrogate = 0),
    red_mem = TRUE, 
    keep_data = TRUE
  )
  
  # get the predictions using the test dataset
  test_dataset$fit <- predict.rforest(train_rf_best, newdata=test_dataset) * as.vector(test_dataset$ClaimNb)
  # get the Poisson deviance for both the test and train dataset (out and in sample)
  
  sev_out_sample_loss_per_fold_rf[i] <- Gamma.Deviance(pred = test_dataset$fit, obs = test_dataset$ClaimAmount, alpha = as.vector(test_dataset$ClaimNb))
  
  predicted_loss_per_fold_rf[i] <- mean(test_dataset$fit)
}

sev_out_sample_loss_rf <- mean(sev_out_sample_loss_per_fold_rf)

predicted_loss_rf <- mean(predicted_loss_per_fold_rf)


### gbm for claim severity 
# devtools::install_github('harrysouthworth/gbm')
library(gbm)
# function for predicting, due to the unavailability in gbm package from github
predict_model.gbm <- function(object, newdata) {
  predict(object, newdata, n.trees = object$n.trees, type = 'response')
}

sev_ntree_grid_gbm <- seq(from = 25, to = 70,  1)
depth <- 3

sev_error_estimates_gbm <- rep(0, times = length(sev_ntree_grid_gbm))


sev_error_estimate_per_fold_gbm <- rep(0, k)

for (j in 1:length(sev_ntree_grid_gbm)) {
  
  ntree_value = sev_ntree_grid_gbm[j]
  
  for(i in 1:k) {
    
    train_dataset = sev_dataset %>% 
      filter(fold != i)
    test_dataset = sev_dataset %>% 
      filter(fold == i)
    
    # fit the gbm
    train_gbm <- gbm::gbm(
      formula = ClaimAmount ~  VehPower + DrivAge
      + BonusMalus + VehBrand + Density + Region,
      data = train_dataset,
      distribution = 'gamma',
      weights = ClaimNb,
      n.trees = ntree_value, 
      interaction.depth =depth, 
      shrinkage = 0.05,
      bag.fraction = 0.75, 
      n.minobsinnode = 0.02 * 0.75 * nrow(train_dataset), 
      verbose = FALSE
    )
    
    # get predictions
    test_dataset$fit <-  predict_model.gbm(train_gbm, newdata = test_dataset)  * as.vector(test_dataset$ClaimNb)
    
    sev_error_estimate_per_fold_gbm[i] <- Gamma.Deviance(pred = test_dataset$fit, obs = test_dataset$ClaimAmount, alpha = as.vector(test_dataset$ClaimNb))
  }
  sev_error_estimates_gbm[j] <- mean(sev_error_estimate_per_fold_gbm)
  
}

# get the best estimate for ntree

sev_table5 <- tibble(
  ntree_value = sev_ntree_grid_gbm,
  error_estimate = sev_error_estimates_gbm
)
ggplot(sev_table5, aes(x = ntree_value, y = error_estimate)) +
  geom_point() +
  labs(x = "number of tree", y = "Gamma Deviance")

sev_best_estimate_ntree <- sev_table5 %>% slice_min(error_estimate) %>% 
  select(ntree_value) %>% 
  as.double()


# tune depth 

depth_grid <- seq(1, 9, 1)
# vector for storing the error for each complexity
# parameter in cp_values_grid.
sev_error_estimates_gbm2 <- rep(0, times = length(depth_grid))

# vector for storing loss for each fold.
sev_error_estimate_per_fold_gbm2 <- rep(0, k)
a <- Sys.time()
for (j in 1:length(depth_grid)) {
  
  depth_value = depth_grid[j]
  
  for(i in 1:k) {
    
    train_dataset = sev_dataset %>% 
      filter(fold != i)
    test_dataset = sev_dataset %>% 
      filter(fold == i)
    
    # a <- Sys.time()
    # fit the regression tree
    train_gbm2 <- gbm::gbm(
      formula = ClaimAmount ~  VehPower + DrivAge
      + BonusMalus + VehBrand + Density + Region,
      data = train_dataset,
      distribution = 'gamma',
      weights = ClaimNb,
      n.trees = sev_best_estimate_ntree, # T in Table 3
      interaction.depth =depth_value, # d in Table 3
      shrinkage = 0.05, # lambda in Table 1
      bag.fraction = 0.75, # delta in Table 1
      n.minobsinnode = 0.0 * 0.75 * nrow(train_dataset), # kappa * delta in Table 1
      verbose = FALSE
    )
    
    
    # get predictions
    test_dataset$fit <- predict_model.gbm(train_gbm2, newdata = test_dataset) * as.vector(test_dataset$ClaimNb)
    
    # get loss
    sev_error_estimate_per_fold_gbm2[i] <- Gamma.Deviance(pred = test_dataset$fit, obs = test_dataset$ClaimAmount, alpha = as.vector(test_dataset$ClaimNb))
  }
  sev_error_estimates_gbm2[j] <- mean(sev_error_estimate_per_fold_gbm2)
  
}

### 2 is the best for depth


# k-fold with best parameters



# GBM optimal parameters already computed 

sev_n_tree_gbm <- sev_best_estimate_ntree
sev_depth <- 2

sev_out_sample_loss_per_fold_gbm <- rep(0, k)

predicted_loss_per_fold_gbm <- rep(0, k)


for(i in 1:k) {
  
  train_dataset = sev_dataset %>% 
    filter(fold != i)
  test_dataset = sev_dataset %>% 
    filter(fold == i)
  
  # fit the gbm
  sev_train_gbm <- gbm::gbm(
    formula = ClaimAmount ~  VehPower +  DrivAge
    + BonusMalus + VehBrand + Density + Region,
    data = train_dataset,
    distribution = 'gamma',
    n.trees = sev_n_tree_gbm, 
    interaction.depth = sev_depth, 
    shrinkage = 0.05,
    bag.fraction = 0.75,
    n.minobsinnode = 0.02 * 0.75 * nrow(train_dataset),
    verbose = FALSE
  )
  
  test_dataset$fit <- predict_model.gbm(sev_train_gbm, newdata = test_dataset)  * as.vector(test_dataset$ClaimNb)
  
  # get loss
  sev_out_sample_loss_per_fold_gbm[i] <- Gamma.Deviance(pred = test_dataset$fit, obs = test_dataset$ClaimAmount, alpha = as.vector(test_dataset$ClaimNb))
  
  predicted_loss_per_fold_gbm[i] <- mean(test_dataset$fit)
  
}

sev_out_sample_loss_gbm <- mean(sev_out_sample_loss_per_fold_gbm)

predicted_loss_gbm <- mean(predicted_loss_per_fold_gbm)



# table with out-of-sample losses
severity_table <- tibble(
  glm = sev_out_sample_loss_glm,
  gam = sev_out_sample_loss_gam ,
  rt = sev_out_sample_loss_rt,
  rf = sev_out_sample_loss_rf,
  gbm = sev_out_sample_loss_gbm
) %>% view()






# graph the out_sample loss for the 5 fold

as <- tibble(loss = sev_out_sample_loss_per_fold_glm, method = "GLM")
bs <- tibble(loss = sev_out_sample_loss_per_fold_gam, method = "GAM")
cs <- tibble(loss = sev_out_sample_loss_per_fold_rt, method = "RT")
ds <- tibble(loss = sev_out_sample_loss_per_fold_rf, method = "RF")
es <- tibble(loss = sev_out_sample_loss_per_fold_gbm, method = "GBM")
oslsev_to_plot <- bind_rows(as, bs, cs, ds, es) %>% 
  group_by(method) %>% 
  mutate(fold = as.character(1:n())) %>% 
  ungroup()

library(ggplot2)


ggplot(data = oslsev_to_plot) +
  geom_line(aes(y = loss, x = fold, color = method, group = method)) + 
  theme_classic()


