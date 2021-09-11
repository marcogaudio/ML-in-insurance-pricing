# GLM implementation
set.seed(123)

#FREQUENCY

# k-fold cross validation

# choose a number of fold (optimal number usually 5 or 10)
k <- 5

# re-sample the data
# set a seed in order to obtain the same resample each time
dataset <- freMTPL2freq %>% 
  sample_frac(size = 1) %>% 
  mutate(fold = rep(1:k, length = n())) %>% 
  arrange(fold)


# overview of the k folds

tot_claim_fold <- dataset %>%  
  group_by(ClaimNb, fold) %>% 
  summarise(n_policy = n()) 

table3.1.3 <- tot_claim_fold %>% 
  group_by(fold) %>% 
  mutate(relative_nc = n_policy/sum(n_policy))

view(table3.1.3)

tot_claim_fold %>% filter(fold == 1) %>% pull(n_policy) %>%  sum()


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
k <- 5 

in_sample_loss_per_fold <- rep(0, k)
out_sample_loss_per_fold <- rep(0, k)

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
  
  in_sample_loss_per_fold[i] <- Poisson.Deviance(train_dataset$fit, train_dataset$ClaimNb)
  out_sample_loss_per_fold[i] <- Poisson.Deviance(test_dataset$fit, test_dataset$ClaimNb)
  
  }

in_sample_loss <- mean(in_sample_loss_per_fold)
out_sample_loss <- mean(out_sample_loss_per_fold)


# to get the regression coefficients, we fit a glm over the whole dataset

glm_model <- glm(ClaimNb ~ VehPowerGLM + VehAgeGLM + DrivAgeGLM + BonusMalusGLM
                 + VehBrand + VehGas + DensityGLM + Region + AreaGLM, 
                 data = processed_dataset, offset=log(Exposure), family=poisson())

summary(glm_model)



# write a function for the Poisson deviance

Poisson.Deviance <- function(pred, obs){
  
  200*(sum(pred)-sum(obs)+sum(log((obs/pred)^(obs))))/length(pred)
  
  }



