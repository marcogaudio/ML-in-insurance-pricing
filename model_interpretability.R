# model explainability

library(tidyverse)
library(MetricsWeighted)
library(flashlight)
library(ggplot2)

# in order to compare the interpretability of the models, these are fitted over
# the whole dataset

# we need the data 

source("utilities.R")

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
processed_dataset$Area <- as.integer(processed_dataset$Area)
processed_dataset$VehPower <- as.factor(pmin(processed_dataset$VehPower,9))
VehAge <- cbind(c(0:110), c(1, rep(2,10), rep(3,100)))
processed_dataset$VehAge <- as.factor(VehAge[processed_dataset$VehAge+1,2])
processed_dataset[,"VehAge"] <-relevel(processed_dataset[,"VehAge"], ref="1")
DrivAge <- cbind(c(18:100), c(rep(1,3), rep(2,5), rep(3,5), rep(4,10), rep(5,10), rep(6,20), rep(7,30)))
processed_dataset$DrivAge <- as.factor(DrivAge[processed_dataset$DrivAge-17,2])
processed_dataset[,"DrivAge"] <-relevel(processed_dataset[,"DrivAge"], ref="1")
processed_dataset$BonusMalus <- as.integer(pmin(processed_dataset$BonusMalus, 150))
processed_dataset$Density <- as.numeric(log(processed_dataset$Density))
processed_dataset[,"Region"] <-relevel(processed_dataset[,"Region"], ref="Centre")

#FREQUENCY

# glm model

glm_model <- glm(ClaimNb ~ VehPower  + DrivAge + BonusMalus
                 + VehBrand  + Density + Region , 
                 data = processed_dataset, offset=log(Exposure), family=poisson())



fl_glm <- flashlight(
   model = glm_model , label = "GLM" ,
   predict_function = function(fit , X) predict(fit , X , type = "response")
  )

# GAM
library(gam)
library(mgcv)


# data needs to be preprocessed before fitting a poisson gam with log link
# (following paper about french case)
processed_dataset2 <- dataset
processed_dataset2$VehGas <- factor(processed_dataset2$VehGas)
# processed_dataset2$VehPowerGLM <- as.factor(pmin(processed_dataset2$VehPower,9))
processed_dataset2$Area <- as.integer(processed_dataset2$Area)
processed_dataset2$BonusMalus <- as.integer(pmin(processed_dataset2$BonusMalus, 150))
processed_dataset2$Density <- as.numeric(log(processed_dataset2$Density))
processed_dataset2[,"Region"] <-relevel(processed_dataset2[,"Region"], ref="Centre")

gam_model <- gam::gam(ClaimNb ~ gam::s(VehPower) + gam::s(VehAge) + gam::s(DrivAge) + gam::s(BonusMalus)
                      + VehBrand + VehGas + gam::s(Density) + Region + Area, 
                      data = processed_dataset2, offset=log(Exposure), family=poisson())

fl_gam <- flashlight (
  model = gam_model , label = "GAM" ,
  predict_function = function(fit , X) gam::predict.Gam(fit, X, "response")
)


# Regression Tree

RT_model <- rpart::rpart(cbind(Exposure , ClaimNb) ~ Area + VehPower + VehAge + DrivAge
                         + BonusMalus + VehBrand + VehGas + Density + Region ,
                         data = dataset , method = "poisson",
                         control = rpart.control( xval = 1,
                                                  minbucket =0.01 * nrow(dataset),
                                                  cp =0.000035))


fl_rt <- flashlight (
  model = RT_model , label = "RT" ,
  predict_function = function(fit , X) predict(fit, X, "vector")
)

# random forrest

RF_model <- distRforest::rforest(
  formula = cbind(Exposure , ClaimNb) ~ Area + VehPower + VehAge + DrivAge
  + BonusMalus + VehBrand + VehGas + Density + Region,
  data = dataset,
  method = 'poisson',
  ntrees = 143, 
  ncand = 7, 
  subsample = 0.75,
  parms = list(shrink = 0.25), 
  control = rpart.control(cp = 0,
                          minbucket = 0.01 * 0.75 * nrow(dataset),
                          xval = 0,
                          maxcompete = 0,
                          maxsurrogate = 0),
  red_mem = TRUE 
)

fl_rf <- flashlight(
  model = RF_model, label = "RF", 
  predict_function = function(fit, X) predict.rforest(fit, X)
)


# GBM

GBM_model <- gbm::gbm(
  formula = ClaimNb ~ offset(log(Exposure)) + Area + VehPower + VehAge + DrivAge
  + BonusMalus + VehBrand + as.factor(VehGas) + Density + Region,
  data = dataset,
  distribution = 'poisson',
  n.trees = 1900, 
  interaction.depth = 3, 
  shrinkage = 0.05, 
  bag.fraction = 0.75, 
  n.minobsinnode = 0.01 * 0.75 * nrow(dataset),
  verbose = FALSE
)


fl_gbm <- flashlight(
  model = GBM_model, label = "GBM", 
  predict_function = function(fit, X) predict(fit, X, n.trees = 1900, "response")
)


metrics <- list("Average deviance" = deviance_poisson ,
                 "Relative deviance reduction" = r_squared_poisson)

fls_glm <- multiflashlight (list(fl_glm), data = processed_dataset ,
                            y = "ClaimNb", w = "Exposure", metrics = metrics)
fls_gam <- multiflashlight (list(fl_gam), data = processed_dataset2 ,
                            y = "ClaimNb", w = "Exposure", metrics = metrics)
fls_rt <- multiflashlight(list(fl_rt), data = dataset,
                          y = "ClaimNb", w = "Exposure", metrics = metrics)
fls_rf <- multiflashlight(list(fl_rf), data = dataset,
                          y = "ClaimNb", w = "Exposure", metrics = metrics)
fls_gbm <- multiflashlight(list(fl_gbm), data = dataset,
                           y = "ClaimNb", w = "Exposure", metrics = metrics)
fls_ml <- multiflashlight(list(fl_rt, fl_rf, fl_gbm), data = dataset,
                          y = "ClaimNb",w = "Exposure", metrics = metrics)

# VARIABLE IMPORTANCE FREQUENCY
x <- c("VehPower", "VehAge", "DrivAge", "BonusMalus",
       "VehBrand", "VehGas", "Density", "Region", "Area")

imp_glm <- light_importance(fls_glm, v = x)
plot(imp_glm, fill = "#0072B2", color = "black")

imp_gam <- light_importance(fls_gam, v = x)
plot(imp_gam, fill = "#0072B2", color = "black")

imp_rt <- light_importance(fls_rt, v = x)
plot(imp_rt, fill = "#0072B2", color = "black")

imp_rf <- light_importance(fls_rf, v = x)
plot(imp_rf, fill = "#0072B2", color = "black")

imp_gbm <- light_importance(fls_gbm, v = x)
plot(imp_gbm, fill = "#0072B2", color = "black")


# SEVERITY 
set.seed(123)


library(CASdatasets)
source("utilities.R")

sev_dataset <- tibble(dplyr::left_join(freMTPL2freq, freMTPL2sev, by = c("IDpol"))) %>% 
  dplyr::filter(!is.na(ClaimAmount)) %>%  as.data.frame()

# resample data and assign number of fold

k <- 5

# re-sample the data
# set a seed in order to obtain the same resample each time
sev_dataset <- sev_dataset %>% 
  sample_frac(size = 1) %>% 
  mutate(fold = rep(1:k, length = n())) %>% 
  arrange(fold)


# GLM severity 
processed_dataset <- sev_dataset
processed_dataset$VehGas <- factor(processed_dataset$VehGas)
processed_dataset$Area <- as.integer(processed_dataset$Area)
processed_dataset$VehPower <- as.factor(pmin(processed_dataset$VehPower,9))
VehAge <- cbind(c(0:110), c(1, rep(2,10), rep(3,100)))
processed_dataset$VehAge <- as.factor(VehAge[processed_dataset$VehAge+1,2])
processed_dataset[,"VehAge"] <-relevel(processed_dataset[,"VehAge"], ref="1")
DrivAge <- cbind(c(18:100), c(rep(1,3), rep(2,5), rep(3,5), rep(4,10), rep(5,10), rep(6,20), rep(7,30)))
processed_dataset$DrivAge <- as.factor(DrivAge[processed_dataset$DrivAge-17,2])
processed_dataset[,"DrivAge"] <-relevel(processed_dataset[,"DrivAge"], ref="1")
processed_dataset$BonusMalus <- as.integer(pmin(processed_dataset$BonusMalus, 150))
processed_dataset$Density <- as.numeric(log(processed_dataset$Density))
processed_dataset[,"Region"] <-relevel(processed_dataset[,"Region"], ref="Centre")

glm_sev_model <- glm(
  ClaimAmount ~ VehPower + DrivAge + BonusMalus
  + VehBrand + Density + Region, weights = ClaimNb, 
  data=processed_dataset, family = Gamma(link = 'log')
)

fl_glm_sev <- flashlight(
  model = glm_sev_model , label = "GLM",
  predict_function = function(fit , X) predict(fit , X , type = "response")
)

# GAM severity 

processed_dataset2_sev <- sev_dataset
processed_dataset2_sev$VehGas <- factor(processed_dataset2_sev$VehGas)
# processed_dataset2_sev$VehPowerGLM <- as.factor(pmin(processed_dataset2_sev$VehPower,9))
processed_dataset2_sev$Area <- as.integer(processed_dataset2_sev$Area)
processed_dataset2_sev$BonusMalus <- as.integer(pmin(processed_dataset2_sev$BonusMalus, 150))
processed_dataset2_sev$Density <- as.numeric(log(processed_dataset2_sev$Density))
processed_dataset2_sev[,"Region"] <-relevel(processed_dataset2_sev[,"Region"], ref="Centre")


gam_sev_model <- mgcv::gam(ClaimAmount ~ s(VehPower)  + s(DrivAge) + s(BonusMalus)
                           + VehBrand  + s(Density) + Region, 
                           data = processed_dataset2_sev, weights = as.vector(ClaimNb), family=Gamma(link = 'log'))
fl_gam_sev <- flashlight(
  model = gam_sev_model , label = "GAM",
  predict_function = function(fit , X) predict(fit , X , type = "response")
)

# REGRESSION TREE

sev_rt <- distRforest::rpart(ClaimAmount ~ VehPower  + DrivAge 
                                     + BonusMalus + VehBrand + Density + Region,
                                     data = sev_dataset, method = "gamma", weights = ClaimNb,
                                     control = rpart.control( xval = 0,
                                                              minbucket = 0.02 * nrow(sev_dataset),
                                                              maxcompete = 0,
                                                              maxsurrogate = 0,
                                                              cp =0.005639279))

fl_rt_sev <- flashlight(
  model = sev_rt, label = "RT",
  predict_function = function(fit, X) predict(fit, X, "vector")
)

# RANDOM FOREST
 
sev_rf <- distRforest::rforest(
  formula = ClaimAmount ~ VehPower + DrivAge
  + BonusMalus + VehBrand  + Density + Region,
  data = sev_dataset,
  method = 'gamma',
  ntrees = 160, 
  weights = ClaimNb,
  ncand = sqrt(6), 
  subsample = 0.75, 
  parms = list(shrink = 0.25), 
  control = rpart.control(cp = 0,
                          minbucket = 0.02 * 0.75 * nrow(sev_dataset),
                          xval = 0,
                          maxcompete = 0,
                          maxsurrogate = 0),
  red_mem = TRUE, 
  keep_data = TRUE
)


fl_rf_sev <- flashlight(
  model = sev_rf, label = "RF",
  predict_function = function(fit, X) predict.rforest(fit, X)
  )

## GBM severity

sev_gbm <- gbm::gbm(
  formula = ClaimAmount ~  VehPower +  DrivAge
  + BonusMalus + VehBrand + Density + Region,
  weights = ClaimNb,
  data = sev_dataset,
  distribution = 'gamma',
  n.trees = 56, 
  interaction.depth = 2, 
  shrinkage = 0.05,
  bag.fraction = 0.75,
  n.minobsinnode = 0.02 * 0.75 * nrow(sev_dataset),
  verbose = FALSE
)


fl_gbm_sev <- flashlight(
  model = sev_gbm, label = "GBM",
  predict_function = function(fit, X) predict(fit, X, n.trees = 56, type = 'response')
    
)

metrics_sev <- list("Average deviance" = deviance_gamma,
                "Relative deviance reduction" = r_squared_poisson )

fls_glm_sev <- multiflashlight (list(fl_glm_sev), data = processed_dataset ,
                            y = "ClaimAmount", w = "ClaimNb", metrics = metrics_sev)
fls_gam_sev <- multiflashlight(list(fl_gam_sev), data = processed_dataset2_sev, 
                               y = "ClaimAmount", w = "ClaimNb", metrics = metrics_sev)
fls_rt_sev <- multiflashlight(list(fl_rt_sev), data = sev_dataset,
                              y = "ClaimAmount", w = "ClaimNb", metrics = metrics_sev)
fls_rf_sev <- multiflashlight(list(fl_rf_sev), data = sev_dataset,
                              y = "ClaimAmount", w = "ClaimNb", metrics = metrics_sev)
fls_gbm_sev <- multiflashlight(list(fl_gbm_sev), data = sev_dataset,
                              y = "ClaimAmount", w = "ClaimNb", metrics = metrics_sev)
fls_ml_sev <-  multiflashlight(list(fl_rt_sev, fl_rf_sev, fl_gbm_sev), data = sev_dataset,
                               y = "ClaimAmount", w = "ClaimNb", metrics = metrics_sev)


x_sev <- c("VehPower", "DrivAge", "BonusMalus",
       "VehBrand", "Density", "Region")
# VARIABLE IMPORTANCE SEVERITY

imp_glm_sev <- light_importance(fls_glm_sev, v = x_sev)
imp_gam_sev <- light_importance(fls_gam_sev, v = x_sev)
imp_rt_sev <- light_importance(fls_rt_sev, v = x_sev)
imp_rf_sev <- light_importance(fls_rf_sev, v = x_sev)
imp_gbm_sev <- light_importance(fls_gbm_sev, v = x_sev)

plot(imp_glm_sev, fill = "#0072B2", color = "black") 
plot(imp_gam_sev, fill = "#0072B2", color = "black")
plot(imp_rt_sev, fill = "#0072B2", color = "black")
plot(imp_rf_sev, fill = "#0072B2", color = "black")
plot(imp_gbm_sev, fill = "#0072B2", color = "black")


### INDIVIDUAL CONDITIONAL EXPECTATION FREQ
plot(light_ice(fls_glm, v = "DrivAge" , n_max = 200 , seed = 3) , alpha = 0.1) +
  ylab("frequency") 
plot(light_ice(fls_gam, v = "DrivAge" , n_max = 200 , seed = 3) , alpha = 0.1) +
  ylab("frequency") 
plot(light_ice(fls_rt, v = "DrivAge" , n_max = 200 , seed = 3) , alpha = 0.1) +
  ylab("frequency") 
plot(light_ice(fls_rf, v = "DrivAge" , n_max = 200 , seed = 3) , alpha = 0.1) +
  ylab("frequency") 
plot(light_ice(fls_gbm, v = "DrivAge" , n_max = 200 , seed = 3) , alpha = 0.1) +
  ylab("frequency") 

### INDIVIDUAL CONDITIONAL EXPECTATION SEV
plot(light_ice(fls_glm_sev , v = "DrivAge" , n_max = 200 , seed = 3) , alpha = 0.1) +
  ylab("severity") 
plot(light_ice(fls_gam_sev, v = "Region" , n_max = 200 , seed = 3) , alpha = 0.1) +
  ylab("severity") 
plot(light_ice(fls_rt_sev, v = "Region" , n_max = 200 , seed = 3) , alpha = 0.1) +
  ylab("severity")
plot(light_ice(fls_rf_sev, v = "Region" , n_max = 200 , seed = 3) , alpha = 0.1) +
  ylab("severity")
plot(light_ice(fls_gbm_sev, v = "Region" , n_max = 200 , seed = 3) , alpha = 0.1) +
  ylab("severity")


### PARTIAL DEPENCENCE CURVES
### freq
### drivage
plot(light_profile(fls_ml, v = "DrivAge", type = "predicted") ) +  
  ylab("frequency")

plot(light_profile( fls_glm, v = "DrivAge" , n_bins = 25, type = "predicted"))+  
  ylab("frequency") 

plot(light_profile( fls_gam_2, v = "DrivAge", type = "predicted" ))+  
  ylab("frequency")

### BonusMalus
plot(light_profile( fls_ml, v = "BonusMalus", type = "predicted" )) +  
  ylab("frequency")

plot(light_profile(fls_glm, v = "BonusMalus", type = "predicted")) +  
  ylab("frequency")

plot(light_profile(fls_gam, v = "BonusMalus", type = "predicted" )) +  
  ylab("frequency")

## VehAge
base::plot(light_profile( fls_ml, v = "VehAge", type = "predicted" , pd_evaluate_at = 0:20)) +  
  ylab("frequency")

plot(light_profile(fls_glm, v = "VehAge", type = "predicted")) +  
  ylab("frequency")

base::plot(light_profile(fls_gam, v = "VehAge", pd_evaluate_at = 0:20)) +  
  ylab("frequency") 

### PARTIAL DEPENCENCE CURVES
#### DrivAge
plot(light_profile(fls_ml_sev, v = "DrivAge", type = "predicted") ) +  
  ylab("severity")

plot(light_profile( fls_glm_sev, v = "DrivAge" , n_bins = 25, type = "predicted"))+  
  ylab("severity") 

plot(light_profile( fls_gam_sev, v = "DrivAge", type = "predicted" ))+  
  ylab("severity")

#### VehBrand

plot(light_profile(fls_ml_sev, v = "VehBrand", type = "predicted") ) +  
  ylab("severity")

plot(light_profile( fls_glm_sev, v = "VehBrand" , n_bins = 25, type = "predicted"))+  
  ylab("severity") 

plot(light_profile( fls_gam_sev, v = "VehBrand", type = "predicted" ))+  
  ylab("severity")

