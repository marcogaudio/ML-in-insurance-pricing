# model explainability

library(tidyverse)
library(MetricsWeighted)
library(flashlight)
library(ggplot2)

# in order to compare the explainability of the models, these are fitted over
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
processed_dataset[,"VehAge"] <-relevel(processed_dataset[,"VehAge"], ref="2")
DrivAge <- cbind(c(18:100), c(rep(1,3), rep(2,5), rep(3,5), rep(4,10), rep(5,10), rep(6,20), rep(7,30)))
processed_dataset$DrivAge <- as.factor(DrivAge[processed_dataset$DrivAge-17,2])
processed_dataset[,"DrivAge"] <-relevel(processed_dataset[,"DrivAge"], ref="5")
processed_dataset$BonusMalus <- as.integer(pmin(processed_dataset$BonusMalus, 150))
processed_dataset$Density <- as.numeric(log(processed_dataset$Density))
processed_dataset[,"Region"] <-relevel(processed_dataset[,"Region"], ref="Centre")

#FREQUENCY

# glm model



glm_model <- glm(ClaimNb ~ VehPower + VehAge + DrivAge + BonusMalus
                 + VehBrand + VehGas + Density + Region + Area, 
                 data = processed_dataset, offset=log(Exposure), family=poisson())



fl_glm <- flashlight(
   model = glm_model , label = "GLM" ,
   predict_function = function(fit , X) predict(fit , X , type = "response")
  )



# GAM

library(gam)
library(mgcv)
library(mgcViz)


# data needs to be preprocessed before fitting a poisson glm with log link
# (following paper about french case)
processed_dataset2 <- dataset
processed_dataset2$VehGas <- factor(processed_dataset2$VehGas)
# processed_dataset2$VehPowerGLM <- as.factor(pmin(processed_dataset2$VehPower,9))
processed_dataset2$AreaGAM <- as.integer(processed_dataset2$Area)
processed_dataset2$BonusMalusGAM <- as.integer(pmin(processed_dataset2$BonusMalus, 150))
processed_dataset2$DensityGAM <- as.numeric(log(processed_dataset2$Density))
processed_dataset2[,"Region"] <-relevel(processed_dataset2[,"Region"], ref="Centre")

gam_model <- gam::gam(ClaimNb ~ gam::s(VehPower) + gam::s(VehAge) + gam::s(DrivAge) + gam::s(BonusMalus)
                      + VehBrand + VehGas + gam::s(Density) + Region + Area, 
                      data = processed_dataset2, offset=log(Exposure), family=poisson())

fl_gam <- flashlight (
  model = gam_model , label = "GAM" ,
  predict_function = function(fit , X) gam::predict.Gam(fit, X, "response")
)


metrics <- list("Average deviance" = deviance_poisson ,
                 "Relative deviance reduction" = r_squared_poisson )

fls_glm <- multiflashlight (list(fl_glm), data = processed_dataset ,
                            y = "ClaimNb", w = "Exposure", metrics = metrics)
fls_gam <- multiflashlight (list(fl_gam), data = processed_dataset2 ,
                            y = "ClaimNb", w = "Exposure", metrics = metrics)


x <- c("VehPower", "VehAge", "DrivAge", "BonusMalus",
       "VehBrand", "VehGas", "Density", "Region", "Area")
imp_glm <- light_importance(fls_glm, v = x)
plot(imp_glm, fill = "#E69F00", color = "black")

imp_gam <- light_importance(fls_gam, v = x)
plot(imp_gam, fill = "#E69F00", color = "black")







