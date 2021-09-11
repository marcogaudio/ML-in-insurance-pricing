# GAM implementation
library(gam)
library(mgcv)
library(mgcViz)

k <- 5 

# data needs to be preprocessed before fitting a poisson glm with log link
# (following paper about french case)
processed_dataset2 <- dataset
processed_dataset2$VehGas <- factor(processed_dataset2$VehGas)
processed_dataset2$AreaGAM <- as.integer(processed_dataset2$Area)
processed_dataset2$BonusMalusGAM <- as.integer(pmin(processed_dataset2$BonusMalus, 150))
processed_dataset2$DensityGAM <- as.numeric(log(processed_dataset2$Density))
processed_dataset2[,"Region"] <-relevel(processed_dataset2[,"Region"], ref="Centre")

gam_model <- mgcv::gam(ClaimNb ~ s(VehPower) + s(VehAge) + s(DrivAge) + s(BonusMalus)
                 + VehBrand + VehGas + Density + Region + Area, 
                 data = processed_dataset2, offset=log(Exposure), family=poisson())

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

