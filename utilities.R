# utilities 

# load the data

library(CASdatasets)
library(distRforest)
library(rpart)
library(tidyverse)
library(DT)


data("freMTPL2freq")
data("freMTPL2sev")


# set exopsure > 1 equal to 1
freMTPL2freq$Exposure <- if_else(freMTPL2freq$Exposure > 1, 1, freMTPL2freq$Exposure)

# set claimNb > 4 equal to 4 
freMTPL2freq[rownames(freMTPL2freq[freMTPL2freq$ClaimNb > 4,] ),2] <- 4


# write a function for the Poisson deviance

Poisson.Deviance <- function(pred, obs){
  
  200*(sum(pred)-sum(obs)+sum(log((obs/pred)^(obs))))/length(pred)
  
}

Gamma.Deviance <- function(obs, pred, alpha){
  
  x <- rep(0, length(obs)) 
  for (i in 1 : length(obs)){
  
    x[i] <- alpha[i] * (((obs[i] - pred[i]) /  pred[i] ) - log(obs[i]/pred[i]))
    
    
  }
  
  return(2 * sum(x))
  
}




