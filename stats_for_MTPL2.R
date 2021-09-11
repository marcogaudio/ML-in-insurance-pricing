library(CASdatasets)
library(distRforest)
library(rpart)
library(tidyverse)
library(DT)

data("freMTPL2freq")
data("freMTPL2sev")

view(freMTPL2freq)


tot_claim <- freMTPL2freq %>%  
  group_by(ClaimNb) %>% 
  summarise(n_policy = n()) 

tot_exposure <- freMTPL2freq %>% 
  mutate(Exposure = if_else(Exposure > 1,1,Exposure)) %>%  # there are many exposure >1 due to errors, correct them
  group_by(ClaimNb) %>%  
  summarise(tot_epxosure = sum(Exposure))


# overview of the portfolio
table1 <- left_join(tot_claim, tot_exposure, by = c("ClaimNb")) %>% 
  mutate(relative_nb = n_policy/sum(n_policy))


# relative Exposure

table(freMTPL2freq$Exposure)/length(freMTPL2freq$Exposure)

# total claim freq

sum(freMTPL2freq$ClaimNb) / sum(freMTPL2freq$Exposure)

# absolute claim number histogram (cambiare y axis)
KULbg <- "#116E8A"
g <- ggplot(freMTPL2freq, aes(ClaimNb)) + theme_bw() + 
  geom_bar(col = KULbg, fill = KULbg) + 
  labs(y = "Abs frequency") +
  ggtitle("MTPL - number of claims") 
g


# relative claim nb


KULbg <- "#116E8A"
g <- ggplot(table1, aes(x = ClaimNb, y = relative_nb)) + theme_bw() + 
  geom_bar(col = KULbg, fill = KULbg, stat = "identity") + 
  labs(y = "relative frequency", x = "Number of claims") +
  ggtitle("French MTPL - number of claims") 
g


#exposure histogram

p <- ggplot(freMTPL2freq, aes(x = Exposure)) + theme_bw() + 
  geom_histogram(col = KULbg, fill = KULbg, aes(y = stat(count)/sum(count)), bins = 15) + 
  scale_y_continuous(labels = scales::percent) +
  labs(y = "relative frequency") + 
  ggtitle("French MTPL - Exposure") 
p

# density plot

d <- ggplot(freMTPL2sev, aes(x = ClaimAmount)) +
  geom_density(col = KULbg, fill = KULbg, alpha=0.8) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  xlim(0, 7000) +
  labs(x = "Claim amount") +
  ggtitle("Claim amounts")
d

#99th quantile 

quantile(freMTPL2sev$ClaimAmount, probs = c(0.99))

# overall severity

sum(freMTPL2sev$ClaimAmount) / sum(freMTPL2freq$ClaimNb)


table(freMTPL2freq$Exposure)/length(freMTPL2freq$Exposure)

# relative frequency CarAge 

p <- ggplot(freMTPL2freq, aes(x = VehAge)) + theme_bw() + 
  geom_histogram(col = KULbg, fill = KULbg, aes(y = stat(count)/sum(count)), stat = "count") + 
  scale_y_continuous(labels = scales::percent) +
  labs(y = "relative frequency") + 
  ggtitle("French MTPL - Vehicle age") 
p

#relative frequency driverAge

p <- ggplot(freMTPL2freq, aes(x = DrivAge)) + theme_bw() + 
  geom_histogram(col = KULbg, fill = KULbg, aes(y = stat(count)/sum(count)), stat = "count") + 
  scale_y_continuous(labels = scales::percent) +
  labs(y = "relative frequency") + 
  ggtitle("French MTPL - Driver age") 
p

# relative freq gas
p <- ggplot(freMTPL2freq, aes(x = VehGas)) + theme_bw() + 
  geom_histogram(col = KULbg, fill = KULbg, aes(y = stat(count)/sum(count)), stat = "count") + 
  scale_y_continuous(labels = scales::percent) +
  labs(y = "relative frequency") + 
  ggtitle("French MTPL - Gas") 
p

# almost 50-50 


# relative freq region
p <- ggplot(freMTPL2freq, aes(x = Region)) + theme_bw() + 
  geom_histogram(col = KULbg, fill = KULbg, aes(y = stat(count)/sum(count)), stat = "count") + 
  scale_y_continuous(labels = scales::percent) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(y = "relative frequency") + 
  ggtitle("French MTPL - Region") 
p


driver_perc <- freMTPL2freq %>% 
  filter(DrivAge <=70,
         DrivAge>=25) %>% 
  count()


tot_driver <- freMTPL2freq %>%  
  count()

driver_perc/tot_driver

# 88.97% of the policyholders are aged between 25 and 70


# car age

car_perc <- freMTPL2freq %>% 
  filter(VehAge <=15,
         VehAge>=0) %>% 
  count()


tot_car <- freMTPL2freq %>%  
  count()

car_perc/tot_car


# percentage people centre
centre <- freMTPL2freq %>% 
  filter(
    Region == "Centre") %>% 
  count()


centre/tot_driver

# gas veh

Diesel <- freMTPL2freq %>% 
  filter(
    VehGas == "Diesel") %>% 
  count()

Diesel/tot_driver


# vehpower relative freq
p <- ggplot(freMTPL2freq, aes(x = VehPower)) + theme_bw() + 
  geom_histogram(col = KULbg, fill = KULbg, aes(y = stat(count)/sum(count)), stat = "count") + 
  scale_y_continuous(labels = scales::percent) +
  labs(y = "relative frequency") + 
  ggtitle("French MTPL - Vehicle power") 
p






