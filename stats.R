library(CASdatasets)
library(distRforest)
library(rpart)
library(tidyverse)
library(DT)

data("freMTPLfreq")
data("freMTPLsev")


# qui ci sono più regioni
data("freMTPL2freq")
data("freMTPL2sev")

view(freMTPLfreq)


tot_claim <- freMTPLfreq %>%  
  group_by(ClaimNb) %>% 
  summarise(n_policy = n()) 

tot_exposure <- freMTPLfreq %>% 
  mutate(Exposure = if_else(Exposure > 1,1,Exposure)) %>%  # there are many exposure >1 due to errors, correct them
  group_by(ClaimNb) %>%  
  summarise(tot_epxosure = sum(Exposure))
    
freMTPLfreq$Exposure <- if_else(freMTPLfreq$Exposure > 1, 1, freMTPLfreq$Exposure)

# overvie of the portfolio
table1 <- left_join(tot_claim, tot_exposure, by = c("ClaimNb")) %>% 
  mutate(relative_nb = n_policy/sum(n_policy))


# absolute claim number histogram (cambiare y axis)
KULbg <- "#116E8A"
g <- ggplot(freMTPLfreq, aes(ClaimNb)) + theme_bw() + 
  geom_bar(col = KULbg, fill = KULbg) + 
  labs(y = "Abs frequency") +
  ggtitle("MTPL - number of claims") 
g


# relative claim nb


KULbg <- "#116E8A"
g <- ggplot(table1, aes(x = ClaimNb, y = relative_nb)) + theme_bw() + 
  geom_bar(col = KULbg, fill = KULbg, stat = "identity") + 
  labs(y = "relative frequency", x = "Number of claims") +
  ggtitle("French MTPL -Relative frequency") 
g


#exposure histogram

p <- ggplot(freMTPLfreq, aes(x = Exposure)) + theme_bw() + 
  geom_histogram(col = KULbg, fill = KULbg, aes(y = stat(count)/sum(count)), bins = 15) + 
  scale_y_continuous(labels = scales::percent) +
  labs(y = "relative frequency") + 
  ggtitle("French MTPL - Exposure") 
p

# density plot

d <- ggplot(freMTPLsev, aes(x = ClaimAmount)) +
  geom_density(col = KULbg, fill = KULbg, alpha=0.8) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  xlim(0, 7000) +
  labs(x = "Claim amount") +
  ggtitle("French MTPL - Claim amount density")
  
d

# overall severity

sum(freMTPLsev$ClaimAmount) / sum(freMTPLfreq$ClaimNb)


table(freMTPLfreq$Exposure)/length(freMTPLfreq$Exposure)

# relative frequency CarAge 

p <- ggplot(freMTPLfreq, aes(x = CarAge)) + theme_bw() + 
  geom_histogram(col = KULbg, fill = KULbg, aes(y = stat(count)/sum(count)), stat = "count") + 
  scale_y_continuous(labels = scales::percent) +
  labs(y = "relative frequency") + 
  ggtitle("French MTPL - CarAge") 
p

#relative frequency driverAge

p <- ggplot(freMTPLfreq, aes(x = DriverAge)) + theme_bw() + 
  geom_histogram(col = KULbg, fill = KULbg, aes(y = stat(count)/sum(count)), stat = "count") + 
  scale_y_continuous(labels = scales::percent) +
  labs(y = "relative frequency") + 
  ggtitle("French MTPL - DriverAge") 
p

# relative freq gas
p <- ggplot(freMTPLfreq, aes(x = Gas)) + theme_bw() + 
  geom_histogram(col = KULbg, fill = KULbg, aes(y = stat(count)/sum(count)), stat = "count") + 
  scale_y_continuous(labels = scales::percent) +
  labs(y = "relative frequency") + 
  ggtitle("French MTPL - Gas") 
p

# almost 50-50 


# relative freq region
p <- ggplot(freMTPLfreq, aes(x = Region)) + theme_bw() + 
  geom_histogram(col = KULbg, fill = KULbg, aes(y = stat(count)/sum(count)), stat = "count") + 
  scale_y_continuous(labels = scales::percent) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(y = "relative frequency") + 
  ggtitle("French MTPL - Region") 
p


driver_perc <- freMTPLfreq %>% 
  filter(DriverAge <=70,
         DriverAge>=25) %>% 
  count()


tot_driver <- freMTPLfreq %>%  
  count()

driver_perc/tot_driver

# 88.97% of the policyholders are aged between 25 and 70


# car age

car_perc <- freMTPLfreq %>% 
  filter(CarAge <=15,
         CarAge>=0) %>% 
  count()


tot_car <- freMTPLfreq %>%  
  count()

car_perc/tot_car


# percentage people centre
centre <- freMTPLfreq %>% 
  filter(
         Region == "Ile-de-France") %>% 
  count()


centre/tot_driver


