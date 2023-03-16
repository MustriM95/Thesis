rm(list = ls())

setwd("~/github*Thesis")
setwd("C:/Users/micho/github/Thesis")

library(dplyr)
library(ggplot2)
library(tidyr)



df3 = read.csv("3x3_newTest.csv")


df = rbind(df2, df3, df4, df5)

hist(df$SMAPE)

hist(df$Eq_SMAPE)

hist(df$NO)

hist(df$domEig)

boxplot(log(df$trc_max))

head(df)

Q <- quantile(df$SMAPE, probs=c(.25, .75), na.rm = FALSE)
iqr <- IQR(df$SMAPE)

df_NOut <- subset(df, df$SMAPE > (Q[1] - 1.5*iqr) & df$SMAPE < (Q[2]+1.5*iqr))


ggplot(data=df_NOut, aes(y=df$SMAPE, x = NO)) + geom_point(aes(color = factor(Group_L))) + 
  labs(title = "Mean Squared Error vs. Niche Overlap", color ="Leakage") + theme(text = element_text(size = 18))

ggplot(data=df_NOut, aes(y=log(Eq_MSE), x = NO)) + geom_point(aes(color = factor(Group_L))) + 
  labs(title = "Equilibrium Mean Squared Error vs. Niche Overlap", color ="Leakage") +
  theme(text = element_text(size = 18))




ggplot(data=df_NOut, aes(y=log(MSE), x=as.factor(Leakage))) + geom_violin() +
  labs(title = "Mean Squared Error vs. Leakage", y = "log(MSE)", x="Leakage", color = "Leakage") + 
  geom_boxplot(width=0.1) + theme(text = element_text(size = 18))

ggplot(data=df_NOut, aes(y=log(Eq_MSE), x=as.factor(Leakage))) + geom_violin() +
  labs(title = "Equilibrium Mean Squared Error vs. Leakage", y = "log(Eq_MSE)", x="Leakage", color = "Leakage") + 
  geom_boxplot(width=0.1) + theme(text = element_text(size = 18))




ggplot(data=df_NOut, aes(y=domEig, x = NO)) + geom_point(aes(color = factor(Group_L))) + 
  labs(title="Dominant Eigenvalue vs Niche Overlap", y="Re(DomEig)", color="Leakage") + 
  theme(text = element_text(size = 18))

## Cooperation analysis

ggplot(data=df_NOut, aes(y=domEig, x = CO)) + geom_point(aes(color = factor(Group_L))) + 
  labs(title="Dominant Eigenvalue vs Cooperation", y="Re(DomEig)", color="Leakage") + 
  theme(text = element_text(size = 18))

ggplot(data=df_NOut, aes(y=log(MSE), x = CO)) + geom_point(aes(color = factor(Group_L))) + 
  labs(title = "Mean Squared Error vs. Cooperation", color ="Leakage") + theme(text = element_text(size = 18))

ggplot(data=df_NOut, aes(y=log(Eq_MSE), x = CO)) + geom_point(aes(color = factor(Group_L))) + 
  labs(title = "Equilibrium Mean Squared Error vs. Cooperation", color ="Leakage") +
  theme(text = element_text(size = 18))






ggplot(data=df_NOut, aes(y=domEigLV, x = NO)) + geom_point(aes(color = factor(Leakage))) + scale_y_continuous(limits = c(-10, 10))

ggplot(data=df_NOut, aes(y=domEig, x = as.factor(N))) + geom_violin() + 
  geom_boxplot(width=0.1)

ggplot(data=df_NOut, aes(y=domEig, x = as.factor(Leakage))) + geom_violin() +
  labs(title="Dominant Eigenvalue vs Leakage", y="Re(DomEig)", x="Leakage") + 
  theme(text = element_text(size = 18)) + geom_boxplot(width=0.1)

ggplot(data=df_NOut, aes(y=log(MSE), x = as.factor(Group_NO))) + geom_violin()

ggplot(data=df_NOut, aes(y=log(trc_max), x = log(MSE))) + geom_point(aes(color = factor(Leakage))) + 
  scale_y_continuous(limits = c(-30, 0))

group_by_L <- function(x){
  
  result <- "NA"
  
  if(x <= 0.2){
    result <- "L1"
  } else if(x <= 0.5){
    result <- "L2"
  } else{
    result <- "L3"
  }
  
  return(result)
}

group_by_NO <- function(x){
  
  result <- "NA"
  
  result <- case_when(
    x <= 0.2 ~ "NO1",
    0.2 < x && x <= 0.4 ~ "NO2",
    0.4 < x && x <= 0.6 ~ "NO3",
    0.6 < x && x <= 0.8 ~ "NO4",
    0.8 < x && x <= 0.9 ~ "NO5"
  )
  
  return(result)
}




df_NOut <- mutate(df_NOut, domEigDiff = abs(domEig - domEigLV))

df_NOut <- mutate(df_NOut, Group_NO = lapply(Leakage, group_by_NO))

df_NOut <- mutate(df_NOut, Group_L = lapply(Leakage, group_by_L))

df_NOut <- as.data.frame(lapply(df_NOut, unlist))

dfL1 <- subset(df_NOut, Leakage == df$Leakage[6])

ggplot(data = dfL1, aes(y=log(MSE), x = NO)) + geom_point()

ggplot(data = df_NOut, aes(y=domEig, x = as.factor(Leakage))) + geom_violin()
