rm(list = ls())


setwd("C:/Users/micho/github/Thesis")

library(dplyr)
library(ggplot2)
library(tidyr)

df = read.csv("3x3_newTest.csv")

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
    x <= 0.1 ~ "NO1",
    0.1 < x && x <= 0.2 ~ "NO2",
    0.2 < x && x <= 0.3 ~ "NO3",
    0.3 < x && x <= 0.4 ~ "NO4",
    0.4 < x && x <= 0.5 ~ "NO5",
    0.5 < x && x <= 0.6 ~ "NO6",
    0.6 < x && x <= 0.7 ~ "NO7",
    0.7 < x && x <= 0.8 ~ "NO8",
    0.8 < x && x <= 0.9 ~ "NO9",
    0.9 < x && x <= 1.0 ~ "NO10"
  )
  
  return(result)
}

group_by_CO <- function(x){
  
  result <- "NA"
  
  result <- case_when(
    x <= 0.05 ~ "CO1",
    0.05 < x && x <= 0.1 ~ "CO2",
    0.1 < x && x <= 0.15 ~ "CO3",
    0.15 < x && x <= 0.2 ~ "CO4",
    0.2 < x && x <= 0.25 ~ "CO5",
    0.25 < x && x <= 0.3 ~ "CO6",
    0.3 < x && x <= 0.35 ~ "CO7",
    0.35 < x && x <= 0.4 ~ "CO8",
    0.4 < x  ~ "CO9"
  )
  
  return(result)
}

## Data wrangling and subsetting

Q <- quantile(df$SMAPE, probs=c(.25, .75), na.rm = FALSE)
iqr <- IQR(df$SMAPE)

df_NOut <- subset(df, df$SMAPE > (Q[1] - 1.5*iqr) & df$SMAPE < (Q[2]+1.5*iqr))

df_NOut <- mutate(df, Group_NO = lapply(NO, group_by_NO))

df_NOut <- mutate(df_NOut, CO_L = CO/(1 + Leakage))

df_NOut <- mutate(df_NOut, Group_CO = lapply(CO_L, group_by_CO))

df_NOut <- mutate(df_NOut, Group_L = lapply(Leakage, group_by_L))


df_NOut$Group_NO <- factor(df_NOut$Group_NO, levels=c("NO1","NO2","NO3", "NO4", "NO5", 
                                                      "NO6", "NO7", "NO8", "NO9", "NO10"))

df_NOut <- as.data.frame(lapply(df_NOut, unlist))

#########################################################################################################


## Niche Overlap analysis

ggplot(data=df_NOut, aes(y=SMAPE, x=as.factor(Group_NO))) + geom_violin() +
  labs(title = "Mean Squared Error vs. Niche Overlap", y = "SMAPE", x="Niche Overlap Group", color = "Leakage") + 
  geom_boxplot(width=0.1) + theme(text = element_text(size = 18))

ggplot(data=df_NOut, aes(y=Eq_SMAPE, x=as.factor(Group_NO))) + geom_violin() +
  labs(title = "Equilibrium MSE vs. Niche Overlap", y = "Eq_SMAPE", x="Niche Overlap Group", color = "Leakage") + 
  geom_boxplot(width=0.1) + theme(text = element_text(size = 18))

ggplot(data=df_NOut, aes(y=-sqrt(abs(domEig)), x=as.factor(Group_NO))) + geom_violin() +
  labs(title = "MiCRM dominant Eigenvalue vs. Niche Overlap", y = "log(-domEig)", x="Niche Overlap Group", color = "Leakage") + 
  geom_boxplot(width=0.1) + theme(text = element_text(size = 18))

ggplot(data=df_NOut, aes(y=log(abs(domEigLV)), x=as.factor(Group_NO))) + geom_violin() +
  labs(title = "LVM dominant Eigenvalue vs. Niche Overlap", y = "log(abs(domEigLV))", x="Niche Overlap Group") + 
  geom_boxplot(width=0.1) + theme(text = element_text(size = 18))

ggplot(data=df_NOut, aes(y=log(abs(domEig - domEigLV)), x=as.factor(Group_NO))) + geom_violin() +
  labs(title = "Dominant eigenvalue difference vs. NO", y = "log(domEigDiff)", x="Niche Overlap Group", color = "Leakage") + 
  geom_boxplot(width=0.1) + theme(text = element_text(size = 18))



#########################################################################################################


## Leakage analysis

ggplot(data=df_NOut, aes(y=SMAPE, x=as.factor(Leakage))) + geom_violin() +
  labs(title = "Mean Squared Error vs. Leakage", y = "log(MSE)", x="Leakage") + 
  geom_boxplot(width=0.1) + theme(text = element_text(size = 18))

ggplot(data=df_NOut, aes(y=Eq_SMAPE, x=as.factor(Leakage))) + geom_violin() +
  labs(title = "Equilibrium Mean Squared Error vs. Leakage", y = "log(Eq_MSE)", x="Leakage") + 
  geom_boxplot(width=0.1) + theme(text = element_text(size = 18))

ggplot(data=df_NOut, aes(y=log(abs(domEig)), x=as.factor(Leakage))) + geom_violin() +
  labs(title = "MiCRM Dominant eigenvalues vs. Leakage", y = "log(-domEig)", x="Leakage") + 
  geom_boxplot(width=0.1) + theme(text = element_text(size = 18))

ggplot(data=df_NOut, aes(y=log(abs(domEigLV)), x=as.factor(Leakage))) + geom_violin() +
  labs(title = "LVM dominant Eigenvalue vs. Leakage", y = "log(abs(domEigLV))", x="Leakage") + 
  geom_boxplot(width=0.1) + theme(text = element_text(size = 18))

ggplot(data=df_NOut, aes(y=log(abs(domEig - domEigLV)), x=as.factor(Leakage))) + geom_violin() +
  labs(title = "Dominant eigenvalue difference vs. Leakage", y = "log(domEigDiff)", x="Leakage") + 
  geom_boxplot(width=0.1) + theme(text = element_text(size = 18))



#########################################################################################################

## Cooperation analysis

ggplot(data=df_NOut, aes(y=SMAPE, x=as.factor(Group_CO))) + geom_violin() +
  labs(title = "Mean Squared Error vs. Cooperation", y = "SMAPE", x="Cooperation") + 
  geom_boxplot(width=0.1) + theme(text = element_text(size = 18))

ggplot(data=df_NOut, aes(y=Eq_SMAPE, x=as.factor(Group_CO))) + geom_violin() +
  labs(title = "Equilibrium MSE vs. Cooperation", y = "Eq_SMAPE", x="Cooperation") + 
  geom_boxplot(width=0.1) + theme(text = element_text(size = 18))

ggplot(data=df_NOut, aes(y=log(abs(domEig)), x=as.factor(Group_CO))) + geom_violin() +
  labs(title = "MiCRM dominant eigenvalue vs. Cooperation", y = "log(-domEig)", x="Cooperation") + 
  geom_boxplot(width=0.1) + theme(text = element_text(size = 18))

ggplot(data=df_NOut, aes(y=log(abs(domEigLV)), x=as.factor(Group_CO))) + geom_violin() +
  labs(title = "LVM dominant Eigenvalue vs. Cooperation", y = "log(abs(domEigLV))", x="Cooperation", color = "Leakage") + 
  geom_boxplot(width=0.1) + theme(text = element_text(size = 18))

ggplot(data=df_NOut, aes(y=log(abs(domEig-domEigLV)), x=as.factor(Group_CO))) + geom_violin() +
  labs(title = "Dominant eigenvalue difference vs. Cooperation", y = "log(domEigDiff)", x="Cooperation") + 
  geom_boxplot(width=0.1) + theme(text = element_text(size = 18))

#########################################################################################################
