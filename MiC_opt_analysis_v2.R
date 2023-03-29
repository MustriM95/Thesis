rm(list = ls())


setwd("C:/Users/micho/github/Thesis/Data_01")

library(dplyr)
library(ggplot2)
library(tidyr)
library(Rmisc)
library(hrbrthemes)


df3 = read.csv("3x3.csv")
df4 = read.csv("4x4.csv")
df5 = read.csv("5x5.csv")
df6 = read.csv("6x6.csv")
df7 = read.csv("7x7.csv")
df8 = read.csv("8x8.csv")
df9 = read.csv("9x9.csv")
df10 = read.csv("10x10.csv")
df11 = read.csv("11x11.csv")
df12 = read.csv("12x12.csv")


df = rbind(df3, df4, df5, df6, df7, df8, df9, df10, df11, df12)

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
    x <= 0.33 ~ "NO1",
    0.33 < x && x <= 0.66 ~ "NO2",
    0.66 < x  ~ "NO3"
  )
  
  return(result)
}

group_by_CO <- function(x){
  
  result <- "NA"
  
  result <- case_when(
    x <= 0.50 ~ "CO1",
    0.50 < x && x <= 0.70 ~ "CO2",
    0.70 < x  ~ "CO3"
  )
  
  return(result)
}

## Data wrangling and subsetting

Q <- quantile(df$Eq_SMAPE, probs=c(.05, .95), na.rm = FALSE)
iqr <- IQR(df$Eq_SMAPE)

df_NOut <- subset(df, df$Eq_SMAPE > (Q[1] - 1.5*iqr) & df$Eq_SMAPE < (Q[2]+1.5*iqr))

df_NOut <- mutate(df_NOut, Group_NO = lapply(NO, group_by_NO))

df_NOut <- mutate(df_NOut, Group_CO = lapply(CO, group_by_CO))

df_NOut <- mutate(df_NOut, Group_L = lapply(Leakage, group_by_L))


df_NOut$Group_NO <- factor(df_NOut$Group_NO, levels=c("NO1","NO2","NO3", "NO4", "NO5", 
                                                      "NO6", "NO7", "NO8", "NO9", "NO10"))

df_NOut <- as.data.frame(lapply(df_NOut, unlist))


#########################################################################################################

hist(df_NOut$NO)

cor.test(df_NOut$CV, df_NOut$Eq_SMAPE)

ggplot(df_NOut, aes(x = SMAPE, fill=as.factor(Group_L))) +
  geom_density(alpha=0.6) + xlim(c(-5, 5))


plots <- list()
elss = unique(df$Leakage)

for(i in 1:length(elss)){
  title<-capture.output(cat("Leakage =", elss[i]))
  plots[[i]] <- ggplot(df_NOut[which(df_NOut$Leakage==elss[i]),], aes(x=SMAPE)) + 
    geom_histogram(alpha=0.6, binwidth=0.25) + xlim(c(-5, 5)) + ggtitle(title)
}

multiplot(plotlist=plots, cols = 2)

#########################################################################################################

ggplot(df_NOut, aes(x = SMAPE, fill=as.factor(Group_NO))) +
  geom_density(alpha = 0.6) + xlim(c(-2.5, 2.5))



#########################################################################################################

ggplot(df_NOut, aes(x = Eq_SMAPE, fill=as.factor(Group_L))) +
  geom_density(alpha = 0.6) 


#########################################################################################################

ggplot(df_NOut, aes(x = Eq_SMAPE, fill=as.factor(Group_NO))) +
  geom_density(alpha = 0.6) 

#########################################################################################################

ggplot(df_NOut, aes(x = Eq_SMAPE, fill=as.factor(Group_CO))) +
  geom_density(alpha = 0.6) 

ggplot(df_NOut, aes(x = SMAPE, fill=as.factor(Group_CO))) +
  geom_density(alpha = 0.6) 

#########################################################################################################

complexify <- function(x){
  sub_z <- gsub(" ", "", sub("im", "i", x))
  
  complex_z <- as.complex(sub_z)
  
  return(complex_z)
}


ggplot(df_NOut, aes(x = Re(complexify(domEigLV)), fill=as.factor(Group_NO))) +
  geom_density(alpha = 0.3) + xlim(-0.1, 0.1)

ggplot(df_NOut, aes(x = Re(complexify(domEig)), fill=as.factor(Group_NO))) +
  geom_density(alpha = 0.3) + xlim(-0.1, 0.1)

hist(df_NOut$gR)
hist(df_NOut$gR_LV)

ggplot(df_NOut, aes(x = gR, fill=as.factor(Group_CO))) +
  geom_density(alpha = 0.3)

ggplot(df_NOut, aes(x = gR_LV, fill=as.factor(Group_CO))) +
  geom_density(alpha = 0.3)


#########################################################################################################

## Playing around with data

plot(Re(complexify(df$domEig)), Im(complexify(df$domEig)))


#########################################################################################################

# Community size

ggplot(data=df_NOut, aes(y=SMAPE, x=as.factor(N))) + geom_violin() +
  geom_boxplot(width=0.1) 

ggplot(data=df_NOut, aes(y=Eq_SMAPE, x=as.factor(N))) + geom_violin() +
  geom_boxplot(width=0.1) 


##########################################################################################

max(df_NOut$gR)
min(df_NOut$gR_LV)

ggplot(df_NOut, aes(Group_CO, Group_L, fill= gR )) + 
  geom_tile() + scale_fill_gradientn(limits = c(-0.2,2.7), colours=c("navyblue", "darkmagenta", "darkorange1"))

ggplot(df_NOut, aes(x = gR, fill=as.factor(Group_NO))) +
  geom_density(alpha = 0.6) 

ggplot(df_NOut, aes(x = gR_LV, fill=as.factor(Group_NO))) +
  geom_density(alpha = 0.6) 

ggplot(df_NOut, aes(x = gR, fill=as.factor(Group_CO))) +
  geom_density(alpha = 0.6) 

ggplot(df_NOut, aes(x = log(trc_max), fill=as.factor(Group_NO))) +
  geom_density(alpha = 0.6) 

ggplot(df_NOut, aes(x = eq_t, fill=as.factor(Group_NO))) +
  geom_density(alpha = 0.6) 

ggplot(df_NOut, aes(x = gR, fill=as.factor(Group_L))) +
  geom_density(alpha = 0.6) 

ggplot(df_NOut, aes(x = gR_LV, fill=as.factor(Group_L))) +
  geom_density(alpha = 0.6)
df$