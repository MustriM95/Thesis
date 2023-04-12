rm(list = ls())


setwd("C:/Users/micho/github/Thesis/Data_01")

library(dplyr)
library(ggplot2)
library(tidyr)
library(Rmisc)
library(hrbrthemes)

# Declare grouping functions
################################################################################
group_by_L <- function(x){
  
  result <- "NA"
  
  result <- case_when(
    x <= 0.2 ~ "L1",
    0.2 < x && x <= 0.4 ~ "L2",
    0.4 < x && x <= 0.6  ~ "L3",
    0.6 < x ~ "L4"
  )
  
  return(result)
}

group_by_NO <- function(x){
  
  result <- "NA"
  
  result <- case_when(
    x <= 0.25 ~ "NO1",
    0.25 < x && x <= 0.5 ~ "NO2",
    0.5 < x && x <= 0.75 ~ "NO3",
    0.75 < x  ~ "NO4"
  )
  
  return(result)
}

group_by_CO <- function(x){
  
  result <- "NA"
  
  result <- case_when(
    x <= 0.25 ~ "CO1",
    0.25 < x && x <= 0.5 ~ "CO2",
    0.5 < x && x <= 0.75 ~ "CO3",
    0.75 < x  ~ "CO4"
  )
  
  return(result)
}

complexify <- function(x){
  sub_z <- gsub(" ", "", sub("im", "i", x))
  
  complex_z <- as.complex(sub_z)
  
  return(complex_z)
}

# Collect and merge data files
###############################################################################

data_files <- list.files(pattern = "\\.csv$")

results = list()

for(i in seq_along(data_files)){
  results[[i]] <- read.csv(file = data_files[i])
}


results_merged <- bind_rows(results, .id = "column_label")


hist(results_merged$Eq_SMAPE)




## Data wrangling and subsetting

Q <- quantile(results_merged$Eq_SMAPE, probs=c(.05, .95), na.rm = FALSE)
iqr <- IQR(results_merged$Eq_SMAPE)

rm_filtered <- subset(results_merged, results_merged$Eq_SMAPE > (Q[1] - 1.5*iqr) & results_merged$Eq_SMAPE < (Q[2]+1.5*iqr))

rm_filtered <- subset(results_merged, results_merged$dtmin_err == "Success")

hist(rm_filtered$Eq_SMAPE)
hist(rm_filtred2$Eq_SMAPE)

rm_filtered <- mutate(rm_filtered, Group_NO = lapply(NO, group_by_NO))

rm_filtered <- mutate(rm_filtered, Group_CO = lapply(CO, group_by_CO))

rm_filtered <- mutate(rm_filtered, Group_L = lapply(Leakage, group_by_L))


df_NOut <- as.data.frame(lapply(rm_filtered, unlist))


#########################################################################################################

hist(df_NOut$NO)

cor.test(df_NOut$NO, df_NOut$SMAPE)

ggplot(df_NOut, aes(x=Re(complexify(domEig)), y=Re(complexify(domEigLV)))) +
  geom_point(aes(color = as.factor(Group_CO) )) + 
  geom_abline(slope=1, intercept = 0) + 
  theme(text = element_text(size = 24))

ggplot(df_NOut, aes(x = SMAPE, fill=as.factor(Group_L))) +
  geom_density(alpha=0.6) + xlim(c(-1, 3)) + theme(text = element_text(size = 24))


#########################################################################################################

ggplot(df_NOut, aes(x = SMAPE, fill=as.factor(Group_NO))) +
  geom_density(alpha = 0.6) + xlim(c(-2.5, 2.5)) + theme(text = element_text(size = 24))

#########################################################################################################

ggplot(df_NOut, aes(x = SMAPE, fill=as.factor(Group_CO))) +
  geom_density(alpha = 0.6) + xlim(c(-2.5, 2.5)) + theme(text = element_text(size = 24))



#########################################################################################################

ggplot(df_NOut, aes(x = Eq_SMAPE, fill=as.factor(dtmin_err))) +
  geom_density(alpha = 0.6) + theme(text = element_text(size = 24))


#########################################################################################################

ggplot(df_NOut, aes(x = Eq_SMAPE, fill=as.factor(Group_NO))) +
  geom_density(alpha = 0.6) + theme(text = element_text(size = 24))

#########################################################################################################

ggplot(df_NOut, aes(x = Eq_SMAPE, fill=as.factor(Group_CO))) +
  geom_density(alpha = 0.6) + theme(text = element_text(size = 24))

#########################################################################################################

ggplot(df_NOut, aes(x = Re(complexify(domEig)), fill=as.factor(Group_L))) +
  geom_density(alpha = 0.3)  + theme(text = element_text(size = 24))

ggplot(df_NOut, aes(x = Re(complexify(domEigLV)), fill=as.factor(Group_L))) +
  geom_density(alpha = 0.3) + theme(text = element_text(size = 24))

#########################################################################################################

ggplot(df_NOut, aes(x = Re(complexify(domEig)), fill=as.factor(Group_NO))) +
  geom_density(alpha = 0.3)  + theme(text = element_text(size = 24))

ggplot(df_NOut, aes(x = Re(complexify(domEigLV)), fill=as.factor(Group_NO))) +
  geom_density(alpha = 0.3) + theme(text = element_text(size = 24))

#########################################################################################################

ggplot(df_NOut, aes(x = Re(complexify(domEig)), fill=as.factor(Group_CO))) +
  geom_density(alpha = 0.3)  + theme(text = element_text(size = 24))

ggplot(df_NOut, aes(x = Re(complexify(domEigLV)), fill=as.factor(Group_CO))) +
  geom_density(alpha = 0.3) + theme(text = element_text(size = 24))

cor.test(df_NOut$NO, Re(complexify(df_NOut$domEig)))
cor.test(df_NOut$NO, Re(complexify(df_NOut$domEigLV)))

#########################################################################################################



#########################################################################################################

ggplot(df_NOut, aes(x = gR, fill=as.factor(dtmin_err))) +
  geom_density(alpha = 0.3) + theme(text = element_text(size = 24))

ggplot(df_NOut, aes(x = gR_LV, fill=as.factor(Group_L))) +
  geom_density(alpha = 0.3) + theme(text = element_text(size = 24))

#########################################################################################################

ggplot(df_NOut, aes(x = gR, fill=as.factor(Group_NO))) +
  geom_density(alpha = 0.3) + theme(text = element_text(size = 24))

ggplot(df_NOut, aes(x = gR_LV, fill=as.factor(Group_NO))) +
  geom_density(alpha = 0.3) + theme(text = element_text(size = 24))

#########################################################################################################

ggplot(df_NOut, aes(x = gR, fill=as.factor(Group_CO))) +
  geom_density(alpha = 0.3) + theme(text = element_text(size = 24))

ggplot(df_NOut, aes(x = gR_LV, fill=as.factor(Group_CO))) +
  geom_density(alpha = 0.3) + theme(text = element_text(size = 24))



##########################################################################################

max(df_NOut$gR)
min(df_NOut$gR_LV)

ggplot(df_NOut, aes(Group_CO, Group_L, fill= gR )) + 
  geom_tile() + scale_fill_gradientn(limits = c(-0.2,2.7), colours=c("navyblue", "darkmagenta", "darkorange1"))

ggplot(df_NOut, aes(Group_CO, Group_L, fill= gR_LV )) + 
  geom_tile() + scale_fill_gradientn(limits = c(-0.2,2.7), colours=c("navyblue", "darkmagenta", "darkorange1"))


###########################################################################################

ggplot(df_NOut, aes(Group_CO, Group_L, fill= gR )) + 
  geom_tile() + scale_fill_gradientn(colours=c("navyblue", "darkmagenta", "darkorange1"))

ggplot(df_NOut, aes(Group_CO, Group_L, fill= gR_LV )) + 
  geom_tile() + scale_fill_gradientn(colours=c("navyblue", "darkmagenta", "darkorange1"))

###########################################################################################

ggplot(df_NOut, aes(Group_NO, Group_L, fill= gR )) + 
  geom_tile() + scale_fill_gradientn(colours=c("navyblue", "darkmagenta", "darkorange1"))

ggplot(df_NOut, aes(Group_NO, Group_L, fill= gR_LV )) + 
  geom_tile() + scale_fill_gradientn(colours=c("navyblue", "darkmagenta", "darkorange1"))

###########################################################################################

ggplot(df_NOut, aes(Group_NO, Group_CO, fill= gR )) + 
  geom_tile() + scale_fill_gradientn(colours=c("navyblue", "darkmagenta", "darkorange1"))

ggplot(df_NOut, aes(Group_NO, Group_CO, fill= gR_LV )) + 
  geom_tile() + scale_fill_gradientn(colours=c("navyblue", "darkmagenta", "darkorange1"))


###########################################################################################

ggplot(df_NOut, aes(Group_CO, Group_L, fill= Eq_SMAPE )) + 
  geom_tile() + scale_fill_gradientn(colours=c("navyblue", "darkmagenta", "darkorange1"))

ggplot(df_NOut, aes(Group_NO, Group_L, fill= Eq_SMAPE )) + 
  geom_tile() + scale_fill_gradientn(colours=c("navyblue", "darkmagenta", "darkorange1"))

ggplot(df_NOut, aes(Group_NO, Group_CO, fill= Eq_SMAPE )) + 
  geom_tile() + scale_fill_gradientn(colours=c("navyblue", "darkmagenta", "darkorange1"))

###########################################################################################

ggplot(df_NOut, aes(Group_CO, Group_L, fill= Re(complexify(domEig)) )) + 
  geom_tile() + scale_fill_gradientn(colours=c("navyblue", "darkmagenta", "darkorange1"))

ggplot(df_NOut, aes(Group_CO, Group_L, fill= Re(complexify(domEigLV)) )) + 
  geom_tile() + scale_fill_gradientn(colours=c("navyblue", "darkmagenta", "darkorange1"))

ggplot(df_NOut, aes(Group_NO, Group_L, fill= Re(complexify(domEig)) )) + 
  geom_tile() + scale_fill_gradientn(colours=c("navyblue", "darkmagenta", "darkorange1"))

ggplot(df_NOut, aes(Group_NO, Group_L, fill= Re(complexify(domEigLV)) )) + 
  geom_tile() + scale_fill_gradientn(colours=c("navyblue", "darkmagenta", "darkorange1"))

ggplot(df_NOut, aes(Group_CO, Group_NO, fill= Re(complexify(domEig)) )) + 
  geom_tile() + scale_fill_gradientn(colours=c("navyblue", "darkmagenta", "darkorange1"))

ggplot(df_NOut, aes(Group_CO, Group_NO, fill= Re(complexify(domEigLV)) )) + 
  geom_tile() + scale_fill_gradientn(colours=c("navyblue", "darkmagenta", "darkorange1"))

###########################################################################################

ggplot(df_NOut, aes(Group_NO, Group_L, fill=eq_t )) + 
  geom_tile() + scale_fill_gradientn(colours=c("navyblue", "darkmagenta", "darkorange1"))

ggplot(df_NOut, aes(Group_NO, Group_L, fill= log(trc_max ))) + 
  geom_tile() + scale_fill_gradientn(colours=c("navyblue", "darkmagenta", "darkorange1"))

###########################################################################################

ggplot(df_NOut, aes(Group_CO, Group_L, fill=eq_t )) + 
  geom_tile() + scale_fill_gradientn(colours=c("navyblue", "darkmagenta", "darkorange1"))

ggplot(df_NOut, aes(Group_CO, Group_L, fill= log10(trc_max ))) + 
  geom_tile() + scale_fill_gradientn(colours=c("navyblue", "darkmagenta", "darkorange1"))

###########################################################################################

ggplot(df_NOut, aes(Group_CO, Group_NO, fill=eq_t )) + 
  geom_tile() + scale_fill_gradientn(colours=c("navyblue", "darkmagenta", "darkorange1"))

ggplot(df_NOut, aes(Group_CO, Group_NO, fill= log(trc_max ))) + 
  geom_tile() + scale_fill_gradientn(colours=c("navyblue", "darkmagenta", "darkorange1"))



###########################################################################################

ggplot(df_NOut, aes(x = log10(trc_max), fill=as.factor(dtmin_err))) +
  geom_density(alpha = 0.6) + theme(text = element_text(size = 24))

ggplot(df_NOut, aes(x = eq_t, fill=as.factor(Group_NO))) +
  geom_density(alpha = 0.6) + theme(text = element_text(size = 24))

###########################################################################################

ggplot(df_NOut, aes(x = log(trc_max), fill=as.factor(Group_CO))) +
  geom_density(alpha = 0.6) + theme(text = element_text(size = 24))

ggplot(df_NOut, aes(x = eq_t, fill=as.factor(Group_CO))) +
  geom_density(alpha = 0.6) + theme(text = element_text(size = 24))

plot(log10(df_NOut$trc_max), df_NOut$Eq_SMAPE)


