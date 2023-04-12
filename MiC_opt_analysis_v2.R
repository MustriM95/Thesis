rm(list = ls())


setwd("C:/Users/micho/github/Thesis/Data_01")

setwd("~/github*Thesis/Data_02")

library(dplyr)
library(ggplot2)
library(tidyr)
library(Rmisc)
library(hrbrthemes)
library(ggforce)

# Declare grouping functions
################################################################################
group_by_L <- function(x){
  
  result <- "NA"
  
  result <- case_when(
    x <= 0.2 ~ "L1",
    0.2 < x && x <= 0.5 ~ "L2",
    0.5 < x  ~ "L3"
  )
  
  return(result)
}

group_by_NO <- function(x){
  
  result <- "NA"
  
  result <- case_when(
    x <= 0.333 ~ "NO1",
    0.333 < x && x <= 0.666 ~ "NO2",
    0.666 < x ~ "NO3",
  )
  
  return(result)
}

group_by_CO <- function(x){
  
  result <- "NA"
  
  result <- case_when(
    x <= 0.309 ~ "CO1",
    0.309 < x && x <= 0.618 ~ "CO2",
    0.618 < x ~ "CO3"
  )
  
  return(result)
}

complexify <- function(x){
  sub_z <- gsub(" ", "", sub("im", "i", x))
  
  complex_z <- as.complex(sub_z)
  
  return(complex_z)
}

merge_zero <- function(x){
  eqsmape_min = -3
  eqsmape = x
  if(eqsmape < eqsmape_min){
    eqsmape = eqsmape_min
  }
  
  return(eqsmape)
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

#Q <- quantile(results_merged$Eq_SMAPE, probs=c(.05, .95), na.rm = FALSE)
#iqr <- IQR(results_merged$Eq_SMAPE)

#rm_filtered <- subset(results_merged, results_merged$Eq_SMAPE > (Q[1] - 1.5*iqr) & results_merged$Eq_SMAPE < (Q[2]+1.5*iqr))
#rm_filtered <- results_merged
#rm_filtered <- subset(results_merged, results_merged$dtmin_err == "Success")

rm_filtered <- mutate(results_merged, Eq_SMAPE_rs = lapply(Eq_SMAPE, merge_zero))

hist(as.numeric(rm_filtered$Eq_SMAPE_rs))

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
  theme(text = element_text(size = 24)) + 
  scale_color_brewer()

ggplot(df_NOut, aes(x = SMAPE, fill=as.factor(Group_L))) +
  geom_density(alpha=0.6) + xlim(c(-1, 4)) + 
  theme(text = element_text(size = 24)) + 
  scale_fill_brewer() +
  guides(fill=guide_legend(title="Leakage"))


#########################################################################################################

ggplot(df_NOut, aes(x = SMAPE, fill=as.factor(Group_NO))) +
  geom_density(alpha = 0.6) + xlim(c(-1, 4)) + 
  theme(text = element_text(size = 24)) +
  scale_fill_brewer() +
  guides(fill=guide_legend(title="Niche Overlap"))

#########################################################################################################

ggplot(df_NOut, aes(x = SMAPE, fill=as.factor(Group_CO))) +
  geom_density(alpha = 0.6) + xlim(c(-1, 4)) + 
  theme(text = element_text(size = 24)) +
  scale_fill_brewer() +
  guides(fill=guide_legend(title="Cross Feeding"))

#########################################################################################################

ggplot(df_NOut, aes(x = Eq_SMAPE_rs, fill=as.factor(Group_L))) +
  geom_density(alpha = 0.6) + xlim(c(-3, 6)) +
  theme(text = element_text(size = 24)) +
  scale_fill_brewer(palette = "YlOrRd") +
  guides(fill=guide_legend(title="Leakage")) +
  facet_zoom(ylim = c(0,3), shrink = TRUE)

ggplot(subset(df_NOut, Group_L %in% c("L2", "L3")),
       aes(x = Eq_SMAPE_rs, fill=as.factor(Group_L))) +
  xlim(c(-3, 6)) +
  geom_density(alpha = 0.6) +
  scale_fill_brewer() +
  guides(fill=guide_legend(title="Leakage")) +
  facet_zoom(ylim = c(0,3), shrink = TRUE)

count(df_NOut$Group_L == "L4")
#########################################################################################################

ggplot(df_NOut, aes(x = Eq_SMAPE_rs, fill=as.factor(Group_NO))) +
  geom_density(alpha = 0.6) + xlim(c(-3, 6)) +
  theme(text = element_text(size = 24)) +
  scale_fill_brewer() +
  guides(fill=guide_legend(title="Niche Overlap")) +
  facet_zoom(ylim = c(0,1.5), shrink = FALSE)

#########################################################################################################

ggplot(df_NOut, aes(x = Eq_SMAPE, fill=as.factor(Group_CO))) +
  geom_density(alpha = 0.6) + theme(text = element_text(size = 24))

#########################################################################################################

ggplot(df_NOut, aes(x = Re(complexify(domEig)), fill=as.factor(Group_L))) +
  geom_density(alpha = 0.3)  + theme(text = element_text(size = 24))

ggplot(df_NOut, aes(x = Re(complexify(domEigLV)), fill=as.factor(Group_L))) +
  geom_density(alpha = 0.3) + theme(text = element_text(size = 24)) + xlim(c(-1, 2))

#########################################################################################################

ggplot(df_NOut, aes(x = Re(complexify(domEig)), fill=as.factor(Group_NO))) +
  geom_density(alpha = 0.3)  + 
  theme(text = element_text(size = 24)) +
  scale_fill_brewer() +
  guides(fill=guide_legend(title="Niche Overlap")) +
  xlab("MiCRM Dominant Eigenvalue Re")

ggplot(df_NOut, aes(x = Re(complexify(domEigLV)), fill=as.factor(Group_NO))) +
  geom_density(alpha = 0.3) + 
  theme(text = element_text(size = 24)) + 
  xlim(c(-0.5, 1.5)) +
  scale_fill_brewer() +
  guides(fill=guide_legend(title="Niche Overlap")) +
  xlab("LV Dominant Eigenvalue Re")

#########################################################################################################

ggplot(df_NOut, aes(x = Re(complexify(domEig)), fill=as.factor(Group_CO))) +
  geom_density(alpha = 0.3)  + 
  theme(text = element_text(size = 24))+
  scale_fill_brewer() +
  guides(fill=guide_legend(title="Cooperation")) +
  xlab("MiCRM Dominant Eigenvalue Re")

ggplot(df_NOut, aes(x = Re(complexify(domEigLV)), fill=as.factor(Group_CO))) +
  geom_density(alpha = 0.3) + xlim(c(-0.5, 2)) +
  theme(text = element_text(size = 24)) +
  scale_fill_brewer() +
  guides(fill=guide_legend(title="Niche Overlap")) +
  xlab("LV Dominant Eigenvalue Re")

cor.test(df_NOut$NO, Re(complexify(df_NOut$domEig)))
cor.test(df_NOut$NO, Re(complexify(df_NOut$domEigLV)))

cor.test(df_NOut$CO, Re(complexify(df_NOut$domEig)))
cor.test(df_NOut$CO, Re(complexify(df_NOut$domEigLV)))


#########################################################################################################



#########################################################################################################

ggplot(df_NOut, aes(x = gR, fill=as.factor(Group_L))) +
  geom_density(alpha = 0.3) +
  theme(text = element_text(size = 24)) +
  scale_fill_brewer() +
  guides(fill=guide_legend(title="Leakage")) +
  xlab("MiCRM Reactivity")
  

ggplot(df_NOut, aes(x = gR_LV, fill=as.factor(Group_L))) +
  geom_density(alpha = 0.3) + xlim(c(-1, 1.5)) +
  theme(text = element_text(size = 24)) +
  scale_fill_brewer() +
  guides(fill=guide_legend(title="Leakage")) +
  xlab("LV Reactivity")

#########################################################################################################

ggplot(df_NOut, aes(x = gR, fill=as.factor(Group_NO))) +
  geom_density(alpha = 0.3) + 
  theme(text = element_text(size = 24)) +
  scale_fill_brewer() +
  guides(fill=guide_legend(title="Niche Overlap")) +
  xlab("MiCRM Reactivity")

ggplot(df_NOut, aes(x = gR_LV, fill=as.factor(Group_NO))) +
  geom_density(alpha = 0.3) + xlim(c(-1, 1)) +
  theme(text = element_text(size = 24)) +
  scale_fill_brewer() +
  guides(fill=guide_legend(title="Niche Overlap")) +
  xlab("LV Reactivity")

#########################################################################################################

ggplot(df_NOut, aes(x = gR, fill=as.factor(Group_CO))) +
  geom_density(alpha = 0.3) +
  theme(text = element_text(size = 24)) +
  scale_fill_brewer() +
  guides(fill=guide_legend(title="Cooperation")) +
  xlab("MiCRM Reactivity")

ggplot(df_NOut, aes(x = gR_LV, fill=as.factor(Group_CO))) +
  geom_density(alpha = 0.3) + theme(text = element_text(size = 24))



##########################################################################################

max(df_NOut$gR)
min(df_NOut$gR_LV)

ggplot(df_NOut, aes(Group_CO, Leakage, fill= gR )) + 
  geom_tile() + scale_fill_gradientn(colours=c("navyblue", "darkmagenta", "darkorange1"))

ggplot(df_NOut, aes(Group_CO, Leakage, fill= gR_LV )) + 
  geom_tile() + scale_fill_gradientn(colours=c("navyblue", "darkmagenta", "darkorange1"))


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

ggplot(df_NOut, aes(Group_CO, Group_L, fill= Eq_SMAPE_rs )) + 
  geom_tile() + scale_fill_gradientn(colours=c("navyblue", "darkmagenta", "darkorange1"))

ggplot(df_NOut, aes(Group_NO, Group_L, fill= Eq_SMAPE_rs )) + 
  geom_tile() + scale_fill_gradientn(colours=c("navyblue", "darkmagenta", "darkorange1"))

ggplot(df_NOut, aes(Group_NO, Group_CO, fill= Eq_SMAPE )) + 
  geom_tile() + scale_fill_gradientn(colours=c("navyblue", "darkmagenta", "darkorange1"))

###########################################################################################

ggplot(df_NOut, aes(Group_CO, Leakage, fill= Re(complexify(domEig)) )) + 
  geom_tile() + scale_fill_gradientn(colours=c("navyblue", "darkmagenta", "darkorange1"))

ggplot(df_NOut, aes(Group_CO, Leakage, fill= Re(complexify(domEigLV)) )) + 
  geom_tile() + scale_fill_gradientn(colours=c("navyblue", "darkmagenta", "darkorange1"))

ggplot(df_NOut, aes(Group_NO, Leakage, fill= Re(complexify(domEig)) )) + 
  geom_tile() + scale_fill_gradientn(colours=c("navyblue", "darkmagenta", "darkorange1"))

ggplot(df_NOut, aes(Group_NO, Leakage, fill= Re(complexify(domEigLV)) )) + 
  geom_tile() + scale_fill_gradientn(colours=c("navyblue", "darkmagenta", "darkorange1"))

ggplot(df_NOut, aes(Group_CO, Group_NO, fill= Re(complexify(domEig)) )) + 
  geom_tile() + scale_fill_gradientn(colours=c("navyblue", "darkmagenta", "darkorange1"))

ggplot(df_NOut, aes(Group_CO, Group_NO, fill= Re(complexify(domEigLV)) )) + 
  geom_tile() + scale_fill_gradientn(colours=c("navyblue", "darkmagenta", "darkorange1"))

###########################################################################################

ggplot(df_NOut, aes(Group_NO, Leakage, fill=log10(eq_t ))) + 
  geom_tile() + scale_fill_gradientn(colours=c("navyblue", "darkmagenta", "darkorange1"))

ggplot(df_NOut, aes(Group_NO, Leakage, fill= log(trc_max ))) + 
  geom_tile() + scale_fill_gradientn(colours=c("navyblue", "darkmagenta", "darkorange1"))

###########################################################################################

ggplot(df_NOut, aes(Group_CO, Leakage, fill=log10(eq_t ) )) + 
  geom_tile() + scale_fill_gradientn(colours=c("navyblue", "darkmagenta", "darkorange1"))

ggplot(df_NOut, aes(Group_CO, Leakage,, fill= log10(trc_max ))) + 
  geom_tile() + scale_fill_gradientn(colours=c("navyblue", "darkmagenta", "darkorange1"))

###########################################################################################

ggplot(df_NOut, aes(Group_CO, Group_NO, fill=log10(eq_t ) )) + 
  geom_tile() + scale_fill_gradientn(colours=c("navyblue", "darkmagenta", "darkorange1"))

ggplot(df_NOut, aes(Group_CO, Group_NO, fill= log(trc_max ))) + 
  geom_tile() + scale_fill_gradientn(colours=c("navyblue", "darkmagenta", "darkorange1"))



###########################################################################################

ggplot(df_NOut, aes(x = log10(trc_max), fill=as.factor(dtmin_err))) +
  geom_density(alpha = 0.6) + theme(text = element_text(size = 24))

ggplot(df_NOut, aes(x = log10(eq_t ), fill=as.factor(Group_NO))) +
  geom_density(alpha = 0.6) + theme(text = element_text(size = 24))

###########################################################################################

ggplot(df_NOut, aes(x = log(trc_max), fill=as.factor(Group_CO))) +
  geom_density(alpha = 0.6) + theme(text = element_text(size = 24))

ggplot(df_NOut, aes(x = log10(eq_t ), fill=as.factor(Group_CO))) +
  geom_density(alpha = 0.6) + theme(text = element_text(size = 24))

plot(log10(df_NOut$trc_max), df_NOut$Eq_SMAPE)


