rm(list = ls())


setwd("C:/Users/micho/github/Thesis/Data_01")

setwd("~/github*Thesis/Data_02")

library(tidyr)
library(dplyr)
library(ggplot2)
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

Re_complexify <- function(x){
  sub_z <- gsub(" ", "", sub("im", "i", x))
  
  complex_z <- as.complex(sub_z)
  
  result <- Re(complex_z)
  
  return(result)
}

merge_zero <- function(x){
  eqsmape_min = -3
  eqsmape = x
  if(eqsmape < eqsmape_min){
    eqsmape = eqsmape_min
  }
  
  return(eqsmape)
}

eqsmape_sat <- function(x){
  eqsmape_low = -2
  eqsmape_high = 2
  eqsmape = x
  if(eqsmape < eqsmape_low){
    eqsmape = eqsmape_low
  }
  else if(eqsmape > eqsmape_high){
    eqsmape = eqsmape_high
  }
  
  return(eqsmape)
}

conf_test <- function(x){
  upper_tol = log(1.20)
  lower_tol = log(0.80)
  result <- case_when(
    x <= lower_tol ~ "U",
    lower_tol < x && x < upper_tol ~ "N",
    upper_tol <= x ~ "O"
  )
  
  return(result)
}

sign_match <- function(x, y){
  case_when(
    x*y > 0 ~ "Match",
    x*y < 0 ~ "Mismatch"
  )
}

DE_conf_test <- function(x){
  upper_tol = 1.20
  lower_tol = 0.80
  result <- case_when(
    x <= lower_tol ~ "U",
    lower_tol < x && x < upper_tol ~ "N",
    upper_tol <= x ~ "O"
  )
  
  return(result)
}


# Collect and merge data files
###############################################################################

data_files <- list.files(pattern = "\\.csv$")

results = list()

for(i in seq_along(data_files)){
  results[[i]] <- read.csv(file = data_files[i])
}


results_merged <- bind_rows(results, .id = "column_label")


## Data wrangling and subsetting

#Q <- quantile(results_merged$Eq_SMAPE, probs=c(.05, .95), na.rm = FALSE)
#iqr <- IQR(results_merged$Eq_SMAPE)

#rm_filtered <- subset(results_merged, results_merged$Eq_SMAPE > (Q[1] - 1.5*iqr) & results_merged$Eq_SMAPE < (Q[2]+1.5*iqr))
#rm_filtered <- results_merged
#rm_filtered <- subset(results_merged, results_merged$dtmin_err == "Success")

rm_filtered <- mutate(results_merged, Eq_SMAPE_rs = lapply(Eq_SMAPE, merge_zero))

rm_filtered <- mutate(rm_filtered, SMAPE_test = lapply(SMAPE, conf_test))
rm_filtered <- mutate(rm_filtered, Eq_SMAPE_test = lapply(Eq_SMAPE, conf_test))
rm_filtered <- mutate(rm_filtered, ReDomEig = lapply(domEig, Re_complexify))
rm_filtered <- mutate(rm_filtered, ReDomEigLV = lapply(domEigLV, Re_complexify))
rm_filtered <- mutate(rm_filtered, ReDomEigRatio =
                        as.numeric(ReDomEigLV)/as.numeric(ReDomEig))

rm_filtered <- mutate(rm_filtered, DE_Scomp = mapply(sign_match, ReDomEig, ReDomEigLV))
rm_filtered <- mutate(rm_filtered, DE_test = lapply(ReDomEigRatio, DE_conf_test))

hist(exp(as.numeric(rm_filtered$Eq_SMAPE_rs)))

rm_filtered <- mutate(rm_filtered, Group_NO = lapply(NO, group_by_NO))

rm_filtered <- mutate(rm_filtered, Group_CO = lapply(CO, group_by_CO))

rm_filtered <- mutate(rm_filtered, Group_L = lapply(Leakage, group_by_L))


df_NOut <- as.data.frame(lapply(rm_filtered, unlist))


#########################################################################################################

barplot(df_NOut$SMAPE_test)

ggplot(df_NOut, aes(x=Re(complexify(domEig)), y=Re(complexify(domEigLV)))) +
  geom_point(aes(color = as.factor(Eq_SMAPE_test) ), alpha=0.6) + 
  geom_abline(slope=1, intercept = 0) + 
  theme(text = element_text(size = 24)) + 
  guides(color=guide_legend(title="EqConf")) +
  ylim(-2.5, 2.5) + ylab("LV-DE") + xlab("DE") +
  scale_color_brewer(palette="YlOrRd")

ggplot(df_NOut, aes(y=log(eq_t), x=log(trc_max))) +
  geom_point(aes(color = as.factor(Eq_SMAPE_test) )) + 
  geom_abline(slope=1, intercept = 0) + 
  theme(text = element_text(size = 24)) + 
  guides(color=guide_legend(title="EqConf")) +
  scale_color_brewer(palette="YlOrRd")

################################################################################
df_NOut %>% group_by(SMAPE_test) %>% dplyr::summarise(count = count(Group_L))
df_NOut <- df_NOut %>% count(Group_L)  %>% group_by(Group_L) 

ggplot(data=df_NOut, aes(x = Group_NO, y =..count.. / sum(..count..), fill=DE_Scomp)) + 
  geom_bar(position = "fill") + ylab("Proportion") + xlab("Leakage Group") +
  scale_fill_brewer(palette = "Set1") +
  guides(fill=guide_legend(title="Trajectory Confidence")) + 
  theme(text = element_text(size = 24))

ggplot(data=df_NOut, aes(x = Group_L, y =..count.. / sum(..count..), fill=Eq_SMAPE_test)) + 
  geom_bar(position = "fill") + ylab("Proportion") + xlab("Leakage Group") +
  scale_fill_brewer(palette = "Set1") +
  guides(fill=guide_legend(title="Equilibrium Confidence")) + 
  theme(text = element_text(size = 24))
  
################################################################################

ggplot(data=df_NOut, aes(x = Group_NO, y =..count.. / sum(..count..), fill=Eq_SMAPE_test)) + 
  geom_bar(position = "fill") + ylab("Proportion") + xlab("Niche Overlap Group") +
  scale_fill_brewer(palette = "Set1") +
  guides(fill=guide_legend(title="Equilibrium Confidence")) + 
  theme(text = element_text(size = 24))

ggplot(data=df_NOut, aes(x = Group_NO, y =..count.. / sum(..count..), fill=SMAPE_test)) + 
  geom_bar(position = "fill") + ylab("Proportion") + xlab("Niche Overlap Group") +
  scale_fill_brewer(palette = "Set1") +
  guides(fill=guide_legend(title="Trajectory Confidence")) + 
  theme(text = element_text(size = 24))
  

###############################################################################

modlm <- glm(Eq_SMAPE_rs ~ Leakage*CO*NO, family = gaussian, 
              data=df_NOut)

modlm <- glm(abs(Re(complexify(domEigLV)) - Re(complexify(domEig))) ~ Leakage, family = gaussian, 
             data=df_NOut)


hist(df_NOut$SMAPE)

summary(modlm)
AOV1 <- aov(Eq_SMAPE_rs ~ Group_L + Group_NO + Group_CO, data=df_NOut)
thsd <- TukeyHSD(AOV1)
print(thsd)
plot(thsd)
plot_model(modlm)

VG

nrow(subset(df_NOut, abs(Eq_SMAPE) > 1e-1)) 

ggplot(df_NOut, aes(x = SMAPE, fill=as.factor(Group_L))) +
  geom_density(alpha=0.6) + xlim(c(-1, 5)) + 
  theme(text = element_text(size = 24)) + 
  scale_fill_brewer(palette="YlOrRd") +
  guides(fill=guide_legend(title="Leakage"))

###########################################################################################

df_melt <- melt(df_NOut, measure.vars = c("gR", "gR_LV"))

ggplot(df_melt, aes(x = value)) +
  geom_density(aes(fill=variable), alpha=0.6) + xlim(c(-1, 1)) + 
  theme(text = element_text(size = 24)) +
  scale_fill_brewer(palette = "Spectral") 

ggplot(df_NOut, aes(x = log10(abs(ReDomEigRatio)) )) +
  geom_density(aes(fill=Group_L), alpha = 0.6) +  xlim(c(-5, 20))
  theme(text = element_text(size = 24)) +
  scale_fill_brewer(palette = "Spectral") 



#########################################################################################################

ggplot(df_NOut, aes(x = (Eq_SMAPE), fill=as.factor(Group_L))) +
  geom_density(alpha = 0.6) + xlim(c(-2.5, 2.5)) + 
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

ggplot(df_NOut, aes(Group_CO, Group_L, fill= Eq_SMAPE_sat )) + 
  geom_tile() + 
  scale_fill_gradientn(limits=c(-2, 2), colours=c("navyblue", "darkmagenta", "darkorange1"))

ggplot(df_NOut, aes(Group_NO, Group_L, fill= Eq_SMAPE_sat )) + 
  geom_tile() + 
  scale_fill_gradientn(limits=c(-2, 2), colours=c("navyblue", "darkmagenta", "darkorange1"))

ggplot(df_NOut, aes(Group_L, c=(Group_CO, Group_NO), fill= Eq_SMAPE_sat )) + 
  geom_tile() + 
  scale_fill_gradientn(limits=c(-2, 2), colours=c("navyblue", "darkmagenta", "darkorange1"))

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

ggplot(df_NOut, aes(x = log10(eq_t ), fill=as.factor(dtmin_err))) +
  geom_density(alpha = 0.6) + theme(text = element_text(size = 24))

###########################################################################################

ggplot(df_NOut, aes(x = log(trc_max), fill=as.factor(SMAPE_test))) +
  geom_density(alpha = 0.6) + theme(text = element_text(size = 24)) + xlim(c(-20, 1))

ggplot(df_NOut, aes(x = log10(eq_t ), fill=as.factor(Eq_SMAPE_test))) +
  geom_density(alpha = 0.6) + theme(text = element_text(size = 24)) 

plot(log10(df_NOut$trc_max), log10(df_NOut$eq_t))


cor.test(df_NOut$trc_max, df_NOut$eq_t)
