rm(list = ls())


setwd("./Thesis/Data_final")


library(purrr)
library(tidyr)
library(dplyr)
library(ggplot2)
library(Rmisc)
library(hrbrthemes)
library(ggforce)
library(reshape2)
library(gridExtra)
library(gtable)
library(grid)
library(gtable)
library(ggpubr)
library(latex2exp)


# Declare grouping functions
################################################################################
group_by_L <- function(x){
  #Groups data into leakage groups
  
  
  result <- "NA"
  
  result <- case_when(
    x <= 0.2 ~ "L1",
    0.2 < x && x <= 0.5 ~ "L2",
    0.5 < x  ~ "L3"
  )
  
  return(result)
}

group_by_NO <- function(x){
  #Groups data into Niche Overlap groups
  
  
  result <- "NA"
  
  result <- case_when(
    x <= 0.333 ~ "NO1",
    0.333 < x && x <= 0.666 ~ "NO2",
    0.666 < x ~ "NO3",
  )
  
  return(result)
}

group_by_CO <- function(x){
  #Groups data into Cooperation groups
  
  
  result <- "NA"
  
  result <- case_when(
    x <= 0.309 ~ "CO1",
    0.309 < x && x <= 0.618 ~ "CO2",
    0.618 < x ~ "CO3"
  )
  
  return(result)
}

complexify <- function(x){
  # Re-formats eigenvalues to be compatible with R
  
  
  sub_z <- gsub(" ", "", sub("im", "i", x))
  
  complex_z <- as.complex(sub_z)
  
  return(complex_z)
}

Re_complexify <- function(x){
  # Re-formats real part of eigenvalues to be compatible with R
  
  
  sub_z <- gsub(" ", "", sub("im", "i", x))
  
  complex_z <- as.complex(sub_z)
  
  result <- Re(complex_z)
  
  return(result)
}

merge_zero <- function(x){
  # Cuts off error values at threshold to avoid overflow 
  
  
  eqsmape_min = -3
  eqsmape = x
  if(eqsmape < eqsmape_min){
    eqsmape = eqsmape_min
  }
  
  return(eqsmape)
}

conf_test <- function(x){
  # Assigns simulations to groups based error tolerances
  
  
  upper_tol = log(1.20)
  lower_tol = log(0.80)
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

# Rescale simulation error
rm_filtered <- mutate(results_merged, Eq_SMAPE_rs = lapply(Eq_SMAPE, merge_zero))

c
rm_filtered <- mutate(rm_filtered, SMAPE_test = lapply(SMAPE, conf_test))
rm_filtered <- mutate(rm_filtered, Eq_SMAPE_test = lapply(Eq_SMAPE, conf_test))

# Reformat egienvalues
rm_filtered <- mutate(rm_filtered, ReDomEig = lapply(domEig, Re_complexify))
rm_filtered <- mutate(rm_filtered, ReDomEigLV = lapply(domEigLV, Re_complexify))

# Groups based off of leakage, Niche overlap, and Cooperation
rm_filtered <- mutate(rm_filtered, Group_NO = lapply(NO, group_by_NO))

rm_filtered <- mutate(rm_filtered, Group_CO = lapply(CO, group_by_CO))

rm_filtered <- mutate(rm_filtered, Group_L = lapply(Leakage, group_by_L))


df_NOut <- as.data.frame(lapply(rm_filtered, unlist))


#########################################################################################################

# Stability and Reactivity plots
########################################################################
# Stability
########################################################################
DE_plot <- ggplot(df_NOut, aes(x=Re(complexify(domEig)), y=Re(complexify(domEigLV)) )) +
  geom_point(aes(color = Leakage, size = NO), alpha=0.3) + 
  geom_abline(slope=1, intercept = 0, color = "black") + 
  theme(text = element_text(size = 25)) + 
  guides(color=guide_legend(title="Leakage"), size=guide_legend(title="Niche Overlap")) +
  ylim(-0.2, 0.005) + ylab(TeX(r'($\lambda_{dom}^{R}(A_{LV})$)')) + 
  xlab(TeX(r'($\lambda_{dom}^{R}(A_{Mi})$)')) +
  scale_color_distiller(
    type = "div",
    palette = "RdYlBu"
  )+ scale_size(range = c(0.7, 7)) 

# Reactivity
########################################################################

gR_plor <- ggplot(df_NOut, aes(x=gR, y=gR_LV )) +
  geom_point(aes(color = Leakage, size = NO), alpha=0.7) + 
  geom_abline(slope=1, intercept = 0, color = "black") + 
  theme(text = element_text(size = 25)) + 
  guides(color=guide_legend(title="Leakage"), size=guide_legend(title="Niche Overlap")) +
  ylim(-0.2, 2.1) + ylab(TeX(r'($\lambda_{dom}^{R}(H[A_{LV}])$)')) + xlab(TeX(r'($\lambda_{dom}^{R}(H[A_{Mi}])$)')) +
  scale_color_distiller(
    type = "div",
    palette = "RdYlBu"
  )+ scale_size(range = c(0.7, 7)) 

##############################################################################

ggpubr::ggarrange(DE_plot, gR_plor, ncol = 2, labels = c("A", "B"), common.legend = TRUE, 
                  legend = "bottom")
##############################################################################
g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}
##############################################################################

legend = g_legend(DE_plot)

lay <- rbind(c(1,1,2,2,3),
             c(1,1,2,2,3),
             c(1,1,2,2,3))

grid.arrange(DE_plot+ theme(legend.position="none"),
             gR_plor+ theme(legend.position="none"), 
             legend, ncol = 4, layout_matrix = lay)



###############################################################################

# Timescale separation plots
###############################################################################
trc_eq_plot <- ggplot(df_NOut, aes(x = log10(1/(eq_t*trc_max)), fill=as.factor(Eq_SMAPE_test))) +
  geom_density(alpha = 0.6) + theme(text = element_text(size = 24)) + 
  scale_fill_brewer(palette = "Set1", labels=unname(c('Within', 'Over', 'Under'))) +
  guides(fill=guide_legend(title="Predicted Equilibrium")) +
  xlab(expression(log[10](epsilon^-1/t[eq]))) + 
  geom_vline(xintercept = 0, color = "magenta", linetype="dashed")

trc_plot <- ggplot(df_NOut, aes(x = log10(trc_max), fill=as.factor(Eq_SMAPE_test))) +
  geom_density(alpha = 0.6) + theme(text = element_text(size = 24)) + 
  scale_fill_brewer(palette = "Set1", labels=unname(c('Within', 'Over', 'Under'))) +
  guides(fill=guide_legend(title="Predicted Equilibrium")) +
  xlab(expression(log[10](epsilon))) + 
  geom_vline(xintercept = -2, color = "red", linetype="dashed") +
  geom_vline(xintercept = 0, color = "magenta", linetype="dashed")

legend = g_legend(trc_eq_plot)

lay <- rbind(c(1,1,1,2,2,2,3),
             c(1,1,1,2,2,2,3),
             c(1,1,1,2,2,2,3))

grid.arrange(trc_plot+ theme(legend.position="none"), 
             trc_eq_plot + theme(legend.position="none"),
             legend, ncol = 3, layout_matrix = lay)

ggarrange(trc_plot, trc_eq_plot, ncol=2, nrow=1, common.legend = TRUE, legend="bottom",
          common.axis.x=TRUE, labels = c("A", "B"), font.label = list(size = 20, color = "black"))
###############################################################################



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
  scale_fill_brewer(palette = "Set1", labels=unname(TeX(c(r'($  \in 80pth$ )', r'($>80pth$ )', r'($<80pth$ )')))) +
  guides(fill=guide_legend(title="predicted Eq.")) + 
  theme(text = element_text(size = 24))
  
################################################################################

plots <- df_NOut %>% 
  # Group and nest data by Species 
  # Creates a dataframe with two columns: 'Species' and 'data' (a list-column)
  group_by(Group_NO) %>% 
  nest() %>% 
  # Add a new column with the ggplot objects
  mutate(plots = pmap(.l = list(data, as.factor(Group_NO)), 
                      ~ ggplot(data = ..1) + # first element of .l
                        aes(x = Group_L, # expose the x and y variables
                            y = ..count.. / sum(..count..),
                            fill=Eq_SMAPE_test) +
                        geom_bar(position = "fill") + ylab("Proportion") + xlab("Leakage Group")
                      + labs(title=..2) + scale_fill_brewer(palette = "Set1", 
                                                            labels=unname(c('Within', 'Over', 'Under'))) + 
                        guides(fill=guide_legend(title="Predicted Equilibrium")) + 
                        theme(text = element_text(size = 24))
                        ))

# Walk through the plots column, printing each ggplot object
walk(.x = plots$plots,
     ~ print(.x))

p1 <- plots$plots[[1]] + labs(title="High Niche Overlap")
p2 <- plots$plots[[2]] + labs(title="Med Niche Overlap")
p3 <- plots$plots[[3]] + labs(title="Low Niche Overlap")

ggarrange(p3, p2, p1, ncol=3, nrow=1, common.legend = TRUE, legend="bottom")

g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  legend
}

legend = g_legend(plots$plots[[1]])

lay <- rbind(c(1,1,2,2,3,3,4),
             c(1,1,2,2,3,3,4),
             c(1,1,2,2,3,3,4))

grid.arrange(p3+ theme(legend.position="none"), 
             p2+ theme(legend.position="none"), 
            p1+ theme(legend.position="none"), 
             legend, ncol = 4, layout_matrix = lay)

################################################################################

plots <- df_NOut %>% 
  # Group and nest data by Species 
  # Creates a dataframe with two columns: 'Species' and 'data' (a list-column)
  group_by(Group_CO) %>% 
  nest() %>% 
  # Add a new column with the ggplot objects
  mutate(plots = pmap(.l = list(data, as.factor(Group_CO)), 
                      ~ ggplot(data = ..1) + # first element of .l
                        aes(x = Group_L, # expose the x and y variables
                            y = ..count.. / sum(..count..),
                            fill=Eq_SMAPE_test) +
                        geom_bar(position = "fill") + ylab("Proportion") + xlab("Leakage Group")
                      + labs(title=..2) + scale_fill_brewer(palette = "Set1")
  ))

# Walk through the plots column, printing each ggplot object
walk(.x = plots$plots,
     ~ print(.x))

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
  
ggplot(data=df_NOut, aes(x = Eq_SMAPE_test, y =..count.. / sum(..count..), fill=DE_Scomp)) + 
  geom_bar(position = "fill") + ylab("Proportion") + xlab("Equilibrium Confidence") +
  scale_fill_brewer(palette = "Spectral") +
  guides(fill=guide_legend(title="Stability Match")) + 
  theme(text = element_text(size = 24))

ggplot(data=df_NOut, aes(x = Eq_SMAPE_test, y =..count.. / sum(..count..), fill=DE_test)) + 
  geom_bar(position = "fill") + ylab("Proportion") + xlab("Equilibrium Confidence") +
  scale_fill_brewer(palette = "Spectral") +
  guides(fill=guide_legend(title="Stability Confidence")) + 
  theme(text = element_text(size = 24))

ggplot(data=df_NOut, aes(x = Eq_SMAPE_test, y =..count.. / sum(..count..), fill=gR_Scomp)) + 
  geom_bar(position = "fill") + ylab("Proportion") + xlab("Equilibrium Confidence") +
  scale_fill_brewer(palette = "YlOrRd") +
  guides(fill=guide_legend(title="Reactivity Match")) + 
  theme(text = element_text(size = 24))

ggplot(data=df_NOut, aes(x = Eq_SMAPE_test, y =..count.. / sum(..count..), fill=gR_test)) + 
  geom_bar(position = "fill") + ylab("Proportion") + xlab("Equilibrium Confidence") +
  scale_fill_brewer(palette = "YlOrRd") +
  guides(fill=guide_legend(title="Reactivity Confidence")) + 
  theme(text = element_text(size = 24))
###############################################################################

modlm <- lm(sqrt(abs(Eq_SMAPE_rs)) ~ CV, data=df_NOut)

modlm <- glm(Re(complexify(domEigLV)) ~ Re(complexify(domEig)), 
             data=subset(df_NOut, Eq_SMAPE_test == "N"))


hist(df_NOut$SMAPE)

summary(modlm)
AOV1 <- aov(Eq_SMAPE_rs ~ Group_L + Group_NO + Group_CO, data=df_NOut)
thsd <- TukeyHSD(AOV1)
print(thsd)
plot(thsd)
plot_model(modlm)
plot(modlm)

VG

nrow(subset(df_NOut, abs(Eq_SMAPE) > 1e-1)) 

ggplot(df_NOut, aes(x = SMAPE, fill=as.factor(Group_L))) +
  geom_density(alpha=0.6) + xlim(c(-1, 5)) + 
  theme(text = element_text(size = 24)) + 
  scale_fill_brewer(palette="YlOrRd") +
  guides(fill=guide_legend(title="Leakage"))








###########################################################################################

ggplot(df_NOut, aes(x = log10(trc_max), fill=as.factor(SMAPE_test))) +
  geom_density(alpha = 0.6) + theme(text = element_text(size = 24)) + xlim(c(-7, 1)) + 
  scale_fill_brewer(palette = "Spectral") +
  guides(fill=guide_legend(title="TrajConf")) +
  xlab("log10(return time ratio)") + 
  geom_vline(xintercept = -2, color = "magenta", linetype="dashed") + 
  geom_vline(xintercept = -1, color = "red", linetype="dashed")

ggplot(df_NOut, aes(x = log10(trc_max), fill=as.factor(gR_Scomp))) +
  geom_density(alpha = 0.6) + theme(text = element_text(size = 24)) + 
  scale_fill_brewer(palette = "Spectral") +
  guides(fill=guide_legend(title="Reactivity Match")) +
  xlab("log10(return time ratio)") +
  geom_vline(xintercept = -2, color = "magenta", linetype="dashed") + 
  geom_vline(xintercept = -1, color = "red", linetype="dashed")

ggplot(df_NOut, aes(y=SMAPE, x=log(trc_max))) +
  geom_point(aes(color = as.factor(Eq_SMAPE_test) ), alpha=0.6) + 
  geom_abline(slope=1, intercept = 0) + 
  theme(text = element_text(size = 24)) + 
  guides(color=guide_legend(title="EqConf")) +
  ylab("LV-DE") + xlab("DE") + 
  scale_color_brewer(palette="YlOrRd")

df_NOut %>% count(gR_LV < 0)
