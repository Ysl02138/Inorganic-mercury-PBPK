library(rstan)
library(tidyverse)
library(bayesplot)
library(GGally)
library(ggpubr)
library(scales)
library(cowplot)
library(magrittr)
library(LaplacesDemon)
library(bayestestR)
library(PowerTOST)
library(foreach)
library(doParallel) 
library(PKNCA)
library(kableExtra)
library(data.table)
library(formattable)
source("MCSim/function.R")

###########################################################################
# Rat -------------------------------------------------------------------
Ratsim1.1 <- fread("outputs/iHgRat_3365.out") 
Ratsim2.1 <- fread("outputs/iHgRat_6734.out") 
Ratsim3.1 <- fread("outputs/iHgRat_4880.out") 
Ratsim4.1 <- fread("outputs/iHgRat_5916.out") 

Rat_x <-mcmc_array(list(Ratsim1.1, Ratsim2.1, Ratsim3.1, Ratsim4.1))
pars_name <- dimnames(Rat_x)[[3]]
str <- which(pars_name == "Ve_Aurine(1)")
end <- which(pars_name == "lnKbrnC(1.4)")
parms <- pars_name[str:end]

j <- seq(50001, 100000, 10) 
sum_chains <- length(j)*4

n = 100 # 100 virtual study
d <- Rat_x[j,,] %>% matrix(nrow = sum_chains) %>% as.data.frame() %>% 
  `colnames<-`(dimnames(Rat_x)[[3]])
i <- which(d[,"LnPosterior"]==max(d[,"LnPosterior"])) %>%
  rep(n) # Pick the sample based on the maximum a posteriori estimation

tmp.Rat_x <- d[i, parms] 
tmp.Rat_x %>% write.table(file="Rat.HED.dat", row.names=T, sep="\t")

vld <- "./mcsim.iHgRatBW.model.exe MCSim/Rat_HED_MAP.in"
system(vld)

Rat_out <- read.delim("Rat.HED.out")

# Tidy data with median and 95% confidence interval
vars <- names(Rat_out)
index <- which(vars == "AUCCL_1.1" | vars == "CKU_1.1")
Rat_Summary <- apply(Rat_out[index[1]:index[2]], 2, quantile,  c(0.5, 0.025, 0.975)) %>% t()
colnames(Rat_Summary) <- c("median", "LCL", "UCL")
df_Rat <- as.data.frame(Rat_Summary)


###########################################################################
# Human -------------------------------------------------------------------
Humansim1.1 <- fread("outputs/iHgHuman_3365.out") 
Humansim2.1 <- fread("outputs/iHgHuman_6734.out") 
Humansim3.1 <- fread("outputs/iHgHuman_4880.out") 
Humansim4.1 <- fread("outputs/iHgHuman_5916.out") 

Human_x <-mcmc_array(list(Humansim1.1, Humansim2.1, Humansim3.1, Humansim4.1))
pars_name <- dimnames(Human_x)[[3]]
str <- which(pars_name == "Ve_Aurine(1)")
end <- which(pars_name == "lnKbrnC(1.2)")
parms <- pars_name[str:end]

j <- seq(150001, 300000, 10) 
sum_chains <- length(j)*4

n = 100 # 100 virtual study
d <- Human_x[j,,] %>% matrix(nrow = sum_chains) %>% as.data.frame() %>% 
  `colnames<-`(dimnames(Human_x)[[3]])
i <- which(d[,"LnPosterior"]==max(d[,"LnPosterior"])) %>%
  rep(n) # Pick the sample based on the maximum a posteriori estimation

tmp.Human_x <- d[i, parms] 
tmp.Human_x %>% write.table(file="Human.HED.dat", row.names=T, sep="\t")

vld <- "./mcsim.iHgHumanBW.model.exe MCSim/Human_HED_MAP_Rat.in"
system(vld)

Human_out <- read.delim("Human.HED.out")
Human_out = as.data.frame(Human_out)
# Tidy data with median and 95% confidence interval
vars <- names(Human_out)
index <- which(vars == "AUCCL_1.1" | vars == "CKU_1.1")
Human_Summary <- apply(Human_out[index[1]:index[2]], 2, quantile,  c(0.5, 0.025, 0.975)) %>% t()
colnames(Human_Summary) <- c("median", "LCL", "UCL")
df_Human <- as.data.frame(Human_Summary)


#####################################################################
## Calculated the HED 

## Rat 
liver330       = ((Rat_out$AUCCL_1.1)/(Human_out$AUCCL_1.1))*330    # Liver,  Rat / Human, NOAEL 330 ug/kg/d
kidney330      = ((Rat_out$AUCCK_1.1)/(Human_out$AUCCK_1.1))*330    # Kidney, Rat / Human, NOAEL 330 ug/kg/d
brain330       = ((Rat_out$AUCCBrn_1.1)/(Human_out$AUCCBrn_1.1))*330# Brain,  Rat / Human, NOAEL 330 ug/kg/d
blood330       = ((Rat_out$AUCCBld_1.1)/(Human_out$AUCCBld_1.1))*330# Blood,  Rat / Human, NOAEL 330 ug/kg/d

# Summary output 
df_Rat <-data.frame(Rat_HED_liver=c(liver330), Rat_HED_kidney= c(kidney330), Rat_HED_brain =c(brain330), Rat_HED_blood = c(blood330)) 
sub_df_Rat <-df_Rat[,c('Rat_HED_liver','Rat_HED_kidney','Rat_HED_brain','Rat_HED_blood')] 
Rat_HED <- apply(sub_df_Rat, 2, function(x) quantile(x,probs = c(0.025, 0.5, 0.975))) %>% as.data.frame
write.csv(Rat_HED, "outputs/Rat_HED.csv", row.names=TRUE)

print(Rat_blood_HED   <- quantile(blood330,  probs = c(0.025, 0.5, 0.975)))
print(Rat_liver_HED   <- quantile(liver330,  probs = c(0.025, 0.5, 0.975)))
print(Rat_kidney_HED  <- quantile(kidney330, probs = c(0.025, 0.5, 0.975)))
print(Rat_brain_HED   <- quantile(brain330,  probs = c(0.025, 0.5, 0.975)))






