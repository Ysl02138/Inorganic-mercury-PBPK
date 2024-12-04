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
# Mouse -------------------------------------------------------------------
Mousesim1.1 <- fread("outputs/iHgMouse_3365.out") 
Mousesim2.1 <- fread("outputs/iHgMouse_6734.out") 
Mousesim3.1 <- fread("outputs/iHgMouse_4880.out") 
Mousesim4.1 <- fread("outputs/iHgMouse_5916.out") 

Mouse_x <-mcmc_array(list(Mousesim1.1, Mousesim2.1, Mousesim3.1, Mousesim4.1))
pars_name <- dimnames(Mouse_x)[[3]]
str <- which(pars_name == "Ve_Aurine(1)")
end <- which(pars_name == "lnKbrnC(1.4)")
parms <- pars_name[str:end]

j <- seq(50001, 100000, 10) 
sum_chains <- length(j)*4

n = 100 # 100 virtual study
d <- Mouse_x[j,,] %>% matrix(nrow = sum_chains) %>% as.data.frame() %>% 
  `colnames<-`(dimnames(Mouse_x)[[3]])
i <- which(d[,"LnPosterior"]==max(d[,"LnPosterior"])) %>%
  rep(n) # Pick the sample based on the maximum a posteriori estimation

tmp.Mouse_x <- d[i, parms] 
tmp.Mouse_x %>% write.table(file="Mouse.HED.dat", row.names=T, sep="\t")

vld <- "./mcsim.iHgMouseBW.model.exe MCSim/Mouse_HED_MAP.in"
system(vld)

Mouse_out <- read.delim("Mouse.HED.out")
Mouse_out = as.data.frame(Mouse_out)

# Tidy data with median and 95% confidence interval
vars <- names(Mouse_out)
index <- which(vars == "AUCCL_1.1" | vars == "CKU_1.1")
Mouse_Summary <- apply(Mouse_out[index[1]:index[2]], 2, quantile,  c(0.5, 0.025, 0.975)) %>% t()
colnames(Mouse_Summary) <- c("median", "LCL", "UCL")
df_Mouse <- as.data.frame(Mouse_Summary)


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

j <- seq(100001, 200000, 10) 
sum_chains <- length(j)*4

n = 100 # 100 virtual study
d <- Human_x[j,,] %>% matrix(nrow = sum_chains) %>% as.data.frame() %>% 
  `colnames<-`(dimnames(Human_x)[[3]])
i <- which(d[,"LnPosterior"]==max(d[,"LnPosterior"])) %>%
  rep(n) # Pick the sample based on the maximum a posteriori estimation

tmp.Human_x <- d[i, parms] 
tmp.Human_x %>% write.table(file="Human.HED.dat", row.names=T, sep="\t")

vld <- "./mcsim.iHgHumanBW.model.exe MCSim/Human_HED_MAP_Mouse.in"
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

## Mouse 
liver118       = ((Mouse_out$AUCCL_1.1)/(Human_out$AUCCL_1.1))*118    # Liver,  Mouse / Human, LOAEL 118 ug/kg/d
kidney118      = ((Mouse_out$AUCCK_1.1)/(Human_out$AUCCK_1.1))*118    # Kidney, Mouse / Human, LOAEL 118 ug/kg/d
brain118       = ((Mouse_out$AUCCBrn_1.1)/(Human_out$AUCCBrn_1.1))*118# Brain,  Mouse / Human, LOAEL 118 ug/kg/d
blood118       = ((Mouse_out$AUCCBld_1.1)/(Human_out$AUCCBld_1.1))*118# Blood,  Mouse / Human, LOAEL 118 ug/kg/d

# Summary output 
df_Mouse <-data.frame(Mouse_HED_liver=c(liver118), Mouse_HED_kidney= c(kidney118), Mouse_HED_brain =c(brain118), Mouse_HED_blood = c(blood118)) 
sub_df_Mouse <-df_Mouse[,c('Mouse_HED_liver','Mouse_HED_kidney','Mouse_HED_brain','Mouse_HED_blood')] 
Mouse_HED <- apply(sub_df_Mouse, 2, function(x) quantile(x,probs = c(0.025, 0.5, 0.975))) %>% as.data.frame
write.csv(Mouse_HED, "outputs/Mouse_HED.csv", row.names=TRUE)

# tidy
print(Mouse_blood_HED   <- quantile(blood118 , probs = c(0.025, 0.5, 0.975)))
print(Mouse_liver_HED   <- quantile(liver118, probs = c(0.025, 0.5, 0.975)))
print(Mouse_kidney_HED  <- quantile(kidney118, probs = c(0.025, 0.5, 0.975)))
print(Mouse_brain_HED   <- quantile(brain118, probs = c(0.025, 0.5, 0.975)))







