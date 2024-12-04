# Load packages ----------------------------------------------------------------

library(rstan, quietly=T)
library(tidyverse, quietly=T)
library(bayesplot, quietly=T)
library(GGally, quietly=T)
library(ggpubr, quietly=T)
library(scales, quietly=T)
library(cowplot, quietly=T)
library(magrittr, quietly=T)
library(LaplacesDemon, quietly=T)
library(bayestestR, quietly=T)
library(PowerTOST, quietly=T)
library(foreach, quietly=T)
library(doParallel, quietly=T) 
library(PKNCA, quietly=T)
library(kableExtra, quietly=T)
library(data.table, quietly=T)
source("MCSim/function.R")


# Validation -------------------------------------------------------------------
Humansim1.1 <- fread("outputs/iHgHuman_3365.out") 
Humansim2.1 <- fread("outputs/iHgHuman_6734.out") 
Humansim3.1 <- fread("outputs/iHgHuman_4880.out") 
Humansim4.1 <- fread("outputs/iHgHuman_5916.out") 

Human_x <-mcmc_array(list(Humansim1.1, Humansim2.1, Humansim3.1, Humansim4.1))
pars_name <- dimnames(Human_x)[[3]]

str <- which(pars_name == "Ve_Aurine(1)")
end <- which(pars_name == "lnKbrnC(1.2)")
parms <- pars_name[str:end]

j <- seq(150001, 300000, 1) 
sum_chains <- length(j)*4

sample_no = 500 # 500 virtual study
d <- Human_x[j,,] %>% matrix(nrow = sum_chains) %>% as.data.frame() %>% 
  `colnames<-`(dimnames(Human_x)[[3]])
i <- which(d[,"LnPosterior"]==max(d[,"LnPosterior"])) %>%
  rep(sample_no) # Pick the sample based on the maximum a posteriori estimation

tmp.Human_x <- d[i, parms] 
tmp.Human_x |> write.table(file="poppredHuman.dat", row.names=T, sep="\t")

vld <- "./mcsim.iHgHumanBW.model.exe MCSim/iHgHuman_validation.in"
system(vld)

df_Human <- fread("poppredHuman.out") %>% as.data.frame()
dim(df_Human)

predHuman_pop <- fread("poppredHuman.out") |> as.data.frame()
dim(predHuman_pop)
names(predHuman_pop)

# Import Yoshida 1997 data for validation
data_c_Yoshida1997urine <- c(409.0 , 934.6 , 1603.7 , 2516.1 , 3439.6 , 4446.1 , 5015.7 , 5585.3 , 6204.6 , 6707.8 , 7365.9 , 7841.5 , 8449.8 , 9107.8 , 9644.2 , 10263.6 , 11065.4 , 11568.7 , 12530.9 , 12934.6 , 13343.8 , 13813.8 , 14234.1 , 14632.3 , 15074.7 , 15257.1 , 15561.3 , 15826.7 , 16119.8 , 16556.7 , 16827.6 , 17087.6 , 17424.9 , 17612.9 , 17812.0 , 17972.4 , 18177.0 , 18436.9 , 18785.3 , 19034.1 , 19249.8 , 19581.6);
validation_data <- data_c_Yoshida1997urine |> tail(1)

str <- which(names(predHuman_pop)=="CKU_1.1")
predHuman_pop_data <- predHuman_pop[,dim(predHuman_pop)[2]]
organs <- "Urine"
study <- "Yoshida"
dose <- "13750 ug/kg"

# check data fitting
reshape2::melt(predHuman_pop_data) |> mutate(obs = rep(na.omit(validation_data), each=sample_no)) |>
  mutate(ratio = value / obs) |> select(ratio) |> 
  mutate(organs = rep(organs, sample_no), study= rep(study, sample_no), dose=rep(dose, sample_no)) |>
  mutate(acceptance = ifelse(ratio > 3, 0, ifelse(ratio < 0.33, 0, 1))) |>
  mutate(accept_rate = sum(acceptance)/sample_no) |>
  unite("label", study:dose) |>
  mutate(species = "Human") |>
  ggplot() + 
  geom_point(aes(x=ratio, y=label)) +
  scale_x_log10() +
  facet_grid(~organs) +
  geom_vline(xintercept = 1, lty=2) +
  labs(x="Prediction / Observation", y="Study") +
  theme(strip.background = element_blank())

validation_Human <- reshape2::melt(predHuman_pop_data) |> mutate(obs = rep(na.omit(validation_data), each=sample_no)) |>
  mutate(ratio = value / obs) |> select(ratio) |> 
  mutate(organs = rep(organs, sample_no), study= rep(study, sample_no), dose=rep(dose, sample_no)) |>
  mutate(acceptance = ifelse(ratio > 3, 0, ifelse(ratio < 0.33, 0, 1))) |>
  mutate(accept_rate = sum(acceptance)/sample_no) |>
  unite("label", study:dose) |> mutate(species = "Human") 
validation_Human  |> distinct(species, organs, accept_rate) 
save(validation_Human, file = "validation_Human.RData")


