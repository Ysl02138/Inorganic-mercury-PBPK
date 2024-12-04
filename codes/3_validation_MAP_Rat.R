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
Ratsim1.1 <- fread("outputs/iHgRat_3365.out") 
Ratsim2.1 <- fread("outputs/iHgRat_6734.out") 
Ratsim3.1 <- fread("outputs/iHgRat_4880.out") 
Ratsim4.1 <- fread("outputs/iHgRat_5916.out") 

Rat_x <-mcmc_array(list(Ratsim1.1, Ratsim2.1, Ratsim3.1, Ratsim4.1))
pars_name <- dimnames(Rat_x)[[3]]

str <- which(pars_name == "Ve_Aurine(1)")
end <- which(pars_name == "lnKbrnC(1.4)")

parms <- pars_name[str:end]

j <- seq(50001, 100000, 1) 
sum_chains <- length(j)*4

sample_no = 500 # 500 virtual study
d <- Rat_x[j,,] %>% matrix(nrow = sum_chains) %>% as.data.frame() %>% 
  `colnames<-`(dimnames(Rat_x)[[3]])
i <- which(d[,"LnPosterior"]==max(d[,"LnPosterior"])) %>%
  rep(sample_no) # Pick the sample based on the maximum a posteriori estimation

tmp.Rat_x <- d[i, parms] 
tmp.Rat_x |> write.table(file="poppredRat.dat", row.names=T, sep="\t")

vld <- "./mcsim.iHgRatBW.model.exe MCSim/iHgRat_validation.in"
system(vld)

df_Rat <- fread("poppredRat.out") %>% as.data.frame()
dim(df_Rat)

predRat_pop <- fread("poppredRat.out") |> as.data.frame()
dim(predRat_pop)
names(predRat_pop)

# validation_data_rat <- c(
# data_c_Ntp1993RatMale230_kidney,
# data_c_Ntp1993RatMale925_liver,
# data_c_Ntp1993RatMale925_kidney,
# data_c_Ntp1993RatMale3695_liver,
# data_c_Ntp1993RatMale3695_kidney,
# data_c_Ntp1993RatMale3695_brain,
# data_c_Ntp1993Male14775_liver,
# data_c_Ntp1993Male14775_kidney,
# data_c_Ntp1993Male14775_brain,
# data_c_Ntp1993Female14775_liver,
# data_c_Ntp1993Female14775_kidney,
# data_c_Ntp1993Female14775_brain,
# data_c_Zhang25560_liver,
# data_c_Zhang25560_kidney)

# Import Rat iHg validation dataset
validation_data <- c(
  22.95 , 33.67 , 47.8 , 0.17 , 0.29 , 0.41 , 60.05 , 64.77 , 89.5 , 0.84 , 
  1.66 , 1.7 , 93.55 , 85.97 , 92.1 , 0.05 , 0.01 , 5.49 , 44.9 , 0.203 , 4.31 , 42.4 , 0.31 , 
  35.36 , 265)

str <- which(names(predRat_pop)=="CKU_1.1")
predRat_pop_data <- predRat_pop[,c(str:dim(predRat_pop)[2])]

organs <- c(
  "Kidney", "Kidney", "Kidney", 
  "Liver", "Liver", "Liver", 
  "Kidney", "Kidney", "Kidney", 
  "Liver", "Liver", "Liver", 
  "Kidney", "Kidney", "Kidney", 
  "Brain", "Brain",
  "Liver", "Kidney", "Brain",
  "Liver", "Kidney", "Brain",
  "Liver", "Kidney") 
study <- c(rep("Ntp, 1993", 23), 
           rep("Zhang 2013", 2))
dose <- c(rep(c("230 ug/kg"), 3),
          rep(c("925 ug/kg", "3695 ug/kg"), 6),
          rep("3695 ug/kg", 2), #NTP 1993 male, brain 
          rep("14775 ug/kg", 6),
          rep("25560 ug/kg", 2))

# check data length
length(na.omit(validation_data))==dim(predRat_pop_data)[2]

reshape2::melt(predRat_pop_data) |> mutate(obs = rep(na.omit(validation_data), each=sample_no)) |>
  mutate(ratio = value / obs) |> select(variable, ratio) |> 
  mutate(organs = rep(organs, each=sample_no), study= rep(study, each=sample_no), dose=rep(dose, each=sample_no)) |> 
  unite("label", study:dose) |>
  mutate(acceptance = ifelse(ratio > 3, 0, ifelse(ratio < 0.33, 0, 1))) |>
  group_by(label) |>
  mutate(accept_rate = sum(acceptance)/n()) |>
  ggplot() + 
  geom_point(aes(x=ratio, y=label)) +
  scale_x_log10() +
  facet_grid(~organs) +
  geom_vline(xintercept = 1, lty=2) +
  labs(x="predRatiction / Observation", y="Study") +
  theme(strip.background = element_blank())

validation_Rat <- reshape2::melt(predRat_pop_data) |> mutate(obs = rep(na.omit(validation_data), each=sample_no)) |>
  mutate(ratio = value / obs) |> select(variable, ratio) |> 
  mutate(organs = rep(organs, each=sample_no), study= rep(study, each=sample_no), dose=rep(dose, each=sample_no)) |> 
  unite("label", study:dose) |>
  mutate(acceptance = ifelse(ratio > 3, 0, ifelse(ratio < 0.33, 0, 1))) |>
  mutate(species = "Rat") |> 
  group_by(species, organs) |>
  mutate(accept_rate = sum(acceptance)/n()) 
validation_Rat |> distinct(species, organs, accept_rate) 
save(validation_Rat, file = "validation_Rat.RData")


