# Load packages ----------------------------------------------------------------

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
source("MCSim/function.R")


# Validation -------------------------------------------------------------------
Micesim1.1 <- fread("outputs/iHgMice_3365.out") 
Micesim2.1 <- fread("outputs/iHgMice_6734.out") 
Micesim3.1 <- fread("outputs/iHgMice_4880.out") 
Micesim4.1 <- fread("outputs/iHgMice_5916.out") 

Mice_x <-mcmc_array(list(Micesim1.1, Micesim2.1, Micesim3.1, Micesim4.1))
pars_name <- dimnames(Mice_x)[[3]]

str <- which(pars_name == "Ve_Aurine(1)")
end <- which(pars_name == "lnKurineC(1.3)")
M_pars <- pars_name[str:end]
Parameter <- c(M_pars) %>%
  strsplit(split =  "\\(1\\)") %>% unlist()

j0 <- 50001:100000
j <- j0[seq(1, length(j0), 10)]
sum_chains <- length(j)*4
mcmc_out <- Mice_x[j,,M_pars] %>% matrix(nrow = sum_chains) %>% `colnames<-`(Parameter[1:49]) 

l <- dim(mcmc_out)[1] # * dim(mcmc_out)[2]
dim(mcmc_out) <- c(l, length(M_pars))
colnames(mcmc_out) <- M_pars

sample_no <- 100
set.seed(111)
i <- sample(dim(mcmc_out)[1], sample_no) # random select 100 draws
tmp.Mice_x <- mcmc_out[i, M_pars] 
tmp.Mice_x |> write.table(file="poppredMice.dat", row.names=T, sep="\t")

vld <- "./mcsim.iHgMiceBW.model.exe modeling/iHgMice_validation.in"
system(vld)

df_Mice <- fread("poppredMice.out") %>% as.data.frame()
dim(df_Mice)

predMice_pop <- fread("poppredMice.out") |> as.data.frame()
dim(predMice_pop)
names(predMice_pop)

# Reserved command for dataset input/QA-QC-check
# validation_data <- c(
  # data_c_Ntp1993Female1250kidney, 
  # data_c_Ntp1993Female1250liver,
  # data_c_Ntp1993Female5000kidney,
  # data_c_Ntp1993Female5000liver,
  # data_c_Ntp1993Female20000kidney,
  # data_c_Ntp1993Female20000liver,
  # data_c_Ntp1993Female20000brain,
  # data_c_Ntp1993Male29550liver,
  # data_c_Ntp1993Male29550kidney,
  # data_c_Ntp1993Male29550brain,
  # data_c_Ntp1993Female29550liver,
  # data_c_Ntp1993Female29550kidney,
  # data_c_Ntp1993Female29550brain,
  # data_c_Sin6000liver,
  # data_c_Sin6000kidney,
  # data_c_Sin6000brain,
  # data_c_Wang10000liver,
  # data_c_Wang10000kidney,
  # data_c_Wang10000brain,
  # data_c_Wang10000blood
# )

# Import Mice iHg validation dataset
validation_data <- c(
     7.290 , 8.320 ,    9.860 ,  0.930 ,  0.840 ,   1.040 , 26.640 , 23.490 , 40.160 , 2.740,
     3.820 , 3.300 ,   87.940 , 97.190 , 87.960 ,   8.700 , 13.460 , 13.320 ,  0.300 , 0.490,
     0.810 , 34.350 , 170.892 ,  0.470 , 29.180 , 115.592 ,  0.493 ,  8.400 , 57.820 , 0.430,
     5.560 , 23.578 ,   0.189 ,  0.400)

str <- which(names(predMice_pop)=="CKU_1.1")
predMice_pop_data <- predMice_pop[,c(str:dim(predMice_pop)[2])]

organs <- c(
  "Kidney", "Kidney", "Kidney", 
  "Liver", "Liver", "Liver", 
  "Kidney", "Kidney", "Kidney", 
  "Liver", "Liver", "Liver", 
  "Kidney", "Kidney", "Kidney", 
  "Liver", "Liver", "Liver", 
  "Brain", "Brain", "Brain",
  "Liver", "Kidney", "Brain", # 29550 ug/kg, male, NTP
  "Liver", "Kidney", "Brain", # 29550 ug/kg, female, NTP
  "Liver", "Kidney", "Brain", # Sin
  "Liver", "Kidney", "Brain", "Blood") # Wang
study <- c(rep("Ntp 1993", 27), 
           rep("Sin 1990", 3),
           rep("Wang 2013", 4))
dose <- c(rep(c("925 ug/kg", "3695 ug/kg", "14775 ug/kg"), 6),
          rep("14775 ug/kg", 3),
          rep("29550 ug/kg", 6),
          rep("6000 ug/kg", 3),
          rep("25560 ug/kg", 4))


# check data length
length(na.omit(validation_data))==dim(predMice_pop_data)[2]

reshape2::melt(predMice_pop_data) |> mutate(obs = rep(na.omit(validation_data), each=sample_no)) |>
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
  labs(x="predMiceiction / Observation", y="Study") +
  theme(strip.background = element_blank())

validation_Mice <- reshape2::melt(predMice_pop_data) |> mutate(obs = rep(na.omit(validation_data), each=sample_no)) |>
  mutate(ratio = value / obs) |> select(variable, ratio) |> 
  mutate(organs = rep(organs, each=sample_no), study= rep(study, each=sample_no), dose=rep(dose, each=sample_no)) |> 
  unite("label", study:dose) |>
  mutate(acceptance = ifelse(ratio > 3, 0, ifelse(ratio < 0.33, 0, 1))) |>
  mutate(species = "Mouse") |> 
  group_by(species, organs) |>
  mutate(accept_rate = sum(acceptance)/n()) 
validation_Mice |> distinct(species, organs, accept_rate) 
save(validation_Mice, file = "validation_Mice.RData")


