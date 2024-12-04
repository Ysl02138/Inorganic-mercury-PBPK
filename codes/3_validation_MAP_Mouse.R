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
Mousesim1.1 <- fread("outputs/iHgMouse_3365.out") 
Mousesim2.1 <- fread("outputs/iHgMouse_6734.out") 
Mousesim3.1 <- fread("outputs/iHgMouse_4880.out") 
Mousesim4.1 <- fread("outputs/iHgMouse_5916.out") 

Mouse_x <-mcmc_array(list(Mousesim1.1, Mousesim2.1, Mousesim3.1, Mousesim4.1))
pars_name <- dimnames(Mouse_x)[[3]]

str <- which(pars_name == "Ve_Aurine(1)")
end <- which(pars_name == "lnKbrnC(1.4)")

parms <- pars_name[str:end]

j <- seq(50001, 100000, 1) 
sum_chains <- length(j)*4

sample_no = 500 # 500 virtual study
d <- Mouse_x[j,,] %>% matrix(nrow = sum_chains) %>% as.data.frame() %>% 
  `colnames<-`(dimnames(Mouse_x)[[3]])
i <- which(d[,"LnPosterior"]==max(d[,"LnPosterior"])) %>%
  rep(sample_no) # Pick the sample based on the maximum a posteriori estimation

tmp.Mouse_x <- d[i, parms] 
tmp.Mouse_x |> write.table(file="poppredMouse.dat", row.names=T, sep="\t")

vld <- "./mcsim.iHgMouseBW.model.exe MCSim/iHgMouse_validation.in"
system(vld)

df_Mouse <- fread("poppredMouse.out") %>% as.data.frame()
dim(df_Mouse)

predMouse_pop <- fread("poppredMouse.out") |> as.data.frame()
dim(predMouse_pop)
names(predMouse_pop)

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

# Import Mouse iHg validation dataset
validation_data <- c(
  7.290 , 8.320 ,    9.860 ,  0.930 ,  0.840 ,   1.040 , 26.640 , 23.490 , 40.160 , 2.740,
  3.820 , 3.300 ,   87.940 , 97.190 , 87.960 ,   8.700 , 13.460 , 13.320 ,  0.300 , 0.490,
  0.810 , 34.350 , 170.892 ,  0.470 , 29.180 , 115.592 ,  0.493 ,  8.400 , 57.820 , 0.430,
  5.560 , 23.578 ,   0.189 ,  0.400)

str <- which(names(predMouse_pop)=="CKU_1.1")
predMouse_pop_data <- predMouse_pop[,c(str:dim(predMouse_pop)[2])]

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
length(na.omit(validation_data))==dim(predMouse_pop_data)[2]

reshape2::melt(predMouse_pop_data) |> mutate(obs = rep(na.omit(validation_data), each=sample_no)) |>
  mutate(Mouseio = value / obs) |> select(variable, Mouseio) |> 
  mutate(organs = rep(organs, each=sample_no), study= rep(study, each=sample_no), dose=rep(dose, each=sample_no)) |> 
  unite("label", study:dose) |>
  mutate(acceptance = ifelse(Mouseio > 3, 0, ifelse(Mouseio < 0.33, 0, 1))) |>
  group_by(label) |>
  mutate(accept_Mousee = sum(acceptance)/n()) |>
  ggplot() + 
  geom_point(aes(x=Mouseio, y=label)) +
  scale_x_log10() +
  facet_grid(~organs) +
  geom_vline(xintercept = 1, lty=2) +
  labs(x="predratiction / Observation", y="Study") +
  theme(strip.background = element_blank())

validation_Mouse <- reshape2::melt(predMouse_pop_data) |> mutate(obs = rep(na.omit(validation_data), each=sample_no)) |>
  mutate(ratio = value / obs) |> select(variable, ratio) |> 
  mutate(organs = rep(organs, each=sample_no), study= rep(study, each=sample_no), dose=rep(dose, each=sample_no)) |> 
  unite("label", study:dose) |>
  mutate(acceptance = ifelse(ratio > 3, 0, ifelse(ratio < 0.33, 0, 1))) |>
  mutate(species = "Mouse") |> 
  group_by(species, organs) |>
  mutate(accept_rate = sum(acceptance)/n()) 
validation_Mouse |> distinct(species, organs, accept_rate) 
save(validation_Mouse, file = "validation_Mouse.RData")


