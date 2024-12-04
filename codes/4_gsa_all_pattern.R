# Tidy
rm(list=ls())

# Load packages ----------------------------------------------------------------
pacman::p_load(pksensi, ggplot2, ggpubr,cowplot, reshape2, dplyr, viridis)

# Create mod.exe ---------------------------------------------------------------
source("MCSim/function.R")
set_PATH()
# makemod()

# Compile model code for Mouse-----------------------------------------------------------
mName <- "iHgMouseBW.model" 

# Parameter range
load("outputs/iHg_Mouse_mcmc.RData")
out_vars <- dimnames(Mouse_mcmc)[[3]] 
nd1 <- dim(Mouse_mcmc)[1] * 4
nd2 <- dim(Mouse_mcmc)[3]
dim(Mouse_mcmc) <- c(nd1,nd2)
colnames(Mouse_mcmc) <- out_vars
str <- which(colnames(Mouse_mcmc) == "M_lnPLC(1)")
end <- which(colnames(Mouse_mcmc) == "M_lnKbrnC(1)")
pop_min <- Mouse_mcmc[,str:end] |> apply(2, min)
pop_max <- Mouse_mcmc[,str:end] |> apply(2, max)
pop_max - pop_min


# Parameter space generating
params <- c("lnPLC", "lnPKC", "lnPBrnC", "lnPRestC", 
            "lnKabsC", "lnKunabsC", "lnKbileC", "lnKurineC", "lnKbrnC")
q <- rep("qunif", 9)
q.arg <- list(list(pop_min[1], pop_max[1]), 
              list(pop_min[2], pop_max[2]), 
              list(pop_min[3], pop_max[3]), 
              list(pop_min[4], pop_max[4]), 
              list(pop_min[5], pop_max[5]), 
              list(pop_min[6], pop_max[6]), 
              list(pop_min[7], pop_max[7]), 
              list(pop_min[8], pop_max[8]), 
              list(pop_min[9], pop_max[9]) 
)
set.seed(1234)
x <- rfast99(params = params, n = 10000, q = q, q.arg = q.arg, replicate = 1)

# Single dose simulation
conditions <- c("BW0 = 64;", 
                "BWgrowth = 0;", 
                "Growthrate = 0;", 
                "sex = 1;", 
                "TChng = 0.5;", 
                "PDose = PerDose(0.0, 24,  0, 0.05);",
                "IVDose = PerDose(0.0 , 24,  0, 0.003);",
                "Drink = 330.0",
                "expowk =  1;", 
                "expodur = 1;", 
                "Drink = 330.0;")
vars <- c("AUCCBld", "AUCCL", "AUCCK", "AUCCBrn")
times <- c(4368)
out <- solve_mcsim(x, mName = mName,
                   params = params, 
                   time = times, 
                   vars = vars,
                   condition = conditions, 
                   rtol = 1e-7, atol = 1e-9)
mSI <- out$mSI
iSI <- out$iSI
dim(mSI) <- dim(iSI) <- c(9,4)
colnames(mSI) <- colnames(iSI) <- vars
rownames(mSI) <- rownames(iSI) <- params
x1 <- melt(mSI) |> `colnames<-`(c("parameter", "variable", "index")) |>
  mutate(order = "main effect", specice = "Mouse", exposure = "Mouse")
x2 <- melt(iSI) |> `colnames<-`(c("parameter", "variable", "index")) |>
  mutate(order = "interaction", specice = "Mouse", exposure = "Mouse")


# Compile model code for rat -----------------------------------------------------------
mName <- "iHgRatBW.model" 

# Parameter range
load("outputs/iHg_Rat_mcmc.RData")
out_vars <- dimnames(Rat_mcmc)[[3]] 
nd1 <- dim(Rat_mcmc)[1] * 4
nd2 <- dim(Rat_mcmc)[3]
dim(Rat_mcmc) <- c(nd1,nd2)
colnames(Rat_mcmc) <- out_vars
str <- which(colnames(Rat_mcmc) == "M_lnPLC(1)")
end <- which(colnames(Rat_mcmc) == "M_lnKbrnC(1)")
pop_min <- Rat_mcmc[,str:end] |> apply(2, min)
pop_max <- Rat_mcmc[,str:end] |> apply(2, max)
pop_max - pop_min


# Parameter space geneRating
params <- c("lnPLC", "lnPKC", "lnPBrnC", "lnPRestC", 
            "lnKabsC", "lnKunabsC", "lnKbileC", "lnKurineC", "lnKbrnC")
q <- rep("qunif", 9)
q.arg <- list(list(pop_min[1], pop_max[1]), 
              list(pop_min[2], pop_max[2]), 
              list(pop_min[3], pop_max[3]), 
              list(pop_min[4], pop_max[4]), 
              list(pop_min[5], pop_max[5]), 
              list(pop_min[6], pop_max[6]), 
              list(pop_min[7], pop_max[7]), 
              list(pop_min[8], pop_max[8]), 
              list(pop_min[9], pop_max[9]) 
)
set.seed(1234)
x <- rfast99(params = params, n = 10000, q = q, q.arg = q.arg, replicate = 1)

# Single dose simulation
conditions <- c("BW0 = 64;", 
                "BWgrowth = 0;", 
                "Growthrate = 0;", 
                "sex = 1;", 
                "TChng = 0.5;", 
                "PDose = PerDose(0.0, 24,  0, 0.05);",
                "IVDose = PerDose(0.0 , 24,  0, 0.003);",
                "Drink = 330.0",
                "expowk =  1;", 
                "expodur = 1;", 
                "Drink = 330.0;")
vars <- c("AUCCBld", "AUCCL", "AUCCK", "AUCCBrn")
times <- c(4368)
out <- solve_mcsim(x, mName = mName,
                   params = params, 
                   time = times, 
                   vars = vars,
                   condition = conditions, 
                   rtol = 1e-7, atol = 1e-9)
mSI <- out$mSI
iSI <- out$iSI
dim(mSI) <- dim(iSI) <- c(9,4)
colnames(mSI) <- colnames(iSI) <- vars
rownames(mSI) <- rownames(iSI) <- params
x3 <- melt(mSI) |> `colnames<-`(c("parameter", "variable", "index")) |>
  mutate(order = "main effect", specice = "rat", exposure = "Rat")
x4 <- melt(iSI) |> `colnames<-`(c("parameter", "variable", "index")) |>
  mutate(order = "interaction", specice = "rat", exposure = "Rat")


# Compile model code for human -----------------------------------------------------------
mName <- "iHgHumanBW.model" 

# Parameter range
load("outputs/iHg_Human_mcmc.RData")
out_vars <- dimnames(Human_mcmc)[[3]] 
nd1 <- dim(Human_mcmc)[1] * 4
nd2 <- dim(Human_mcmc)[3]
dim(Human_mcmc) <- c(nd1,nd2)
colnames(Human_mcmc) <- out_vars
str <- which(colnames(Human_mcmc) == "M_lnPLC(1)")
end <- which(colnames(Human_mcmc) == "M_lnKbrnC(1)")
pop_min <- Human_mcmc[,str:end] |> apply(2, min)
pop_max <- Human_mcmc[,str:end] |> apply(2, max)
pop_max - pop_min


# Parameter space generating
params <- c("lnPLC", "lnPKC", "lnPBrnC", "lnPRestC", 
  "lnKabsC", "lnKunabsC", "lnKbileC", "lnKurineC", "lnKbrnC")
q <- rep("qunif", 9)
q.arg <- list(list(pop_min[1], pop_max[1]), 
  list(pop_min[2], pop_max[2]), 
  list(pop_min[3], pop_max[3]), 
  list(pop_min[4], pop_max[4]), 
  list(pop_min[5], pop_max[5]), 
  list(pop_min[6], pop_max[6]), 
  list(pop_min[7], pop_max[7]), 
  list(pop_min[8], pop_max[8]), 
  list(pop_min[9], pop_max[9]) 
)
set.seed(1234)
x <- rfast99(params = params, n = 10000, q = q, q.arg = q.arg, replicate = 1)

# Single dose simulation
conditions <- c("BW0 = 64;", 
                "BWgrowth = 0;", 
                "Growthrate = 0;", 
                "sex = 1;", 
                "TChng = 0.5;", 
                "PDose = PerDose(0.0, 24,  0, 0.05);",
                "IVDose = PerDose(0.0 , 24,  0, 0.003);",
                "Drink = 330.0",
                "expowk =  1;", 
                "expodur = 1;", 
                "Drink = 330.0;")
vars <- c("AUCCBld", "AUCCL", "AUCCK", "AUCCBrn")
times <- c(4368)
out <- solve_mcsim(x, mName = mName,
                   params = params, 
                   time = times, 
                   vars = vars,
                   condition = conditions, 
                   rtol = 1e-7, atol = 1e-9)
mSI <- out$mSI
iSI <- out$iSI
dim(mSI) <- dim(iSI) <- c(9,4)
colnames(mSI) <- colnames(iSI) <- vars
rownames(mSI) <- rownames(iSI) <- params
x5 <- melt(mSI) |> `colnames<-`(c("parameter", "variable", "index")) |>
  mutate(order = "main effect", specice = "human", exposure = "Human")
x6 <- melt(iSI) |> `colnames<-`(c("parameter", "variable", "index")) |>
  mutate(order = "interaction", specice = "human", exposure = "Human")


x_total <- do.call(rbind, list(x1, x2, x3 ,x4, x5, x6))
x_total$order <- factor(x_total$order, level=c("interaction", "main effect"))

set_theme <- theme(
  legend.position  = "top",
  legend.title     = element_blank(),
  axis.title       = element_blank(),
)

pdf("plots/Figure_5_GSA.pdf",width=11)
ggplot(x_total, aes(x = parameter, y = index, fill=order)) +
  geom_bar_pattern(position = "stack", stat = "identity", pattern_fill="white", aes(pattern=order)) + 
  coord_flip() +
  facet_grid(variable ~ factor (exposure, levels=c("Mouse", "Rat", "Human"))) + 
  theme_bw() +
  ggtitle(" ") +
  scale_fill_viridis(discrete = TRUE, direction = -1, end = 0.8) + 
  set_theme
dev.off()  
# 