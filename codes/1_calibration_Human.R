
# Set working directory (need to be customized) -------------------------------------
if(!dir.exists("outputs")) dir.create("outputs")
if(!dir.exists("plots")) dir.create("plots")
if(!dir.exists("plots/suppl")) dir.create("plots/suppl")

# set file path ----------------------------------------------
source("MCSim/function.R")
set_PATH()
##############################################################

# package
library(bayesplot, quietly=T)
library(bayestestR, quietly=T)
library(coda, quietly=T)
library(cowplot, quietly=T)
library(data.table, quietly=T)
library(doParallel, quietly=T)
library(dplyr, quietly=T)
library(foreach, quietly=T)
library(GGally, quietly=T)
library(ggplot2, quietly=T)
library(ggpubr, quietly=T)
library(kableExtra, quietly=T)
library(LaplacesDemon, quietly=T)
library(magrittr, quietly=T)
library(PKNCA, quietly=T)
library(pksensi, quietly=T)
library(PowerTOST, quietly=T)
library(purrr, quietly=T)
library(rstan, quietly=T)
library(scales, quietly=T)
library(sensitivity, quietly=T)
library(tidyverse, quietly=T)
source("MCSim/function.R")

# posterior check
out <- c("outputs/iHgHuman_3365.out",
         "outputs/iHgHuman_4880.out",
         "outputs/iHgHuman_5916.out",
         "outputs/iHgHuman_6734.out")
data <- out |> map(fread) |> map(as.data.frame)
n_chains <- length(data)
sample_number <- dim(data[[1]])[1]
dim <- c(sample_number, n_chains, dim(data[[1]])[2])
n_iter <- dim(data[[1]])[1]
n_param <- dim(data[[1]])[2]
Human_x <- array(sample_number:(n_iter * n_chains * n_param), dim = dim)
for (i in 1:n_chains) {
  Human_x[, i, ] <- as.matrix(data[[i]][1:n_iter, ])
}
dimnames(Human_x)[[3]] <- names(data[[1]])
dim(Human_x)

# Save to RData
Human_mcmc <- Human_x[seq(50001, 100000, 1), , ]
save(Human_mcmc, file = "outputs/iHg_Human_mcmc.RData")

# data manipulate (random sample 125 iterations from 4 chains)
load("outputs/iHg_Human_mcmc.RData")
no_sample <- 125
set.seed(12345)
sample_iters <- sample(seq_len(dim(Human_mcmc)[1]), no_sample)
sample_Human_mcmc <- Human_mcmc[sample_iters, , ]
nd2 <- dim(sample_Human_mcmc)[3]
dim(sample_Human_mcmc) <- c(4 * no_sample, nd2)
dim(sample_Human_mcmc)

# posterior predictive simulation
model <- "iHgHumanBW.model"
if (!file.exists("mcsim.iHgHumanBW.model.exe")) {
  RMCSim::makemcsim(model, dir = "MCSim")
}
for (iter in seq(dim(sample_Human_mcmc)[1])){
  head(sample_Human_mcmc, iter) |> tail(1) |>
  write.table(file = "MCMC.Humancheck.dat", row.names = FALSE, sep = "\t")
  vld <- "./mcsim.iHgHumanBW.model.exe MCSim/iHgHuman.MCMC.check.in"
  system(vld)
  out <- read.delim("MCMC.Humancheck.out")
  out$iter <- iter
  if (iter == 1) Human_xx <- out
  else Human_xx <- rbind(Human_xx, out)
}
Human_xx$Output_Var |> unique()

# output manipulate
Human_xx <- Human_xx |>
  mutate(conc = ifelse(Output_Var == "Aurine", "Urine",
                          ifelse(Output_Var == "ABld", "BloodAmt",
                                 ifelse(Output_Var == "CKU", "Kidney",
                                        ifelse(Output_Var == "CBrnU", "Brain",
                                               ifelse(Output_Var == "CBldU", "Blood", ifelse(Output_Var == "CLU", "Liver", "Feces")))))))
Human_xx <- Human_xx |>
  mutate(tag = ifelse(Simulation == 1, "IV: 0.025 ug Hg/kg (Hall, 1995)","Oral: 0.09375 ug Hg/kg (Rahola, 1973)"))

Human_xx$Data[Human_xx$Data == -1] <- NA
adj_level <- Human_xx$tag |> unique()
Human_xx$tag <- factor(Human_xx$tag, level = adj_level)
Human_xx |> tail()

# define plotting element
set_theme <- theme(
  axis.text.y      = element_text(color = "black"),
  axis.ticks.y     = element_line(color = "black"),
  axis.text.x      = element_text(color = "black"),
  axis.ticks.x     = element_line(color = "black"),
  axis.line.x      = element_line(color = "black"),
  axis.line.y      = element_line(color = "black"),
  legend.key       = element_blank(),
  axis.title       = element_blank(),
  panel.background = element_blank()
)
options(warn=-1)

p1median <- Human_xx |>
  filter(Simulation == 1  & Time > 0) %>% group_by(Output_Var, Time, conc, tag)%>% 
  summarise(Data = median(Data), median = median(Prediction),
            UCL = quantile(Prediction, 0.95), LCL = quantile(Prediction, 0.05)) 
p1 <- Human_xx |> filter(Simulation == 1 & Time > 0) |>
  ggplot() +
  scale_y_log10(lim = c(10^-4, 10^1),
                breaks = trans_breaks("log10", function(x) 10^x, n = 4),
                labels = trans_format("log10", scales::math_format(10^.x))) +
  geom_line(aes(x = Time, y = Prediction, group = iter), color = "grey") +
  geom_line(data = p1median, aes(Time, median), color = "black") +
  geom_point(aes(x = Time, y = Data)) +
  facet_grid(conc ~ tag, scales = "free") +
  theme_bw() +
  set_theme


p2median <- Human_xx |>
  filter(Simulation == 2  & Time > 0) %>% group_by(Output_Var, Time, conc, tag)%>% 
  summarise(Data = median(Data), median = median(Prediction),
            UCL = quantile(Prediction, 0.95), LCL = quantile(Prediction, 0.05)) 
p2 <- Human_xx |> filter(Simulation == 2 & Time > 0) |>
  ggplot() +
  scale_y_log10(lim = c(10^-3, 10^2),
                breaks = trans_breaks("log10", function(x) 10^x, n = 2),
                labels = trans_format("log10", scales::math_format(10^.x))) +
  geom_line(aes(x = Time, y = Prediction, group = iter), color = "grey") +
  geom_line(data = p2median, aes(Time, median), color = "black") +
  geom_point(aes(x = Time, y = Data)) +
  facet_grid(conc ~ tag, scales = "free") +
  theme_bw() +
  set_theme

# add the title and axis label
title <- ggdraw() +
  draw_label(
    "Human",
    fontface = "bold",
    x = 0,
    size = 18,
    hjust = 0
  ) +
  theme(
    plot.margin = margin(0, 0, 0, 1)
  )
xlab <- ggdraw() +
  draw_label(
    "Time (hr)",
    fontface = "bold", size = 14, hjust = 0,
  ) + theme(
    plot.margin = margin(0, 0, 0, 1)
  )
ylab <- ggdraw() +
  draw_label(
    "Hg in blood, or excretion in urine or feces (ug)",
    fontface = "bold", size = 14, vjust = 0, angle = 90
  ) + theme(
    plot.margin = margin(0, 0, 0, 1)
  )

# plot
pdf(file = "plots/Figure_3C_calibration_Human.pdf", height = 6, width = 18)
plot_grid(
  ylab,
  plot_grid(
    title,
    plot_grid(p1, p2, nrow = 1, labels = c("I", "II"), rel_widths = c(0.5, 0.5)),
    xlab, nrow = 3, rel_heights = c(0.05, 1, 0.05)
  ),
  nrow = 1, rel_widths = c(0.02, 1)
)
dev.off()


