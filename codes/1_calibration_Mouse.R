
# Set working directory (need to be customized) -------------------------------------
if(!dir.exists("outputs")) dir.create("outputs")
if(!dir.exists("plots")) dir.create("plots")
if(!dir.exists("plots/suppl")) dir.create("plots/suppl")

# set file path ----------------------------------------------
source("MCSim/function.R")
set_PATH()
##############################################################

# package
library(bayesplot, quietly =T)
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
out <- c("outputs/iHgMouse_3365.out",
         "outputs/iHgMouse_4880.out",
         "outputs/iHgMouse_5916.out",
         "outputs/iHgMouse_6734.out")
data <- out |> map(fread) |> map(as.data.frame)
n_chains <- length(data)
sample_number <- dim(data[[1]])[1]
dim <- c(sample_number, n_chains, dim(data[[1]])[2])
n_iter <- dim(data[[1]])[1]
n_param <- dim(data[[1]])[2]
Mouse_x <- array(sample_number:(n_iter * n_chains * n_param), dim = dim)
for (i in 1:n_chains) {
  Mouse_x[, i, ] <- as.matrix(data[[i]][1:n_iter, ])
}
dimnames(Mouse_x)[[3]] <- names(data[[1]])
dim(Mouse_x)

# Save to RData
Mouse_mcmc <- Mouse_x[seq(50001, 100000, 1), , ]
save(Mouse_mcmc, file = "outputs/iHg_Mouse_mcmc.RData")

# data manipulate (random sample 125 iterations from 4 chains)
load("outputs/iHg_Mouse_mcmc.RData")
no_sample <- 125
set.seed(12345)
sample_iters <- sample(seq_len(dim(Mouse_mcmc)[1]), no_sample)
sample_Mouse_mcmc <- Mouse_mcmc[sample_iters, , ]
nd2 <- dim(sample_Mouse_mcmc)[3]
dim(sample_Mouse_mcmc) <- c(4 * no_sample, nd2)
dim(sample_Mouse_mcmc)

# posterior predictive simulation
model <- "iHgMouseBW.model"
if (!file.exists("mcsim.iHgMouseBW.model.exe")) {
  RMCSim::makemcsim(model, dir = "MCSim")
}
for (iter in seq(dim(sample_Mouse_mcmc)[1])){
  head(sample_Mouse_mcmc, iter) |> tail(1) |>
  write.table(file = "MCMC.Mousecheck.dat", row.names = FALSE, sep = "\t")
  vld <- "./mcsim.iHgMouseBW.model.exe MCSim/iHgMouse.MCMC.check.in"
  system(vld)
  out <- read.delim("MCMC.Mousecheck.out")
  out$iter <- iter
  if (iter == 1) Mouse_xx <- out
  else Mouse_xx <- rbind(Mouse_xx, out)
}
Mouse_xx$Output_Var |> unique()



# output manipulate
Mouse_xx <- Mouse_xx |>
  mutate(conc = ifelse(Output_Var == "Aurine", "Urine",
                       ifelse(Output_Var == "CKU", "Kidney",
                              ifelse(Output_Var == "CBrnU", "Brain",
                                     ifelse(Output_Var == "CBldU", "Blood", "Liver")))))
Mouse_xx <- Mouse_xx |>
  mutate(tag = ifelse(Simulation == 1, "IV: 400 ug Hg/kg, (Aaseth, 1982)",
                        ifelse(Simulation == 2, "Oral water: 1,000 ug Hg/kg (Nielsen, 1992)",
                               ifelse(Simulation == 3, "Oral gavage: 3,000 ug Hg/kg/d (Dieter, 1983)",
                                      ifelse(Simulation == 4, "Oral gavage: 15,000 ug Hg/kg/d (Dieter, 1983)",
                                             ifelse(Simulation == 5, "Oral gavage: 75,000 ug Hg/kg/d (Dieter, 1983)",
                                                    ifelse(Simulation == 6, "Oral gavage: 925 ug Hg/kg/d (Ntp, 1993)",  
                                                           ifelse(Simulation == 7, "Oral gavage: 3,695 ug Hg/kg/d (Ntp, 1993)",
                                                               "Oral gavage: 14,775 ug/kg/d (Ntp, 1993)"))))))))

Mouse_xx$Data[Mouse_xx$Data == -1] <- NA
adj_level <- Mouse_xx$tag |> unique()
Mouse_xx$tag <- factor(Mouse_xx$tag, level = adj_level)
Mouse_xx |> tail()

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

Mouse_xx$tag = factor(Mouse_xx$tag, levels=c('IV: 400 ug Hg/kg, (Aaseth, 1982)','Oral water: 1,000 ug Hg/kg (Nielsen, 1992)','Oral gavage: 3,000 ug Hg/kg/d (Dieter, 1983)'
                                         ,'Oral gavage: 15,000 ug Hg/kg/d (Dieter, 1983)' ,'Oral gavage: 75,000 ug Hg/kg/d (Dieter, 1983)' ,'Oral gavage: 925 ug Hg/kg/d (Ntp, 1993)' ,'Oral gavage: 3,695 ug Hg/kg/d (Ntp, 1993)','Oral gavage: 14,775 ug/kg/d (Ntp, 1993)'))
plmedian <- Mouse_xx |>
  filter(Simulation == 1 & Time > 0) %>% group_by(Output_Var, Time, conc, tag)  %>% 
  summarise(Data = median(Data), median = median(Prediction),
            UCL = quantile(Prediction, 0.95), LCL = quantile(Prediction, 0.05))  
p1 <- Mouse_xx |>
  filter(Simulation == 1 & Time > 0) |>
  ggplot() +
  scale_y_log10(lim = c(10^-3, 10^2),
                breaks = trans_breaks("log10", function(x) 10^x, n = 3),
                labels = trans_format("log10", scales::math_format(10^.x))) +
  geom_line(aes(x = Time, y = Prediction, group = iter), color = "grey") +
  geom_line(data = plmedian, aes(Time, median), color = "black") +
  geom_point(aes(x = Time, y = Data)) +
  facet_grid(conc ~ tag, scales = "free") +
  theme_bw() +
  set_theme


p2median <- Mouse_xx |>
  filter(Simulation == 2 & Time > 0 & conc == "Blood") %>% group_by(Output_Var, Time, conc, tag)%>% 
  summarise(Data = median(Data), median = median(Prediction),
            UCL = quantile(Prediction, 0.95), LCL = quantile(Prediction, 0.05)) 
p2 <- Mouse_xx |>
  filter(Simulation == 2 & Time > 0) |>
  ggplot() +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x, n = 3),
                labels = trans_format("log10", scales::math_format(10^.x))) +
  geom_line(aes(x = Time, y = Prediction, group = iter), color = "grey") +
  geom_line(data = p2median, aes(Time, median), color = "black") +
  geom_point(aes(x = Time, y = Data)) +
  facet_grid(conc ~ tag, scales = "free") +
  theme_bw() +
  set_theme


p3median <- Mouse_xx |>
  filter(Simulation %in% c(3:5) & Time > 0) %>% group_by(Output_Var, Time, conc, tag)%>% 
  summarise(Data = median(Data), median = median(Prediction),
            UCL = quantile(Prediction, 0.95), LCL = quantile(Prediction, 0.05)) 
p3 <- Mouse_xx |>
  filter(Simulation %in% c(3:5) & Time > 0) |>
  ggplot() +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", scales::math_format(10^.x))) +
  geom_line(aes(x = Time, y = Prediction, group = iter), color = "grey") +
  geom_line(data = p3median, aes(Time, median), color = "black") +
  geom_point(aes(x = Time, y = Data)) +
  facet_grid(conc ~ tag, scales = "free") +
  theme_bw() +
  set_theme

f_label4 <- data.frame(tag = c("Oral gavage: 925 ug Hg/kg/d (Ntp, 1993)"), label = c("data not available"), conc = c("Brain"))
f_label5 <- data.frame(tag = c("Oral gavage: 3,695 ug Hg/kg/d (Ntp, 1993)"), label = c("data not available"), conc = c("Brain"))
p4median <- Mouse_xx |>
  filter(Simulation %in% c(6:8) & Time > 0) %>% group_by(Output_Var, Time, conc, tag)%>% 
  summarise(Data = median(Data), median = median(Prediction),
            UCL = quantile(Prediction, 0.95), LCL = quantile(Prediction, 0.05)) 
p4 <- Mouse_xx |>
  filter(Simulation %in% c(6:8) & Time > 0) |>
  ggplot() +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x, n = 5),
                labels = trans_format("log10", scales::math_format(10^.x))) +
  geom_line(aes(x = Time, y = Prediction, group = iter), color = "grey") +
  geom_line(data = p4median, aes(Time, median), color = "black") +
  geom_point(aes(x = Time, y = Data)) +
  facet_grid(conc ~ factor(tag), scales = "free") +
  geom_text(aes(x = 2200, y = 0.001, label = label), data = f_label4, size=8, color = "grey", position = position_stack(vjust = .5)) +
  geom_text(aes(x = 2200, y = 0.001, label = label), data = f_label5, size=8, color = "grey", position = position_stack(vjust = .5)) +
  theme_bw() +
  set_theme


# add the title and axis label
title <- ggdraw() +
  draw_label(
    "Mouse",
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
    "Tisse Hg concentration (ug/mL) / Hg excreted in urine or feces (ug)",
    fontface = "bold", size = 14, vjust = 0, angle = 90
  ) + theme(
    plot.margin = margin(0, 0, 0, 1)
  )

# plot
pdf(file = "plots/Figure_3A_calibration_Mouse.pdf", height = 12, width = 20)
plot_grid(
  ylab,
  plot_grid(
    title,
    plot_grid(
      plot_grid(p1, p2, nrow = 2, labels = c("I", "II"),
                rel_heights = c(3 / 4, 1 / 4)),
      plot_grid(
        p3, p4, nrow = 2,
        labels = c("III", "IV")
      ),
      nrow = 1, rel_widths = c(0.33, 0.66)
    ),
    xlab, nrow = 3, rel_heights = c(0.05, 1, 0.05)),
  nrow = 1, rel_widths = c(0.02, 1)
)
dev.off()

