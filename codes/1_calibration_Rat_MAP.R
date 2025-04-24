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
library(formattable)
library(GGally, quietly=T)
library(ggplot2, quietly=T)
library(ggpubr, quietly=T)
library(ggforce, quietly=T)
library(grid, quietly=T)
library(gridExtra, quietly=T)
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


# posterior check with MAP
Ratsim1.1 <- fread("outputs/iHgRat_3365.out") 
Ratsim2.1 <- fread("outputs/iHgRat_6734.out") 
Ratsim3.1 <- fread("outputs/iHgRat_4880.out") 
Ratsim4.1 <- fread("outputs/iHgRat_5916.out") 

Rat_x <-mcmc_array(list(Ratsim1.1, Ratsim2.1, Ratsim3.1, Ratsim4.1))

pars_name <- dimnames(Rat_x)[[3]]
str <- which(pars_name == "M_lnPLC(1)")
end <- which(pars_name == "M_lnKbrnC(1)") 
M_pars <- pars_name[str:end]

str <- which(pars_name == "V_lnPLC(1)")
end <- which(pars_name == "V_lnKbrnC(1)") 
V_pars <- pars_name[str:end]

str <- which(pars_name == "Ve_Aurine(1)")
end <- which(pars_name == "Ve_CBrnU(1)") 
Ve_pars <- pars_name[str:end]

Rat_mnt <- monitor(Rat_x[,,c(M_pars, V_pars, Ve_pars)], digit=4)

str <- which(pars_name == "Ve_Aurine(1)")
end <- which(pars_name == "V_lnKbrnC(1)")
parms <- pars_name[str:end]

j <- seq(50001, 100000, 1) 
sum_chains <- length(j)*4 #

n = 500 # 500 virtual study
d <- Rat_x[j,,] %>% matrix(nrow = sum_chains) %>% as.data.frame() %>% 
  `colnames<-`(dimnames(Rat_x)[[3]])
i <- which(d[,"LnPosterior"]==max(d[,"LnPosterior"])) %>%
  rep(n) # Pick the sample based on the maximum a posteriori estimation

tmp.Rat_x <- d[i, parms] 
tmp.Rat_x %>% write.table(file="iHg_calibration_MAP_Rat.dat", row.names=T, sep="\t")
vld <- "./mcsim.iHgRatBW.model.exe MCSim/iHg_calibration_MAP_Rat.in"
system(vld)

df <- fread("iHg_calibration_MAP_Rat.out") %>% as.data.frame()
file.remove(c("iHg_calibration_MAP_Rat.dat", "iHg_calibration_MAP_Rat.out"))

df

str_1.1 <- which(names(df) == "Aurine_1.1")
end_8.10 <- which(names(df) == "CBrnU_8.10")


#####################################################################
#####################################################################
#####################################################################

# Import Rat iHg calibration dataset --- same datasets as "iHg_calibration_MAP_Rat.in"
Time <- c(0.5, 1, 2, 4, 8, 12, 24,	48,	144,	360, 
          
          0.5, 1, 2, 4, 8, 12, 24, 48,	144,	360, 
          
          0.5, 1, 2, 4, 8, 12, 24, 48, 672 , 1344 , 2016, 
          
          0.5, 1, 2, 4, 8, 12, 24, 48, 336,	672,	1008,	1344, 
          
          0.5, 1, 2, 4, 8, 12, 24, 48, 	336,	672,	1008,	1344, 
          
          0.5, 1, 2, 4, 8, 12, 24, 48, 168 , 336 , 504 , 672 , 840 , 1008 , 1176 , 1344, 
          
          0.5, 1, 2, 4, 8, 12, 24, 48, 	336,	672,	1008,	1344,  
          
          0.5, 1, 2, 4, 8, 12, 24, 48, 336,	672,	1008,	1344, 
          
          0.5, 1, 2, 4, 8, 12, 24, 48, 168 , 336 , 504 , 672 , 840 , 1008 , 1176 , 1344, 
          
          0.5, 1, 2, 4, 8, 12, 24, 48, 336 , 672 , 1008 , 1344, 
          
          0.5, 1, 2, 4, 8, 12, 24, 48, 	336,	672,	1008,	1344, 
          
          0.5, 1, 2, 4, 8, 12, 24, 48, 168 , 336 , 504 , 672 , 840 , 1008 , 1176 , 1344,
          
          0.5, 1, 2, 4, 8, 12, 24, 48, 1428, 2868,	4308,
          
          0.5, 1, 2, 4, 8, 12, 24, 48, 1428, 2868,	4308, 
          
          0.5, 1, 2, 4, 8, 12, 24, 48, 1428, 2868,	4308,
          
          0.5, 1, 2, 4, 8, 12, 24, 48, 1428, 2868,	4308,
          
          0.5, 1, 2, 4, 8, 12, 24, 48, 1428, 2868,	4308,
          
          0.5, 1, 2, 4, 8, 12, 24, 48, 1428,  4308) # Observation time 

Conc <- c(-1, -1, -1, -1, -1, -1, 1.306962025,	2.458860759,	3.743670886,	6.091772152, 
          
          -1, -1, -1, -1, -1, -1, 3.344936709, 6.82278481, 10.65506329, 13.18037975, 
          
          -1, -1, -1, -1, -1, -1, -1, -1, 0.042074074 , 0.047703704 , 0.057777778, 
          
          -1, -1, -1, -1, -1, -1, -1, -1, 0.603448276,	3.448275862,	6.637931034,	9.827586207, 
          
          -1, -1, -1, -1, -1, -1, -1, -1, 206.8965517,	448.2758621,	724.137931,	984.615385, 
          
          
          -1, -1, -1, -1, -1, -1, -1, -1, 0.00848,	0.01066,	0.01045,	0.0107,	0.01127,	0.00811,	0.00659,	0.00525, 
          
          -1, -1, -1, -1, -1, -1, -1, -1, 8.720689655,	29.31034483,	87.93103448,	163.7931034, 
          
          -1, -1, -1, -1, -1, -1, -1, -1, 2086.956522,	5113.043478,	7513.043478,	10121.73913, 
          
          -1, -1, -1, -1, -1, -1, -1, -1, 0.07167,	0.05471,	0.06394,	0.06063,	0.04959,	0.06375,	0.06978,	0.05874, 
          
          -1, -1, -1, -1, -1, -1, -1, -1, 818.9655172 , 2068.965517 , 3405.172414 , 4827.586207, 
          
          -1, -1, -1, -1, -1, -1, -1, -1, 14102.5641,	28846.15385,	45512.82051,	60897.4359,  
          
          -1, -1, -1, -1, -1, -1, -1, -1,  0.29892 , 0.1852 , 0.21552 , 0.34776 , 0.28396 , 0.32497 , 0.2738 , 0.25756,
          
          -1, -1, -1, -1, -1, -1, -1, -1, 29.37, 55.59, 46.98, 
          
          -1, -1, -1, -1, -1, -1, -1, -1, 100.87, 97.19, 86.08,
          
          -1, -1, -1, -1, -1, -1, -1, -1, 0.17, 0.29,	0.41,
          
          -1, -1, -1, -1, -1, -1, -1, -1, 101.87, 122.19, 92.78,
          
          -1, -1, -1, -1, -1, -1, -1, -1, 1.00, 1.62, 1.86,
          
          -1, -1, -1, -1, -1, -1, -1, -1, 0.06,  0.01) # Observation concentration

formattable(Time , digits = 3, format = "f")
cal_data <- data.frame(Conc, Time)
# study <- c("Rothstein, 1960", "Oriquat, 2013", "Morcillo, 1995", "NTP, 1993")

dose <- c("IV: 250 ug Hg/kg", "Oral: 2,770 ug Hg/kg", "Oral water: 100 ug Hg/kg/d","Oral water: 1,000 ug Hg/kg/d", "Oral water: 7,200 ug Hg/kg/d", "Oral gavage: 230 ug Hg/kg/d", "Oral gavage: 925 ug Hg/kg/d", "Oral gavage: 3,695 ug/kg/d")

organs <- c(rep("Urine", 10), rep("Feces", 10), rep("Blood", 11), rep("Urine", 12), rep("Feces", 12),
            rep("Blood", 16), rep("Urine", 12), rep("Feces", 12), rep("Blood", 16),  rep("Urine", 12), rep("Feces", 12), 
            rep("Blood", 16), rep("Kidney", 11),  rep("Kidney", 11),  rep("Liver", 11),   rep("Kidney", 11),  rep("Liver", 11), rep("Brain", 10))

Output_Var <- c(rep("Aurine", 10), rep("Afeces", 10), rep("CBldU", 11), rep("Aurine", 12), rep("Afeces", 12),
                rep("CBldU", 16), rep("Aurine", 12), rep("Afeces", 12), rep("CBldU", 16),  rep("Aurine", 12), rep("Afeces", 12), 
                rep("CBldU", 16),  rep("CKU", 11), rep("CKU", 11),  rep("CLU", 11), rep("CKU", 11),  rep("CLU", 11), rep("CBrnU", 10))

study <- c(rep("(Rothstein, 1960)", 20),
           rep("(Oriquat, 2013)", 11), 
           rep("(Morcillo, 1995)", 120), 
           rep("(NTP, 1993)", 65))

dose <- c(rep("IV: 250 ug Hg/kg", 20),
          rep("Oral: 2,770 ug Hg/kg", 11),
          rep("Oral water: 100 ug Hg/kg/d", 40),
          rep("Oral water: 1000 ug Hg/kg/d", 40),
          rep("Oral water: 7720 ug Hg/kg/d", 40),
          rep("Oral gavage: 230 ug Hg/kg/d", 11),
          rep("Oral gavage: 925 ug Hg/kg/d", 22),
          rep("Oral gavage: 3695 ug Hg/kg/d", 32))

Simulation <- c(
  rep("1", 20),
  rep("2", 11),
  rep("3", 40),
  rep("4", 40),
  rep("5", 40),
  rep("6", 11),
  rep("7", 22),
  rep("8", 32))

# ========================================
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
# ========================================

cal_data$organs <- organs
cal_data$study  <- study
cal_data$dose <- dose
cal_data$Simulation <- Simulation
cal_data  %<>% unite(label, dose, study, sep = " ", remove = FALSE)
label <- c(cal_data$label)
cal_data <- select(cal_data, -dose, -study)

qt.line <- df[,c(str_1.1:end_8.10)] %>% gather() %>% 
  `colnames<-`(c("var","Prediction")) %>%
  add_column(Time =  rep(rep(Time, each = n),1)) %>%
  add_column(Conc =  rep(rep(Conc, each = n),1)) %>%
  add_column(label = rep(rep(label, each = n),1)) %>%
  add_column(Simulation =  rep(rep(Simulation, each = n),1)) %>%
  add_column(Output_Var =  rep(rep(Output_Var, each = n),1)) %>%
  add_column(organs = rep(rep(organs, each = n),1))

###################################################################################

p1_summary <- qt.line |>
  filter(Simulation %in% c(1), Time > 0) |>
  group_by(label, organs, Time, Simulation) |>
  summarise(
    med = median(Prediction),
    UCL = quantile(Prediction, 0.975),
    LCL = quantile(Prediction, 0.025),
    .groups = "drop"
  )

# Filter observed data for non-missing points
obs_data <- cal_data %>% filter(Conc > 0, Simulation == "1")

# Plot
p1 <- ggplot() +
  # 95% CI ribbon from simulation
  geom_ribbon(
    data = p1_summary,
    aes(x = Time, ymin = LCL, ymax = UCL, group = label),
    fill = "grey", alpha = 0.3
  ) +
  # Predicted median line in blue, no legend
  geom_line(
    data = p1_summary,
    aes(x = Time, y = med, group = label),
    color = "blue",
    size = 0.75
  ) +
  # Observed data points
  geom_point(
    data = obs_data,
    aes(x = Time, y = Conc),
    size = 1.5
  ) +
  facet_grid(organs ~ label, scales = "free_y") +
  scale_y_log10(lim = c(10^-1, 10^2),
                breaks = trans_breaks("log10", function(x) 10^x, n = 4),
                labels = trans_format("log10", scales::math_format(10^.x)))+
  labs(
    x = "Time (hr)",
    y = "Concentration / Prediction (log scale)"
  ) +
  theme_bw() +
  theme(
    legend.position = "none",
    strip.text = element_text(size = 8)
  ) +
  set_theme


###################################################################################

p2_summary <- qt.line |>
  filter(Simulation %in% c(2), Time > 0) |>
  group_by(label, organs, Time, Simulation) |>
  summarise(
    med = median(Prediction),
    UCL = quantile(Prediction, 0.975),
    LCL = quantile(Prediction, 0.025),
    .groups = "drop"
  )

# Filter observed data for non-missing points
obs_data <- cal_data %>% filter(Conc > 0, Simulation == "2")

# Plot
p2 <- ggplot() +
  # 95% CI ribbon from simulation
  geom_ribbon(
    data = p2_summary,
    aes(x = Time, ymin = LCL, ymax = UCL, group = label),
    fill = "grey", alpha = 0.3
  ) +
  # Predicted median line in blue, no legend
  geom_line(
    data = p2_summary,
    aes(x = Time, y = med, group = label),
    color = "blue",
    size = 0.75
  ) +
  # Observed data points
  geom_point(
    data = obs_data,
    aes(x = Time, y = Conc),
    size = 1.5
  ) +
  facet_grid(organs ~ label, scales = "free_y") +
  scale_y_log10(lim = c(10^-3, 10^0),
                breaks = trans_breaks("log10", function(x) 10^x, n = 4),
                labels = trans_format("log10", scales::math_format(10^.x))) +
  labs(
    x = "Time (hr)",
    y = "Concentration / Prediction (log scale)"
  ) +
  theme_bw() +
  theme(
    legend.position = "none",
    strip.text = element_text(size = 8)
  ) +
  set_theme


###################################################################################

p3_summary <- qt.line |>
  filter(Simulation %in% c(3:5), Time > 0) |>
  group_by(label, organs, Time, Simulation) |>
  summarise(
    med = median(Prediction),
    UCL = quantile(Prediction, 0.975),
    LCL = quantile(Prediction, 0.025),
    .groups = "drop"
  )

# Filter observed data for non-missing points
obs_data <- cal_data %>% filter(Conc > 0, Simulation %in% c(3:5))

# Plot
p3 <- ggplot() +
  # 95% CI ribbon from simulation
  geom_ribbon(
    data = p3_summary,
    aes(x = Time, ymin = LCL, ymax = UCL, group = label),
    fill = "grey40", alpha = 0.3
  ) +
  # Predicted median line in blue, no legend
  geom_line(
    data = p3_summary,
    aes(x = Time, y = med, group = label),
    color = "blue",
    size = 0.75
  ) +
  # Observed data points
  geom_point(
    data = obs_data,
    aes(x = Time, y = Conc),
    size = 1.5
  ) +
  facet_grid(organs ~ label, scales = "free_y") +
  scale_y_log10(# lim = c(10^-6, 10^3),
    breaks = trans_breaks("log10", function(x) 10^x, n = 4),
    labels = trans_format("log10", scales::math_format(10^.x))) +
  labs(
    x = "Time (hr)",
    y = "Concentration / Prediction (log scale)"
  ) +
  theme_bw() +
  theme(
    legend.position = "none",
    strip.text = element_text(size = 8)
  ) +
  set_theme

###################################################################################

# Summarize simulation data
p4_summary <- qt.line |>
  filter(Simulation %in% c(6:8), Time > 0) |>
  group_by(label, organs, Time, Simulation) |>
  summarise(
    med = median(Prediction),
    UCL = quantile(Prediction, 0.975),
    LCL = quantile(Prediction, 0.025),
    .groups = "drop"
  )

# Filter observed data for non-missing points
obs_data <- cal_data %>% 
  filter(Conc > 0, Simulation %in% c(6:8))

# Set custom factor levels to control facet order
label_order <- unique(cal_data$label)
organs_order <- c("Brain", "Kidney", "Liver")

p4_summary$label <- factor(p4_summary$label, levels = label_order)
p4_summary$organs <- factor(p4_summary$organs, levels = organs_order)

obs_data$label <- factor(obs_data$label, levels = label_order)
obs_data$organs <- factor(obs_data$organs, levels = organs_order)

# Create annotation for missing data in organ = "Brain" and Simulation 6 or 7
labels_sim_6_7 <- cal_data %>%
  filter(Simulation %in% c("6", "7")) %>%
  distinct(label) %>%
  pull(label)

labels_sim_6 <- cal_data %>%
  filter(Simulation == "6") %>%
  distinct(label) %>%
  pull(label)

annot_df <- bind_rows(
  expand.grid(label = labels_sim_6, organs = "Liver"),
  expand.grid(label = labels_sim_6_7, organs = "Brain")
) %>%
  mutate(
    x = 2200,  # position for the text (log scale)
    y = 1,
    text = "data not available"
  )

# Set factor levels to match plot
annot_df$label <- factor(annot_df$label, levels = label_order)
annot_df$organs <- factor(annot_df$organs, levels = organs_order)

# Plot
p4 <- ggplot() +
  # 95% CI ribbon from simulation
  geom_ribbon(
    data = p4_summary,
    aes(x = Time, ymin = LCL, ymax = UCL, group = label),
    fill = "grey", alpha = 0.3
  ) +
  # Predicted median line in blue
  geom_line(
    data = p4_summary,
    aes(x = Time, y = med, group = label),
    color = "blue",
    size = 0.75
  ) +
  # Observed data points
  geom_point(
    data = obs_data,
    aes(x = Time, y = Conc),
    size = 1.5
  ) +
  # Text annotations for missing data
  geom_text(
    data = annot_df,
    aes(x = x, y = y, label = text),
    size = 5,
    color = "grey80"
  ) +
  facet_grid(organs ~ label, scales = "free_y") +
  scale_y_log10(
    limits = c(10^-3, 10^3),
    breaks = trans_breaks("log10", function(x) 10^x, n = 3),
    labels = trans_format("log10", math_format(10^.x))
  ) +
  labs(
    x = "Time (hr)",
    y = "Concentration / Prediction (log scale)"
  ) +
  theme_bw() +
  theme(
    legend.position = "none",
    strip.text = element_text(size = 8)
  ) +
  set_theme  

###################################################################################

# add the title and axis label
title <- ggdraw() +
  draw_label(
    "Rat: Population model fit to calibration data",
    fontface = "plain", #bold
    x = 0,
    size = 15,
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

plot_grid(
  ylab,
  plot_grid(
    title,
    plot_grid(
      plot_grid(p1, p2, nrow = 2, labels = c("I", "II"),
                rel_heights = c(2 / 3, 1 / 3)),
      plot_grid(
        p3, p4, nrow = 2,
        labels = c("III", "IV")
      ),
      nrow = 1, rel_widths = c(0.33, 0.66)
    ),
    xlab, nrow = 3, rel_heights = c(0.05, 1, 0.05)),
  nrow = 1, rel_widths = c(0.02, 1)
)
ggsave(file = "plots/suppl/Supplemental_Figure_S12_population_calibration_Rat.jpg", height = 12, width = 20, dpi = 600)
dev.off()