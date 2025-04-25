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
Mousesim1.1 <- fread("outputs/iHgMouse_3365.out") 
Mousesim2.1 <- fread("outputs/iHgMouse_6734.out") 
Mousesim3.1 <- fread("outputs/iHgMouse_4880.out") 
Mousesim4.1 <- fread("outputs/iHgMouse_5916.out") 

Mouse_x <-mcmc_array(list(Mousesim1.1, Mousesim2.1, Mousesim3.1, Mousesim4.1))

pars_name <- dimnames(Mouse_x)[[3]]
str <- which(pars_name == "M_lnPLC(1)")
end <- which(pars_name == "M_lnKbrnC(1)") 
M_pars <- pars_name[str:end]

str <- which(pars_name == "V_lnPLC(1)")
end <- which(pars_name == "V_lnKbrnC(1)") 
V_pars <- pars_name[str:end]

str <- which(pars_name == "Ve_Aurine(1)")
end <- which(pars_name == "Ve_CBrnU(1)") 
Ve_pars <- pars_name[str:end]

Mouse_mnt <- monitor(Mouse_x[,,c(M_pars, V_pars, Ve_pars)], digit=4)

str <- which(pars_name == "Ve_Aurine(1)")
end <- which(pars_name == "V_lnKbrnC(1)")
parms <- pars_name[str:end]

j <- seq(50001, 100000, 1) 
sum_chains <- length(j)*4 #

n = 500 # 500 virtual study
d <- Mouse_x[j,,] %>% matrix(nrow = sum_chains) %>% as.data.frame() %>% 
  `colnames<-`(dimnames(Mouse_x)[[3]])
i <- which(d[,"LnPosterior"]==max(d[,"LnPosterior"])) %>%
  rep(n) # Pick the sample based on the maximum a posteriori estimation

tmp.Mouse_x <- d[i, parms] 
tmp.Mouse_x %>% write.table(file="iHg_calibration_MAP_Mouse.dat", row.names=T, sep="\t")
vld <- "./mcsim.iHgMouseBW.model.exe MCSim/iHg_calibration_MAP_Mouse.in"
system(vld)

df <- fread("iHg_calibration_MAP_Mouse.out") %>% as.data.frame()
file.remove(c("iHg_calibration_MAP_Mouse.dat", "iHg_calibration_MAP_Mouse.out"))

df

str_1.1 <- which(names(df) == "Aurine_1.1")
end_8.11 <- which(names(df) == "CBrnU_8.11")


#####################################################################
#####################################################################
#####################################################################
# Import Mouse iHg calibration dataset --- same datasets as "iHg_calibration_MAP_Mouse.in"
Time <- c(0.5, 1, 2, 4, 8, 12, 24, 48, 72, 96, 
          
          0.5, 1, 2, 4, 8, 12, 24, 48, 72, 96, 
          
          0.5, 1, 2, 4, 8, 12, 24, 48, 72, 96, 
          
          0.5, 1, 2, 4, 8, 12, 24 , 72 , 240 , 480, 
          
          0.5, 1, 2, 4, 8, 12, 24, 48, 168, 504, 1176, 
          
          0.5, 1, 2, 4, 8, 12, 24, 48, 168, 504, 840, 1176, 
          
          0.5, 1, 2, 4, 8, 12, 24, 48, 168, 504, 840, 1176,  
          
          0.5, 1, 2, 4, 8, 12, 24, 48, 168, 504, 840, 1176, 
          
          0.5, 1, 2, 4, 8, 12, 24, 48, 168, 504, 840, 1176, 
          
          0.5, 1, 2, 4, 8, 12, 24, 48, 168, 504, 840, 1176, 
          
          0.5, 1, 2, 4, 8, 12, 24, 48, 1428, 2868,	4308, 
          
          0.5, 1, 2, 4, 8, 12, 24, 48, 1428, 2868,	4308,
          
          0.5, 1, 2, 4, 8, 12, 24, 48, 1428, 2868,	4308,
          
          0.5, 1, 2, 4, 8, 12, 24, 48, 1428, 2868,	4308, 
          
          0.5, 1, 2, 4, 8, 12, 24, 48, 1428, 2868,	4308,
          
          0.5, 1, 2, 4, 8, 12, 24, 48, 1428, 2868,	4308,
          
          0.5, 1, 2, 4, 8, 12, 24, 48, 1428, 2868,	4308) # Observation time 

Conc <- c(-1, -1, -1, -1, -1, -1, 1.60056272 , 2.220332931 , 2.602229033 , 2.772038576, 
          
          -1,-1,-1,-1,-1,-1,-1,-1, 3.319728,	2.560555, 
          
          -1, -1, -1,-1,-1,-1,-1,-1 ,0.13055565,0.03, 
          
          -1, -1, -1, -1, -1, -1,0.014794515 , 0.005133886 , 0.001263642 , 0.000373918, 
          
          -1, -1, -1, -1, -1, -1, -1, -1, 0.1, 0.1, 0.4, 
          
          -1, -1, -1, -1, -1, -1, -1, -1, 3.9, 3.1, 1, 4, 
          
          -1, -1, -1, -1, -1, -1, -1, -1, 0.2, 0.5, 0.7, 0.6, 
          
          -1, -1, -1, -1, -1, -1, -1, -1, 14.9, 16.9, 22.7, 20.7, 
          
          -1, -1, -1, -1, -1, -1, -1, -1, 0.3, 0.5, 1.1, 1.4, 
          
          -1, -1, -1, -1, -1, -1, -1, -1, 74.9, 68.9, 57.7, 57.7, 
          
          -1, -1, -1, -1, -1, -1, -1, -1, 6.98,	7.02,	6.97,  
          
          -1, -1, -1, -1, -1, -1, -1, -1,1.2,	0.93,	0.8,
          
          -1, -1, -1, -1, -1, -1, -1, -1, 36.43,	26.43,	35.68, 
          
          -1, -1, -1, -1, -1, -1, -1, -1,2.66,	9.71,	2.93,
          
          -1, -1, -1, -1, -1, -1, -1, -1, 111.73,	104.73, 86.88,
          
          -1, -1, -1, -1, -1, -1, -1, -1, 9.88,	9.71,	10.55,
          
          -1, -1, -1, -1, -1, -1, -1, -1, 0.24,	0.35,	0.19) # Observation concentration

formattable(Time , digits = 3, format = "f")
cal_data <- data.frame(Conc, Time)

dose <- c("IV: 400 ug Hg/kg", "Oral: 1000 ug Hg/kg", "oral gavage: 3,000 ug Hg/kg","oral gavage: 15,000 ug Hg/kg", "oral gavage: 75,000 ug Hg/kg", "Oral gavage: 925 ug Hg/kg/d", "Oral gavage: 3695 ug Hg/kg/d", "Oral gavage: 14775 ug/kg/d")

organs <- c(rep("Urine", 10), rep("Kidney", 10),  rep("Brain", 10) , rep("Blood", 10), rep("Blood", 11), rep("Kidney", 12),  rep("Blood", 12), rep("Kidney", 12),  rep("Blood", 12), rep("Kidney", 12),  rep("Kidney", 11),  rep("Liver", 11),   rep("Kidney", 11),  rep("Liver", 11),  rep("Kidney", 11),  rep("Liver", 11), rep("Brain", 11))

Output_Var <- c(rep("Aurine", 10), rep("CKU", 10), rep("CBrnU", 10), rep("CBldU", 10), rep("CBldU", 11),  rep("CKU", 12), rep("CBldU", 12), rep("CKU", 12),  rep("CBldU", 12), rep("CKU", 12), rep("CKU", 11),  rep("CLU", 11), rep("CKU", 11),  rep("CLU", 11), rep("CKU", 11),  rep("CLU", 11),  rep("CBrnU", 11))

study <- c(rep("(Aaseth, 1982))", 30),
           rep("(Nielsen, 1992)", 10), 
           rep("(Dieter, 1983)", 71), 
           rep("(NTP, 1993)",77))

dose <- c(rep("IV: 400 ug Hg/kg", 30),
          rep("Oral: 1000 ug Hg/kg water", 10),
          rep("oral gavage of 3,000 ug Hg/kg", 23),
          rep("oral gavage of 15,000 ug Hg/kg", 24),
          rep("oral gavage of 75,000 ug Hg/kg", 24),
          rep("Oral gavage: 925 ug Hg/kg/d", 22),
          rep("Oral gavage: 3695 ug Hg/kg/d",22),
          rep("Oral gavage: 14775 ug Hg/kg/d",33))

Simulation <- c(
  rep("1", 30),
  rep("2", 10),
  rep("3",23),
  rep("4", 24),
  rep("5",24),
  rep("6", 22),
  rep("7",22),
  rep("8", 33))

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

qt.line <- df[,c(str_1.1:end_8.11)] %>% gather() %>% 
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
  scale_y_log10(
    limits = c(10^-3, 10^2),
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

###################################

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
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x, n = 3),
                labels = trans_format("log10", scales::math_format(10^.x))
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

###################################

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

# Set custom factor levels to control facet order
label_order <- unique(cal_data$label)  

organs_order <- c("Brain","Kidney", "Liver")  

p3_summary$label <- factor(p3_summary$label, levels = label_order)
p3_summary$organs <- factor(p3_summary$organs, levels = organs_order)

obs_data$label <- factor(obs_data$label, levels = label_order)
obs_data$organs <- factor(obs_data$organs, levels = organs_order)


# Plot
p3 <- ggplot() +
  # 95% CI ribbon from simulation
  geom_ribbon(
    data = p3_summary,
    aes(x = Time, ymin = LCL, ymax = UCL, group = label),
    fill = "grey", alpha = 0.3
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
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x, n = 3),
                labels = trans_format("log10", scales::math_format(10^.x))
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

###################################

# f_label4 <- data.frame(tag = c("Oral gavage: 925 ug Hg/kg/d (Ntp, 1993)"), label = c("data not available"), organs = c("Brain"))
# f_label5 <- data.frame(tag = c("Oral gavage: 3,695 ug Hg/kg/d (Ntp, 1993)"), label = c("data not available"), organs = c("Brain"))
# geom_text(aes(x = 2200, y = 0.001, label = label), data = f_label4, size=8, color = "grey", position = position_stack(vjust = .5)) +
# geom_text(aes(x = 2200, y = 0.001, label = label), data = f_label5, size=8, color = "grey", position = position_stack(vjust = .5)) +

    
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
  obs_data <- cal_data %>% filter(Conc > 0, Simulation %in% c(6:8))
  
  # Set custom factor levels to control facet order
  label_order <- unique(cal_data$label)
  organs_order <- c("Brain", "Kidney", "Liver")
  
  p4_summary$label <- factor(p4_summary$label, levels = label_order)
  p4_summary$organs <- factor(p4_summary$organs, levels = organs_order)
  
  obs_data$label <- factor(obs_data$label, levels = label_order)
  obs_data$organs <- factor(obs_data$organs, levels = organs_order)
  
  # Create annotation for missing data in organ = "Brain" and Simulation 6 or 7
  # Identify corresponding labels for Simulations 6 and 7
  labels_sim_6_7 <- cal_data %>%
  filter(Simulation %in% c("6", "7")) %>%
  distinct(label) %>%
  pull(label)
  
  annot_df <- expand.grid(
   label = labels_sim_6_7,
   organs = "Brain"
  ) %>%
   mutate(
     x = 2200,  # position for the text (log scale)
     y = 0.01,
     text = "data not available"
   )

  # Set factor levels to match plot
  annot_df$label <- factor(annot_df$label, levels = label_order)
  annot_df$organs <- factor(annot_df$organs, levels = organs_order)
  
  # Plot
p4 <-   ggplot() +
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
     size = 5, # fontface = "normal",
     color = "grey80"
   ) +
   facet_grid(organs ~ label, scales = "free_y") +
   scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x, n = 5),
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
    
# add the title and axis label
title <- ggdraw() +
  draw_label(
    "Figure S11 (Mouse calibration under the population model)",
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
ggsave(file = "plots/suppl/Supplemental_Figure_S11_population_calibration_Mouse.jpg", height = 12, width = 20, dpi = 600)
dev.off()

