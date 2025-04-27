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
Humansim1.1 <- fread("outputs/iHgHuman_3365.out") 
Humansim2.1 <- fread("outputs/iHgHuman_6734.out") 
Humansim3.1 <- fread("outputs/iHgHuman_4880.out") 
Humansim4.1 <- fread("outputs/iHgHuman_5916.out") 

Human_x <-mcmc_array(list(Humansim1.1, Humansim2.1, Humansim3.1, Humansim4.1))

pars_name <- dimnames(Human_x)[[3]]
str <- which(pars_name == "M_lnPLC(1)")
end <- which(pars_name == "M_lnKbrnC(1)") 
M_pars <- pars_name[str:end]

str <- which(pars_name == "V_lnPLC(1)")
end <- which(pars_name == "V_lnKbrnC(1)") 
V_pars <- pars_name[str:end]

str <- which(pars_name == "Ve_Aurine(1)")
end <- which(pars_name == "Ve_ABld(1)") 
Ve_pars <- pars_name[str:end]

Human_mnt <- monitor(Human_x[,,c(M_pars, V_pars, Ve_pars)], digit=4)

str <- which(pars_name == "Ve_Aurine(1)")
end <- which(pars_name == "V_lnKbrnC(1)")
parms <- pars_name[str:end]

j <- seq(150001, 300000, 1) 
sum_chains <- length(j)*4 #

n = 500 # 500 virtual study
d <- Human_x[j,,] %>% matrix(nrow = sum_chains) %>% as.data.frame() %>% 
  `colnames<-`(dimnames(Human_x)[[3]])
i <- which(d[,"LnPosterior"]==max(d[,"LnPosterior"])) %>%
  rep(n) # Pick the sample based on the maximum a posteriori estimation

tmp.Human_x <- d[i, parms] 
tmp.Human_x %>% write.table(file="iHg_calibration_MAP_Human.dat", row.names=T, sep="\t")
vld <- "./mcsim.iHgHumanBW.model.exe MCSim/iHg_calibration_MAP_Human.in"
system(vld)

df <- fread("iHg_calibration_MAP_Human.out") %>% as.data.frame()
file.remove(c("iHg_calibration_MAP_Human.dat", "iHg_calibration_MAP_Human.out"))

df

str_1.1 <- which(names(df) == "Afeces_1.1")
end_2.23 <- which(names(df) == "Aurine_2.23")

#####################################################################
#####################################################################
#####################################################################
# Import Human iHg calibration dataset --- same datasets as "iHg_calibration_MAP_Human.in"
Time <- c( 24 , 72 , 96 , 120 , 144 , 168 , 216 , 240 , 264 , 288 , 336 , 384 , 408 , 456 , 480 , 504 , 528 , 600 , 624 , 696 , 768 , 840 , 888 , 912 , 1032 , 1080 , 1272 , 1296 , 1368 , 1392 , 1488 , 1512 , 1608 , 1632 , 1680 , 1704, 
          
           24 , 48 , 72 , 120 , 144 , 168 , 192 , 216 , 240 , 288 , 312 , 336 , 360 , 384 , 432 , 504 , 552 , 600 , 648 , 696 , 720 , 744 , 768 , 840 , 912 , 960 , 1008 , 1032 , 1056 , 1080 , 1128 , 1152 , 1200 , 1224 , 1272 , 1296 , 1344 , 1392 , 1488 , 1536 , 1560 , 1656 , 1680 , 1704, 
          
           12 , 24 , 48 , 72 , 96 , 120 , 144 , 168 , 192 , 240 , 264 , 336 , 600 , 792 , 888 , 1032 , 1152 , 1320 , 1680, 

           24 , 48 , 72 , 96 , 120 , 144 , 192 , 216 , 240 , 288 , 336 , 384 , 408 , 432 , 456 , 528 , 672 , 720 , 1224 , 1248 , 1272,
           
           24 , 48 , 72 , 96 , 120 , 144 , 192 , 216 , 240 , 288 , 336 , 384 , 408 , 432 , 456 , 528 , 672 , 720 , 1224 , 1248 , 1272 , 1800 , 1824) # Observation time 

Conc <-  c(0.0226975 , 0.1563625 , 0.20055 , 0.235025 , 0.24927 , 0.267085 , 0.2848825 , 0.2991275 , 0.3101525 , 0.32053 , 0.332115 , 0.34923 , 0.358925 , 0.3651725 , 0.3696525 , 0.3706325 , 0.3952725 , 0.3813075 , 0.386505 , 0.4273325 , 0.4228875 , 0.4131225 , 0.42378 , 0.434595 , 0.4587275 , 0.5030025 , 0.512785 , 0.4764375 , 0.498925 , 0.4865875 , 0.491785 , 0.4970525 , 0.5020925 , 0.544845 , 0.5245975 , 0.53382, 
          
          0.0035525 , 0.0123725 , 0.02919 , 0.0380975 , 0.0454125 , 0.0574175 , 0.05908 , 0.06489 , 0.06615 , 0.077805 , 0.100975 , 0.1026375 , 0.10675 , 0.1381975 , 0.143605 , 0.14931 , 0.144025 , 0.14945 , 0.2083375 , 0.21938 , 0.2232475 , 0.2345875 , 0.2416575 , 0.2433025 , 0.222985 , 0.23198 , 0.2997925 , 0.306495 , 0.306005 , 0.31682 , 0.3254125 , 0.3291925 , 0.33838 , 0.3493525 , 0.351435 , 0.35567 , 0.3709125 , 0.3670275 , 0.382165 , 0.3973025 , 0.4137 , 0.3668875 , 0.42251, 0.4321625, 
          
          0.188265 , 0.0868 , 0.0531825 , 0.03969 , 0.026355 , 0.0270375 , 0.0196 , 0.0168 , 0.02051 , 0.0112175 , 0.0131425 , 0.0079275 , 0.0074375 , 0.0065625 , 0.0037625 , 0.0069475 , 0.0054775 , 0.00301 , 0.003045, 
          
          1.474756451 , 3.463299951 , 4.714515346 , 5.261447416 , 5.508576199 , 5.654103023 , 5.762256363 , 5.782233664 , 5.799170469 , 5.814025431 , 5.823694268 , 5.830198949 , 5.835715647 , 5.843403469 , 5.84711481 , 5.851088374 , 5.852221182 , 5.853751437 , 5.854550353 , 5.855493808 , 5.856571646,
          
          0.001803294 , 0.004636915 , 0.006745647 , 0.009028727 , 0.011016002 , 0.012919823 , 0.016757153 , 0.019013174 , 0.020972077 , 0.023108392 , 0.028847792 , 0.031639119 , 0.03636045 , 0.04322589 , 0.046753491 , 0.05153253 , 0.057786306 , 0.065201937 , 0.066729829 , 0.06898874 , 0.070779686 , 0.075422294 , 0.078354255) # Observation concentration

formattable(Time , digits = 3, format = "f")
cal_data <- data.frame(Conc, Time)

# study <- c("Hall 1995", "Rahola, 1973")

dose <- c("IV: 0.025 ug Hg/kg/d", "Oral: 0.09375 ug Hg/kg")

organs <- c(rep("Feces", 36),  rep("Urine", 44), rep("Blood", 19),  rep("Feces", 21),  rep("Urine", 23))

Output_Var <- c(rep("Afeces", 36), rep("Aurine", 44), rep("ABld", 19), rep("Afeces", 21), rep("Aurine", 23))

study <- c(rep("(Hall 1995)", 99),
           rep("(Rahola, 1973)",44))

dose <- c(rep("single IV dose of 0.025 ug Hg/kg/d", 99),
          rep("single oral of 0.09375 ug Hg/kg: 3695 ug Hg/kg/d",44))

Simulation <- c(
  rep("1", 99),
  rep("2",44))

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

qt.line <- df[,c(str_1.1:end_2.23)] %>% gather() %>% 
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
  scale_y_log10(lim = c(10^-4, 10^1),
                breaks = trans_breaks("log10", function(x) 10^x, n = 3),
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
  scale_y_log10(lim = c(10^-3, 10^1.5),
                breaks = trans_breaks("log10", function(x) 10^x, n = 4),
                labels = trans_format("log10", scales::math_format(10^.x)))  +
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
    "Figure S13 (Human calibration under the population model)",
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
    "Hg in blood, or excretion in urine or feces (ug)",
    fontface = "bold", size = 14, vjust = 0, angle = 90
  ) + theme(
    plot.margin = margin(0, 0, 0, 1)
  )

# plot

plot_grid(
  ylab,
  plot_grid(
    title,
    plot_grid(p1, p2, nrow = 1, labels = c("I", "II"), rel_widths = c(0.5, 0.5)),
    xlab, nrow = 3, rel_heights = c(0.05, 1, 0.05)
  ),
  nrow = 1, rel_widths = c(0.02, 1)
)
ggsave(file = "plots/suppl/Supplemental_Figure_S13_population_calibration_Human.jpg", height = 12, width = 20, dpi = 600)
dev.off()

