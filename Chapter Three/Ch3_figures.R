library(ggpubr)
library(janitor)
library(ggrepel)
library(tidyverse)
library(dplyr)
library(ggplot2)
library(growthcurver)
library(purrr)
library(here)

#################
## Figure 3.2. ##
#################

##Import and clean data
Growthcurve <- read.csv(here("Fig3.2_data.csv"))
Growthcurve <- janitor::clean_names(Broth_compare)

## Plot simple scatter plot for colony counts over time
Growthcurve_plot <- Growthcurve %>% 
  ggplot(aes(x= time_hours, y= log_cfu, col= broth)) +
  scale_color_manual(values= c("DM" = "#89d2d4", "BHI" = "#f7867e", "BPW" = "#b7cd88", "TS"= "#d677cc")) +
  geom_point(size=2, alpha= 0.5) +
  labs(x = "Incubation Time (Hours)", y = "cfu/mL")

## Create subset data plots and predicted growth curves by fitting population logistic curves to each pre-enrichment broth

bpw.predicted <- data.frame(
  time_hours =
    subset(Broth_compare, broth == "BPW")$time_hours,
  pred.bpw = predict(bpw_fit$model))

dm.predicted <- data.frame(
  time_hours =
    subset(Broth_compare, broth == "DM")$time_hours,
  pred.dm = predict(dm_fit$model))

ts.predicted <- data.frame(
  time_hours =
    subset(Broth_compare, broth == "TS")$time_hours,
  pred.ts = predict(ts_fit$model))

bhi.predicted <- data.frame(
  time_hours =
    subset(Broth_compare, broth == "BHI")$time_hours,
  pred.bhi = predict(bhi_fit$model))

## Fit predicted population logistic growth curves to simple scatter plot
Growthcurve_plot <- 
  Growthcurve_plot + 
  geom_line(data=bpw.predicted, aes(time_hours, pred.bpw), color="red")+
  geom_line(data=dm.predicted, aes(time_hours, pred.dm), color="forestgreen")+
  geom_line(data=ts.predicted, aes(time_hours, pred.ts), color="magenta")+
  geom_line(data=bhi.predicted, aes(time_hours, pred.bhi), color="steelblue3")

ggsave(here("Figure3.2_Growthcurve.png"), plot = Growthcurve_plot, device = "png", scale =1, width = 25, height = 20, units = "cm", dpi = 300)

#################
## Figure 3.3. ##
#################

## Import and clean data
ODvLOG <- read.csv(here("Fig3.3_data.csv"))
ODvLOG <- janitor::clean_names(All_broths) 

## Linear regression
ODvLOG_plot <- ODvLOG %>% ggplot(aes(od600, log_cfu_m_l, col = broth)) +
  geom_point() +
  guides(col = guide_legend("Broth")) +
  scale_color_manual(values= c("DM" = "#89d2d4", "BHI" = "#f7867e", "BPW" = "#b7cd88", "TS"= "#d677cc")) +
  geom_smooth(method = "lm", formula = y ~ x + 0 , se= FALSE)+
  stat_regline_equation(cex= 8, show.legend = FALSE) +
  labs(x = "OD600", y = "log(CFU/mL)") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme_bw() +
  theme(text = element_text(size=20))

All_broth_plot

ggsave(here("Figure3.3_ODvLOG.png"), plot = All_broth_plot, device = "png", scale =1, width = 25, height = 20, units = "cm", dpi = 300)

#################
## Figure 3.4. ##
#################

## Import and clean data
Broth_compare <- read.csv(here("Fig3.4_data.csv"))
Broth_compare <- janitor::clean_names(Broth_compare) 

## Plot ESBL-EC recovery on agar after pre-enrichment
Broth_compare_plot <- Broth_compare %>% ggplot(aes(x= broth, y= log_cfu, col = broth)) + 
  scale_color_manual(values= c("DM" = "#89d2d4", "BHI" = "#f7867e", "BPW" = "#b7cd88", "TS"= "#d677cc"))+
  geom_jitter()+
  facet_wrap(~ time, ncol = 2) +
  geom_boxplot(alpha = 0.5, lwd = 0.2, notch = FALSE)+
  labs(x= "Broth", y = "logCFU/mL") + 
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), text = element_text(size=13))

Broth_compare_plot

## Save plot
ggsave(here("Figure3.4_Broth_compare.pdf"), plot = Broth_compare_plot, device = "pdf", scale =1, width = 20, height = 25, units = "cm", dpi = 300)

#################
## Figure 3.5. ##
#################

## Import and clean data
Qubit_yield <- read.csv(here("Fig3.5_data.csv"))
Qubit_yield <- janitor::clean_names(R_analysis)
Qubit_yield <- Qubit_yield %>% mutate(extraction_kit = factor(extraction_kit, levels = c("Boiling", "NEB Monarch", "Lucigen MasterPure", "Promega Wizard", "Qiagen Dneasy", "Zymo Miniprep")))

## Plot DNA yield
Qubit_yield_plot <- Qubit_yield %>% ggplot(aes(extraction_kit, dna_yield, color= extraction_kit)) + scale_color_manual(values= c("Boiling" = "#5A5A5A", "Qiagen Dneasy" = "#c4567c", "Promega Wizard" = "#6ccff6", "Lucigen MasterPure" = "#b5838d", "Zymo Miniprep" = "#B5F58C", "NEB Monarch"= "#6bb4d1")) +
  geom_boxplot(alpha = 0.5, notch = FALSE, outlier.shape = NA) + 
  geom_jitter(shape=16, position=position_jitter(0.2)) +
  labs(x = "Extraction Kit", y = "DNA yield (ng/uL)") + 
  guides(color = 'none') +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), text = element_text(size=20))

Qubit_yield_plot

## Save plot
ggsave("Figure3.5_Qubit_yield_plot.tiff", plot = Qubit_yield_plot, device = "tiff", scale =1, width = 25, height = 20, units = "cm", dpi = 300)

#################
## Figure 3.6. ##
#################

## Import and clean data
Nanodrop.260.230 <- read.csv(here("Fig3.6A_data.csv"))
Nanodrop.260.230 <- janitor::clean_names(R_analysis)
Nanodrop.260.230 <- Nanodrop.260.230 %>% mutate(extraction_kit = factor(extraction_kit, levels = c("Boiling", "Lucigen MasterPure", "NEB Monarch", "Promega Wizard", "Qiagen Dneasy", "Zymo Miniprep")))

## Plot DNA quality (260/230)
Nanodrop_260_230_plot <- Nanodrop.260.230 %>% ggplot(aes(extraction_kit, purity, color = extraction_kit)) + 
  scale_color_manual(values= c("Boiling" = "#5A5A5A", "Qiagen Dneasy" = "#c4567c", "Promega Wizard" = "#6ccff6", "Lucigen MasterPure" = "#b5838d", "Zymo Miniprep" = "#B5F58C", "NEB Monarch"= "#6bb4d1")) +
  geom_rect(aes(xmin=-Inf,xmax=Inf,ymin=2.0,ymax=3.0), color= NA, fill="grey", alpha = 0.3) +
  geom_boxplot(alpha = 0.5, notch = FALSE, outlier.shape = NA) +
  geom_jitter(shape=16, position=position_jitter(0.2)) +
  labs(x= "Extraction Kit", y = "A260/A230")+ 
  guides(color = 'none') +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), text = element_text(size=15))

## Import and clean data
Nanodrop.260.280 <- read.csv(here("Fig3.6B_data.csv"))
Nanodrop.260.280 <- janitor::clean_names(R_analysis)
Nanodrop.260.280 <- Nanodrop.260.280 %>% mutate(extraction_kit = factor(extraction_kit, levels = c("Boiling", "Lucigen MasterPure", "NEB Monarch", "Promega Wizard", "Qiagen Dneasy", "Zymo Miniprep")))

## Plot DNA quality (260/280)
Nanodrop_260_280_plot <- Nanodrop.260.280 %>% ggplot(aes(extraction_kit, contam, colour = extraction_kit)) + 
  scale_color_manual(values= c("Boiling" = "#5A5A5A", "Qiagen Dneasy" = "#c4567c", "Promega Wizard" = "#6ccff6", "Lucigen MasterPure" = "#b5838d", "Zymo Miniprep" = "#B5F58C", "NEB Monarch"= "#6bb4d1")) +
  geom_rect(aes(xmin=-Inf,xmax=Inf,ymin=1.8,ymax=2.1),color=NA, fill="grey", alpha = 0.3) + 
  geom_boxplot(alpha = 0.5, notch = FALSE, outlier.shape = NA)+ 
  geom_jitter(shape=16, position=position_jitter(0.2))+
  labs(x= "Extraction Kit", y = "A260/A280")+ 
  guides(color = 'none') +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), text = element_text(size=15))

## Combine plots
DNA_quality_plots <- ggarrange(Nanodrop_260_280_plot, Nanodrop_260_230_plot, labels = c("A", "B"), ncol = 1, nrow = 2)

DNA_quality_plots

## Save plot
ggsave("Figure3.6_DNA_quality_plots.png", plot = DNA_quality_plots, device = "png", scale =0.5, width = 30, height = 50, units = "cm", dpi = 300)

#################
## Figure 3.7. ##
#################

## Save colours for plots
colours <- c("Qiagen Dneasy" = "#c4567c", 
             "Lucigen MasterPure" = "#b5838d", 
             "Promega Wizard" = "#6ccff6", 
             "NEB Monarch"= "navyblue", 
             "Zymo Miniprep" = "#B5F58C",  
             "Boiling" = "#5A5A5A")

## Set radar chart format
radarchart(Kit_comparison_scoring_matrix, cglty = 1, cglcol = "gray", 
           pcol = colours,plwd = 1, plty = 1, cglwd = 1, vlcex = 0.8) 

Extraction_score_plot <- plot_ly(type = 'scatterpolar', fill='toself')

Extraction_score_plot <- 
  Extraction_score_plot %>% layout(polar = list(radialaxis = list(visible = T, range = c(0,10), 
                                                                  tickfont = list(size = 20)), 
                                                                  angularaxis = list(tickfont = list(size = 20))))

## Plot radar charts for extraction kit scores
Dneasy_score_plot <- 
  Extraction_score_plot %>%
  add_trace(r = c(1, 7, 1, 5, 7),
            theta = c('Cost','Time','Yield', 'Quality', 'Difficulty'),name = 'Qiagen Dneasy', color = I("green"))

Dneasy_score_plot

MasterPure_score_plot <- 
  Extraction_score_plot %>%
  add_trace(r = c(7, 6, 10, 7, 5),
            theta = c('Cost','Time','Yield', 'Quality', 'Difficulty'),name = 'Lucigen MasterPure', color = I("orange"))

MasterPure_score_plot

Promega_score_plot <- 
  Extraction_score_plot %>%
  add_trace(r = c(1, 6, 1, 7, 6),
            theta = c('Cost','Time','Yield', 'Quality', 'Difficulty'),name = 'Promega Wizard', color = I("#6ccff6"))

Promega_score_plot

Monarch_score_plot <- 
  Extraction_score_plot %>%
  add_trace(r = c(5, 9, 1, 6, 9),
            theta = c('Cost','Time','Yield', 'Quality', 'Difficulty'),name ='NEB Monarch', color = I("navyblue"))

Monarch_score_plot

Miniprep_score_plot <- 
  Extraction_score_plot %>%
  add_trace(r = c(2, 7, 1, 7, 7),
            theta = c('Cost','Time','Yield', 'Quality', 'Difficulty'),name = 'Zymo Miniprep', color = I("yellow"))

Miniprep_score_plot

Boiling_score_plot <- 
  Extraction_score_plot %>%
  add_trace(r = c(10, 10, 10, 1, 10),
            theta = c('Cost','Time','Yield', 'Quality', 'Difficulty'),name = 'Boiling', color = I("darkgrey"))

Boiling_score_plot

## Save plot
ggsave("Figure3.7_Dneasy_score_plot.png", plot = Dneasy_score_plot, device = "png", scale =0.5, width = 30, height = 50, units = "cm", dpi = 300)
ggsave("Figure3.7_MasterPure_score_plot.png", plot = MasterPure_score_plot, device = "png", scale =0.5, width = 30, height = 50, units = "cm", dpi = 300)
ggsave("Figure3.7_Promega_score_plot.png", plot = Promega_score_plot, device = "png", scale =0.5, width = 30, height = 50, units = "cm", dpi = 300)
ggsave("Figure3.7_Monarch_score_plot.png", plot = Monarch_score_plot, device = "png", scale =0.5, width = 30, height = 50, units = "cm", dpi = 300)
ggsave("Figure3.7_Miniprep_score_plot.png", plot = Miniprep_score_plot, device = "png", scale =0.5, width = 30, height = 50, units = "cm", dpi = 300)
ggsave("Figure3.7_Boiling_score_plot.png", plot = Boiling_score_plot, device = "png", scale =0.5, width = 30, height = 50, units = "cm", dpi = 300)