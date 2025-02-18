library(ggpubr)
library(janitor)
library(tidyverse)
library(ggplot2)
library(here)
library(scales)

#################
## Figure 4.2. ##
#################

##Import and clean data
Spike_autoclave_NCTC13441 <- janitor::clean_names(Spiked_experiments_2023)
Spike_autoclave_NCTC13441 <- Spike_autoclave_NCTC13441 %>% mutate(incubation_time = factor(incubation_time, levels = c("4 hrs", "18 hrs")))

## Create plot of stool spiked with ESBL-EC
Spike_autoclave_NCTC13441_plot <- Spike_autoclave_NCTC13441 %>% 
  ggplot(aes(stool, log_cfu, col= stool)) +
  scale_color_manual(values=c('orange','navy'))+
  scale_y_continuous(breaks = seq(5,11,by = 1),
                     limits = c(5,11),
                     labels = math_format()) +
  theme_bw()+
  geom_jitter(shape=16, size=3, position=position_jitter(0.3)) +
  facet_grid (broth ~ incubation_time) + 
  geom_pwc(method = "t_test", label = "p.signif", p.adjust.method = "bonferroni", p.adjust.by = "panel", hide.ns = TRUE, size = 1, label.size = 10, bracket.nudge.y = 0.1) +
  labs(x = NULL, y = "CFU/mL", title = "NCTC13441 spiked into stool") +
  guides(col='none') +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), text = element_text(size=20))

plot(Spike_autoclave_NCTC13441_plot)

ggsave("Spike_autoclave_NCTC13441_plot.png", plot = Spike_autoclave_NCTC13441_plot, device = "png", scale =1, width = 20, height = 20, units = "cm", dpi = 300)

## Import and clean data
Spike_autoclave_CAB17W <- janitor::clean_names(Spiked_experiments_2023)
Spike_autoclave_CAB17W  <- Spike_autoclave_CAB17W %>% mutate(incubation_time = factor(incubation_time, levels = c("4 hrs", "18 hrs")))

## Create plot of autoclaved stool spiked with ESBL-EC
Spike_autoclave_CAB17W_plot <- Spike_autoclave_CAB17W %>% 
  ggplot(aes(stool, log_cfu, col= stool)) +
  scale_color_manual(values=c('orange','navy'))+
  scale_y_continuous(breaks = seq(5,11,by = 1),
                     limits = c(5,11),
                     labels = math_format()) +
  theme_bw()+
  geom_jitter(shape=16, size=3, position=position_jitter(0.3)) +
  geom_pwc(method = "t_test", label = "p.signif", p.adjust.method = "bonferroni", p.adjust.by = "panel", hide.ns = TRUE, size = 1,   label.size = 10, bracket.nudge.y = 0.1) +
  facet_grid (broth ~ incubation_time) + 
  labs(x= NULL, y = "CFU/mL", title = "CAB17W spiked into stool") +
  guides(col='none') +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), text = element_text(size=20))

plot(Spike_autoclave_CAB17W_plot)

ggsave("Spike_autoclave_CAB17W_plot.png", plot = Spike_autoclave_CAB17W_plot, device = "png", scale =1, width = 20, height = 20, units = "cm", dpi = 300)

#################
## Figure 4.3. ##
#################

## Create plot for spiked stool samples incubated for 4 hours
Spike_CAB17W <- janitor::clean_names(Spiked_experiments_2023) 

Spike_CAB17W_plot <- Spike_CAB17W %>% 
  ggplot(aes(isolate, log_cfu, col= broth)) +
  scale_color_manual(values=c('forestgreen','purple'))+
  ylim(6, 9.5) +
  theme_bw()+
  geom_jitter(shape=16, size=2, position=position_jitter(0.1)) +
  facet_grid (antibiotic ~ incubation_time) + 
  labs(y = "log(CFU/mL)") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), text = element_text(size=13))

Spike_CAB17W_plot

ggsave(here("Figure4.2_Spike_CAB17W_plot.png"), plot = Spike_CAB17W_plot, device = "png", scale =1, width = 25, height = 20, units = "cm", dpi = 300)

## Statistical comparison between growth of spiked isolates within two different pre-enrichment with/without antibiotic for spiked stool samples incubated for 4 hours

lm(log_cfu ~ broth + antibiotic, data = Spike_four) -> model_results
summary(model_results_four)

plot(lm(log_cfu ~ broth + antibiotic, data = Spike_four))

coef(model_results_four)
confint(model_results_four)

## Create plot for spiked stool samples incubated for 18 hours
Spike_NCTC13441 <- janitor::clean_names(Spiked_experiments_2023) 

Spike_NCTC13441_plot <- Spike_NCTC13441 %>% 
  ggplot(aes(isolate, log_cfu, col= broth)) +
  scale_color_manual(values=c('forestgreen','purple'))+
  ylim(6, 9.5) +
  theme_bw()+
  geom_jitter(shape=16, size=2, position=position_jitter(0.1)) +
  facet_grid (antibiotic ~ incubation_time) + 
  labs(y = "log(CFU/mL)") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), text = element_text(size=13))

Spike_NCTC13441_plot

ggsave(here("Figure4.2_Spike_NCTC13441_plot.png"), plot = Spike_NCTC13441_plot, device = "png", scale =1, width = 25, height = 20, units = "cm", dpi = 300)

##Statistical comparison between growth of spiked isolates within two different pre-enrichment with/without antibiotic for spiked stool samples incubated for 18 hours

spike_stats <- aov(log_cfu ~ broth + antibiotic, Spike_eighteen)

summary(spike_stats)

plot(lm(log_cfu ~ broth + antibiotic, data = Spike_eighteen))

Post_hoc <- TukeyHSD(spike_stats)
Post_hoc

coef(model_results)
confint(model_results)

#################
## Figure 4.4. ##
#################

##Import and clean data
LOD <- janitor::clean_names(Spiked_experiments_2023)

## Create plot for lowest detectable limit
LOD_plot <- LOD %>% ggplot(aes(x= input_log_cfu, y= log_cfu, col = method)) +
  geom_point(cex = 3) +
  geom_smooth(method = "lm", formula = y ~ x, se=FALSE, fullrange = TRUE) +
  scale_color_manual(values= c("Direct plating" = "darkgray", "Pre-enrichment" = 'forestgreen')) +
  stat_regline_equation(cex= 8, show.legend = FALSE) +
  scale_x_continuous(breaks = seq(-11,7,by = 1),
                     limits = c(-11,7),
                     labels = math_format()) +
  scale_y_continuous(breaks = seq(-11,6,by = 1),
                     limits = c(-11,6),
                     labels = math_format()) +
  geom_hline(yintercept=0, color = "gray") +
  guides(col = guide_legend("Method")) +
  labs(x = "Input (CFU/mL)", y = "Colony count (CFU/mL)") +
  theme_bw() +
  theme(text = element_text(size=20))

plot(LOD_plot)

ggsave("LOD_plot.tiff", plot = LOD_plot, device = "tiff", scale =1, width = 25, height = 20, units = "cm", dpi = 300)
