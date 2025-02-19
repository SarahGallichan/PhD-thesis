library(ggpubr)
library(janitor)
library(tidyverse)
library(ggplot2)
library(here)

#################
## Figure 5.2. ##
#################

## Import and clean data
Stool1_species <- read.csv(here("Fig5.2_data.csv"))
Stool1_species <- janitor::clean_names(Stool1_species) 
Stool1_species <- Stool1_species %>% mutate(species = factor(species, levels = c("Other", "Escherichia coli")))


Stool1_species_plot <- Stool1_species %>% ggplot(aes(x = condition, y = relative_abundance, fill = species, label = reads)) +
  geom_col() +
  geom_text(size = 4, position = position_stack(vjust = 0.7)) +
  guides(fill=guide_legend(title="Species")) +
  scale_fill_manual(values= c("Other" = "#8d96a3", "Escherichia coli" = "#d1495b")) +
  labs(title = "Stool 1", x = "Protocol", y = "Relative Abundance (%)") +
  scale_x_discrete(limits = c('Untreated', '4hrs BPW', '4hrs TS', '18hrs BPW', '18hrs TS' )) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), text = element_text(size=20))

plot(Stool1_species_plot)
ggsave("Stool1_species_plot.png", plot = Stool1_species_plot, device = "png", scale =1, width = 25, height = 20, units = "cm", dpi = 300)

## Import and clean data
Stool2 <- janitor::clean_names(Braken_species) 
Stool2 <- Stool2 %>% mutate(species = factor(species, levels = c("Other", "Raoultella planticola", "Klebsiella pneumoniae", "Escherichia coli")))

Stool_plot2 <- Stool2 %>% ggplot(aes(x = condition, y = relative_abundance, fill = species, label = reads)) +
  geom_col() +
  geom_text(size = 4, position = position_stack(vjust = 0.7)) +
  guides(fill=guide_legend(title="Species")) +
  scale_fill_manual(values= c("Other" = "#8d96a3", "Raoultella planticola" = "#66a182", "Klebsiella pneumoniae" = "#edae49", "Escherichia coli" = "#d1495b")) +
  labs(title = 'Stool 2', x = "Protocol", y = "Relative Abundance (%)") +
  scale_x_discrete(limits = c('Untreated', '4hrs BPW', '4hrs TS', '18hrs BPW', '18hrs TS' )) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), text = element_text(size=20))

plot(Stool_plot2)
ggsave("Stool_plot2.png", plot = Stool_plot2, device = "png", scale =1, width = 25, height = 20, units = "cm", dpi = 300)

## Import and clean data
Stool3 <- janitor::clean_names(Braken_species) 
Stool3 <- Stool3 %>% mutate(species = factor(species, levels = c("Other", "Escherichia coli")))

Stool_plot3 <- Stool3 %>% ggplot(aes(x = condition, y = relative_abundance, fill = species, label = reads)) +
  geom_col() +
  geom_text(size = 4, position = position_stack(vjust = 0.7)) +
  guides(fill=guide_legend(title="Species")) +
  scale_fill_manual(values= c("Other" = "#8d96a3", "Escherichia coli" = "#d1495b")) +
  labs(title = "Stool 3", x = "Protocol", y = "Relative Abundance (%)") +
  scale_x_discrete(limits = c('Untreated', '4hrs BPW', '4hrs TS', '18hrs BPW', '18hrs TS' )) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), text = element_text(size=20))

plot(Stool_plot3)
ggsave("Stool_plot3.png", plot = Stool_plot3, device = "png", scale =1, width = 25, height = 20, units = "cm", dpi = 300)

#################
## Figure 5.3. ##
#################

Stool1_6mil <- read.csv(here("Fig5.3_data.csv"))
Stool1_6mil <- janitor::clean_names(Braken_species) 
Stool1_6mil <- Stool1_6mil %>% mutate(species = factor(species, levels = c("Other", "Escherichia coli")))

Stool_plot1_6mil <- Stool1_6mil %>% ggplot(aes(x = condition, y = relative_abundance, fill = species, label = reads)) +
  geom_col() +
  geom_text(size = 4, position = position_stack(vjust = 0.7)) +
  guides(fill=guide_legend(title="Species")) +
  scale_fill_manual(values= c("Other" = "#8d96a3", "Escherichia coli" = "#d1495b")) +
  labs(title = "Stool 1", x = "Protocol", y = "Relative Abundance (%)") +
  scale_x_discrete(limits = c('Untreated', '4hrs BPW', '18hrs BPW' )) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), text = element_text(size=20))

plot(Stool_plot1_6mil)
ggsave("Stool_plot1_6mil.png", plot = Stool_plot1_6mil, device = "png", scale =1, width = 25, height = 20, units = "cm", dpi = 300)

#################
## Figure 5.4. ##
#################

##Import and clean data
Stool1_strains <- read.csv(here("Fig5.4_Stool1_data.csv"))
Stool1_strains <- janitor::clean_names(Strain_results) 
Stool1_strains <- Stool1_strains %>% mutate(duration = factor(duration, levels = c("0", "4", "18")))

Stool1_strains_plot <- Stool1_strains %>% ggplot(aes(x = duration, y = relative_abundance, fill = e_coli_strain)) +
  geom_col() +
  scale_fill_manual(values= c("Non E. coli" = "#d3d3d3", "ST399" = "#cae7b9", "ST216" = "#1c5253", "ST10" = "#6ccff6", "ST131" = "#a23e48", "ST543" = "#dd7230", "ST906" = "#f4c95d", "ST28"= "#7871aa")) +
  facet_wrap('broth', ncol = 3) +
  guides(fill = guide_legend(title="E.coli Sequence type")) +
  labs(x = "Pre-enrichment duration (hrs)", y = "Relative Abundance (%)") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), text = element_text(size=20))

plot(Stool1_strains_plot)

##Import and clean data
Stool3_species <- read.csv(here("Fig5.4_Stool3_data.csv"))
Stool3_strains <- janitor::clean_names(Strain_results) 
Stool3_strains <- Stool3_strains %>% mutate(duration = factor(duration, levels = c("0", "4", "18")))

Stool3_strains_plot <- Stool3_strains %>% ggplot(aes(x = duration, y = relative_abundance, fill = e_coli_strain)) +
  geom_col() +
  scale_fill_manual(values= c("Non E. coli" = "#d3d3d3", "ST472" = "#061148", "ST635" = "#689695", "ST162" = "#b5838d", "ST491" = "#a74a43", "ST1193" = "#8d94c5", "ST69" = "#f4c95d", "ST28"= "#6bb4d1")) +
  facet_wrap('broth', ncol = 2) +
  guides(fill = guide_legend(title="E.coli Sequence type")) +
  labs(x = "Pre-enrichment duration (hrs)", y = "Relative Abundance (%)") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), text = element_text(size=20))

plot(Stool3_strains_plot)

#################
## Figure 5.5. ##
#################

##Import and clean data
Stool1_CoveragevCondition <- read.csv(here("Fig5.5_Stool1_data.csv"))
Stool1_CoveragevCondition <- janitor::clean_names(Strain_results)
Stool1_CoveragevCondition <- Stool1_CoveragevCondition %>% mutate(incubation_time = factor(incubation_time, levels = c("4 hours", "18 hours")))

Stool1_Coverage <- 
  Stool1_CoveragevCondition %>% ggplot(aes(reference_strain, coverage, col = incubation_time)) +
  facet_wrap(~ pre_enrichment, ncol = 2) +
  geom_point(aes(y=coverage),size = 3) +
  guides(col = guide_legend("Incubation time")) +
  scale_color_manual(values= c("4 hours" = "#8d94c5", "18 hours" = "#f4c95d")) +
  labs(x = "Reference Strain", y = "Coverage") +
  ylim(0.0, 2.0) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 70, vjust = 1, hjust=1), text = element_text(size=15))

Stool1_Coverage

Stool1_Callable <- 
  Stool1_CoveragevCondition %>% ggplot(aes(reference_strain, callable, fill = incubation_time)) +
  facet_wrap(~ pre_enrichment, ncol = 2) +
  geom_col(aes(fill = incubation_time), alpha = 0.5) +
  guides(fill = guide_legend("Incubation time")) +
  scale_fill_manual(values= c("4 hours" = "#8d94c5", "18 hours" = "#f4c95d")) +
  ylim(0, 50) +
  labs(x = "Reference Strain", y = "Callable (%)", Title = "Stool 1") +
  theme_bw() +
  theme(axis.title.x = element_blank(), axis.text.x = element_blank(), text = element_text(size=15))

Stool1_Callable

Stool1_Call_v_cov <- ggarrange(Stool1_Callable, Stool1_Coverage, ncol = 1, nrow = 2, common.legend = TRUE, legend = "right")

## Import and clean data
Stool3_CoveragevCondition <- read.csv(here("Fig5.5_Stool3_data.csv"))
Stool3_CoveragevCondition <- janitor::clean_names(Strain_results)
Stool3_CoveragevCondition <- Stool3_CoveragevCondition %>% mutate(incubation_time = factor(incubation_time, levels = c("4 hours", "18 hours")))

Stool3_Coverage <- 
  Stool3_CoveragevCondition %>% ggplot(aes(reference_strain, coverage, col = incubation_time)) +
  facet_wrap(~ pre_enrichment, ncol = 2) +
  geom_point(aes(y=coverage),size = 3) +
  guides(col = guide_legend("Incubation time")) +
  scale_color_manual(values= c("4 hours" = "#8d94c5", "18 hours" = "#f4c95d")) +
  labs(x = "Reference Strain", y = "Coverage", Title = "Stool 3") +
  ylim(0.0, 2.0) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 70, vjust = 1, hjust=1), text = element_text(size=15))

plot(Stool3_Coverage)

Stool3_Callable <- 
  Stool3_CoveragevCondition %>% ggplot(aes(reference_strain, callable, fill = incubation_time)) +
  facet_wrap(~ pre_enrichment, ncol = 2) +
  geom_col(aes(fill = incubation_time), alpha = 0.5) +
  guides(fill = guide_legend("Incubation time")) +
  scale_fill_manual(values= c("4 hours" = "#8d94c5", "18 hours" = "#f4c95d")) +
  labs(x = "Reference Strain", y = "Callable (%)") +
  ylim(0, 50) +
  theme_bw() +
  theme(axis.title.x = element_blank(), axis.text.x = element_blank(), text = element_text(size=15))

plot(Stool3_Callable)

Stool3_Call_v_cov <- ggarrange(Stool3_Callable, Stool3_Coverage, ncol = 1, nrow = 2, common.legend = TRUE, legend = "right")

#################
## Figure 5.6. ##
#################

## Import and clean data
Reference_compare <- read.csv(here("Fig5.6_Refs_data.csv"))
ref_compare_matrix <- Reference_compare
as.data.frame(ref_compare_matrix) -> ref_compare_matrix
rownames(ref_compare_matrix) <- ref_compare_matrix[,1]

Phylogroup <- read.csv(here("Fig5.6_Phylogroup_data.csv"))
phylogroup <- Phylogroup
as.data.frame(phylogroup) -> phylogroup
rownames(phylogroup) <- phylogroup$Reference

ref_compare_matrix <- ref_compare_matrix[,-1]

my_colours <- phylogroup %>% list(Phylogroup = c(A = "#216461", B1="#d1e2a9", B2="#66883f", D = "#8dba24"))

ref_heatmap <- pheatmap(ref_compare_matrix,  annotation_row = phylogroup[1], annotation_colors = my_colours, clustering_method = 'ward.D', angle_col = "90", treeheight_col = 0, treeheight_row = 0, fontsize = 15)

#################
## Figure 5.7. ##
#################

Phylogroup <- read.csv(here("Fig5.7_data.csv"))
Stool2Klebs_strains <- janitor::clean_names(Strain_results) 

Stool2Klebs_strains <- Stool2Klebs_strains %>% 
         mutate(duration = factor(duration, levels = c("0", "4", "18")), 
         sequence_type_strain = factor(sequence_type_strain, 
         levels = c("Other species",  "Klebsiella pneumonia ST1117", "Klebsiella pneumonia ST35", 
                    "Klebsiella pneumonia ST37", "Klebsiella pneumonia ST3750", "Klebsiella pneumonia ST268", 
                    "Raoultella planticola S25")))

Stool2Klebs_strains_plot <- Stool2Klebs_strains %>% ggplot(aes(x = duration, y = relative_abundance, fill = sequence_type_strain)) +
  geom_col() +
  scale_fill_manual(values= c("Other species" = "#d3d3d3", "Klebsiella pneumonia ST1117" = "#F9F871", "Klebsiella pneumonia ST35" = "#D65DB1", 
                              "Klebsiella pneumonia ST37" = "#B5F58C", "Klebsiella pneumonia ST3750" = "#8d94c5", "Klebsiella pneumonia ST268" = "#6bb4d1", 
                              "Raoultella planticola S25" = "#FF9671")) +
  facet_wrap('broth', ncol = 2) +
  guides(fill = guide_legend(title="Sequence type/ Strain")) +
  labs(x = "Pre-enrichment duration (hrs)", y = "Relative Abundance (%)") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), text = element_text(size=20))

plot(Stool2Klebs_strains_plot)
