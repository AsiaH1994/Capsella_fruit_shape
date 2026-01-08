# Silicula fruit shape analysis
# Asia Hightower and Emily Josehps
# September 2025 - December 2025 

rm(list=ls())

#### packages #### 
library(tidyverse)
library(ggplot2)
library(dplyr) 
library(ggplot2)
library(emmeans)
library(AICcmodavg)
library(ggpubr)
library(mapdata)
library(maps)
library(rnaturalearth)
library(rnaturalearthdata)
library(Momocs)
library(ggsignif)
library(car)
library(GGally)
library(jtools)
library(broom)
library(patchwork)
library(wesanderson)
library(heplots)
library(factoextra)
#install.packages("geodata")
library(geodata)
library(geosphere)
library(reshape2)
library(mapview)
library(sf)

setwd("/Users/asiahightower/Documents/fruit_shape_paper_data_010725")

#### Reading in data ####

# fruit_data = 83 invidiuals collected as part of dataset, including seed number, weight, and lateral shoots
fruit_data <- read.csv('temp_fruit_seed_project_010726.csv')

# count data of 83 invidiuals 
indv_count <- read.csv("genotype_temperature_count.csv")

# shape_angles = shape data including PC1/PC2, shape descriptors, and shape measurements
shape_angles <- read.csv('shape_angles_01072026.csv')

# Growth data 
growth_clean <- read.csv('growth_clean_010726.csv')

#### count and survivorship analysis ####
count_p0 <- ggplot(indv_count, aes(x = factor(temperature), y = count)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_text(aes(label = count), 
            position = position_dodge(width = 0.9), 
            vjust = -0.3, size = 8) +
  geom_signif(comparisons = list(c("16", "20"), c("16", "30"), c("20", "30")),
              map_signif_level = TRUE,step_increase = 0.1, textsize = 6) +
  facet_wrap(~paper_genotype)

png(filename = 'survivorship_count_withingeno_010726.png', res = 300, width= 3800, height = 4000)
count_p0 + xlab('Temperature') + ylab("Count of individuals at end of experiment") +
  theme(panel.grid.major = element_line(linewidth = 0.3, linetype = 'solid',
                                        colour = "lightgray"), 
        panel.grid.minor = element_line(linewidth = 0.3, linetype = 'solid',
                                        colour = "lightgray"), 
        panel.background = element_rect(fill = "white", colour = "black",
                                        linewidth = 0.3, linetype = "solid"), 
        axis.line = element_line(colour = "black"), 
        plot.title = element_text(size = 14), 
        legend.position = "right",
        legend.title = element_text(size = 30),
        legend.text = element_text(size = 30), 
        axis.title = element_text(size = 30), 
        axis.text = element_text(size = 30),
        strip.text = element_text(size=30),
        axis.text.x = element_text(angle = 90,vjust = 0.5))
dev.off()

count_p1 <- ggplot(indv_count, aes(x = paper_genotype, y = count, fill = paper_genotype, colour = genotype)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_text(aes(label = count), 
            position = position_dodge(width = 0.9), 
            vjust = -0.3, size = 8) +
  facet_wrap(~temperature)

cfill <- c("white", "white", "white", "white", "white")

png(filename = 'survivorship_count_010726.png', res = 300, width= 3800, height = 4000)
count_p1 + ylim(0,10) + xlab('Genotype') + ylab("Count of individuals at end of experiment") +
  guides(size = "none") +
  guides(colour = "none") +
  guides(fill = "none") + 
  scale_fill_manual(values = cfill) + 
  scale_colour_manual(values = c("black", "black", "black", "black", "black")) + 
  theme(panel.grid.major = element_line(linewidth = 0.3, linetype = 'solid',
                                        colour = "lightgray"), 
        panel.grid.minor = element_line(linewidth = 0.3, linetype = 'solid',
                                        colour = "lightgray"), 
        panel.background = element_rect(fill = "white", colour = "black",
                                        linewidth = 0.3, linetype = "solid"), 
        axis.line = element_line(colour = "black"), 
        plot.title = element_text(size = 14), 
        legend.position = "right",
        legend.title = element_text(size = 30),
        legend.text = element_text(size = 30), 
        axis.title = element_text(size = 30), 
        axis.text = element_text(size = 30),
        strip.text = element_text(size=30),
        axis.text.x = element_text(angle = 90,vjust = 0.5))
dev.off()

lm_count <- lm(count ~ paper_genotype + temperature + paper_genotype:temperature, data = indv_count)
summary(lm_count)
aov_count <- Anova(lm_count)

summ(lm_count)

tmodel <- tidy(lm_count)
tmodel2 <- tidy(aov_count)

lm_sur <- lm(survived ~ paper_genotype + temperature + paper_genotype:temperature, data = indv_count)
summary(lm_sur)
aov_sur <- Anova(lm_sur)

sur_p1 <- ggplot(indv_count, aes(x = paper_genotype, y = survived, fill = paper_genotype, colour = paper_genotype)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_text(aes(label = survived), 
            position = position_dodge(width = 0.9), 
            vjust = -0.3, size = 8) +
  facet_wrap(~temperature)

png(filename = 'survivorship_percent_010726.png', res = 300, width= 3800, height = 4000)
sur_p1 + ylim(0,0.7) + xlab('Genotype') + ylab("Survivorship (# grew/ # planted)") +
  guides(size = "none") +
  guides(colour = "none") +
  guides(fill = "none") + 
  scale_fill_manual(values = cfill) + 
  scale_colour_manual(values = c("black", "black", "black", "black", "black")) + 
  theme(panel.grid.major = element_line(linewidth = 0.3, linetype = 'solid',
                                        colour = "lightgray"), 
        panel.grid.minor = element_line(linewidth = 0.3, linetype = 'solid',
                                        colour = "lightgray"), 
        panel.background = element_rect(fill = "white", colour = "black",
                                        linewidth = 0.3, linetype = "solid"), 
        axis.line = element_line(colour = "black"), 
        plot.title = element_text(size = 14), 
        legend.position = "right",
        legend.title = element_text(size = 30),
        legend.text = element_text(size = 30), 
        axis.title = element_text(size = 30), 
        axis.text = element_text(size = 30),
        strip.text = element_text(size=30),
        axis.text.x = element_text(angle = 90,vjust = 0.5))
dev.off()

#### Analysis of reproductive traits ####

fruit_whole <- fruit_data %>%
  group_by(paper_genotype, temp_condition) %>%
  mutate(mean_SN = mean(seed_per_fruit, na.rm = TRUE), 
         mean_SW = mean(seed_weight_fruit, na.rm = TRUE), 
         mean_LS = mean(lateral_shoots, na.rm = TRUE),,
         sd_SN = sd(seed_per_fruit), 
         sd_SW = sd(seed_weight_fruit), 
         sd_LS = sd(lateral_shoots)) %>%
  reframe(paper_genotype = paper_genotype,
          temp_condition = temp_condition,
          mean_SN = mean_SN,
          mean_SW = mean_SW,
          mean_LS = mean_LS,
          sd_SN = sd_SN,
          sd_SW = sd_SW,
          sd_LS = sd_LS) %>% 
  distinct(paper_genotype, temp_condition, .keep_all = TRUE)

# seed number
lm_SC <- lm(number_seeds ~ paper_genotype + temp_condition + paper_genotype:temp_condition + (1|set_number), data = fruit_data)
summary(lm_SC)
aov_SC <- Anova(lm_SC)
summ(lm_SC)

fruit_seedp1 <- ggplot(fruit_whole, aes(x = factor(temp_condition), y = mean_SN)) +
  geom_point(aes(size = 4)) + geom_line(aes(group= factor(paper_genotype))) +
  geom_errorbar(aes(ymin=mean_SN-sd_SN, ymax=mean_SN+sd_SN), width=0.1) + 
  geom_signif(comparisons = list(c("16", "20"), c("16", "30"), c("20", "30")),
              map_signif_level = TRUE,step_increase = 0.1, textsize = 6)

png(filename = 'fruit_seed_reaction_norm_010726.png', res = 300, width= 2600, height = 3200)
fruit_seedp1 + xlab('Temperature Condition') + ylab('Average seed number per silicula') +
  guides(size = "none") + 
  theme(panel.grid.major = element_line(linewidth = 0.3, linetype = 'solid',
                                        colour = "lightgray"), 
        panel.grid.minor = element_line(linewidth = 0.3, linetype = 'solid',
                                        colour = "lightgray"), 
        panel.background = element_rect(fill = "white", colour = "black",
                                        linewidth = 0.3, linetype = "solid"), 
        axis.line = element_line(colour = "black"), 
        plot.title = element_text(size = 30), 
        legend.text = element_text(size = 9), 
        axis.title = element_text(size = 30), 
        axis.text = element_text(size = 30),
        strip.text = element_text(size=25))
dev.off()

# seed weight 

lm_SW <- lm(see_weight_mg ~ paper_genotype + temp_condition + paper_genotype:temp_condition + (1|round), data = fruit_data)
summary(lm_SW)
aov_SW <- Anova(lm_SW)
summ(lm_SW)

fruit_weightp1 <- ggplot(fruit_whole, aes(x = factor(temp_condition), y = mean_SW)) +
  geom_point(aes(size = 4)) + geom_line(aes(group= factor(paper_genotype))) +
  geom_errorbar(aes(ymin=mean_SW-sd_SW, ymax=mean_SW+sd_SW), width=0.1) + 
  geom_signif(test =  ,comparisons = list(c("16", "20"), c("16", "30"), c("20", "30")),
              map_signif_level = TRUE,step_increase = 0.1, textsize = 6)

png(filename = 'seed_weight_reaction_norm_010726.png', res = 300, width= 2600, height = 3200)
fruit_weightp1 + xlab('Temperature Condition') + ylab('Average seed weight per fruit (mg)') +
  guides(size = "none") + 
  theme(panel.grid.major = element_line(linewidth = 0.3, linetype = 'solid',
                                        colour = "lightgray"), 
        panel.grid.minor = element_line(linewidth = 0.3, linetype = 'solid',
                                        colour = "lightgray"), 
        panel.background = element_rect(fill = "white", colour = "black",
                                        linewidth = 0.3, linetype = "solid"), 
        axis.line = element_line(colour = "black"), 
        plot.title = element_text(size = 30), 
        legend.text = element_text(size = 9), 
        axis.title = element_text(size = 30), 
        axis.text = element_text(size = 30),
        strip.text = element_text(size=25))
dev.off()

# lateral shoots
lm_LS <- lm(lateral_shoots ~ paper_genotype + temp_condition + paper_genotype:temp_condition + (1|round), data = fruit_data)
summary(lm_LS)
summ(lm_LS)
aov_LS <- Anova(lm_LS)

lateral_shootp1 <- ggplot(fruit_whole, aes(x = factor(temp_condition), y = mean_LS)) +
  geom_point(aes(size = 4)) + geom_line(aes(group= factor(paper_genotype))) +
  geom_errorbar(aes(ymin=mean_LS-sd_LS, ymax=mean_LS+sd_LS), width=0.1) + 
  geom_signif(test =  ,comparisons = list(c("16", "20"), c("16", "30"), c("20", "30")),
              map_signif_level = TRUE,step_increase = 0.1, textsize = 6)

png(filename = 'lateral_shoot_reaction_norm_010726.png', res = 300, width= 2600, height = 3200)
lateral_shootp1 + xlab('Temperature Condition') + ylab('Count of Lateral Shoots') +
  guides(size = "none") + 
  theme(panel.grid.major = element_line(linewidth = 0.3, linetype = 'solid',
                                        colour = "lightgray"), 
        panel.grid.minor = element_line(linewidth = 0.3, linetype = 'solid',
                                        colour = "lightgray"), 
        panel.background = element_rect(fill = "white", colour = "black",
                                        linewidth = 0.3, linetype = "solid"), 
        axis.line = element_line(colour = "black"), 
        plot.title = element_text(size = 30), 
        legend.text = element_text(size = 9), 
        axis.title = element_text(size = 30), 
        axis.text = element_text(size = 30),
        strip.text = element_text(size=25))
dev.off()

# correlation between reproductive traits

cor1 <- cor(fruit_data[, c('number_seeds', 'see_weight_mg', 'lateral_shoots')], method = "pearson")
cor2 <- cor(fruit_data[, c('number_seeds', 'see_weight_mg', 'lateral_shoots')], method = "spearman")

cor_plot1 <- ggpairs(fruit_data,                 
                     columns = c('seed_per_fruit', 'seed_weight_fruit', 'lateral_shoots'),        
                     aes(color = paper_genotype,
                         alpha = 0.5),
                     columnLabels = c("SN", "SW (mg)", "LS"))

png(filename = 'cor_plot_010726.png', res = 300, width= 4500, height = 4800)
cor_plot1 + scale_colour_manual(values = wes_palette("AsteroidCity1")) +
  scale_fill_manual(values = wes_palette("AsteroidCity1")) + 
  theme(panel.grid.major = element_line(linewidth = 0.3, linetype = 'solid',
                                        colour = "lightgray"), 
        panel.grid.minor = element_line(linewidth = 0.3, linetype = 'solid',
                                        colour = "lightgray"), 
        panel.background = element_rect(fill = "white", colour = "black",
                                        linewidth = 0.3, linetype = "solid"), 
        axis.line = element_line(colour = "black"), 
        plot.title = element_text(size = 30), 
        axis.title = element_text(size = 30), 
        axis.text = element_text(size = 30),
        strip.text = element_text(size=30))
dev.off()

SW_SW <- ggplot(fruit_data, aes(x = seed_per_fruit, y = seed_weight_fruit, col = paper_genotype)) +
  geom_point(aes(size = 4)) + geom_smooth(method = "lm", formula = y ~ x, se = FALSE, color = "blue")

png(filename = 'SW_SN_010726.png', res = 300, width= 3200, height = 3200)
SW_SW + xlab('Seed number per fruit') + ylab('Seed weight per fruit (mg)') +
  guides(size = "none") + 
  stat_cor(method = "pearson",
           aes(label = paste(after_stat(rr.label), after_stat(p.label), sep = "*`,`~")), 
           label.x = 5,
           position = "jitter",
           size = 7.5) +
  scale_colour_manual(values = wes_palette("AsteroidCity1")) + 
  theme(panel.grid.major = element_line(linewidth = 0.3, linetype = 'solid',
                                        colour = "lightgray"), 
        panel.grid.minor = element_line(linewidth = 0.3, linetype = 'solid',
                                        colour = "lightgray"), 
        panel.background = element_rect(fill = "white", colour = "black",
                                        linewidth = 0.3, linetype = "solid"), 
        axis.line = element_line(colour = "black"), 
        plot.title = element_text(size = 30), 
        legend.text = element_text(size = 10), 
        axis.title = element_text(size = 30), 
        axis.text = element_text(size = 30))
dev.off()

g2 <- filter(fruit_data, genotype == "2.16")
g20 <- filter(fruit_data, genotype == "20.295")
gM6 <- filter(fruit_data, genotype == "Magnolia6")
g7 <- filter(fruit_data, genotype == "7.99")
g14 <- filter(fruit_data, genotype == "14.2")

g2p1 <- ggplot(g2, aes(x = seed_per_fruit, y = lateral_shoots)) +
  geom_point(aes(size = 4)) + 
  geom_smooth(method = "lm", formula = y ~ x + I(x^2), se = FALSE, color = "blue") +
  stat_cor(
    aes(label = paste(after_stat(rr.label), after_stat(p.label), sep = "*`,`~")), 
    label.x = 5, 
    size = 7.5)

g2p2 <- g2p1 + xlim(0,35) + ylim(0,14) + 
  xlab('Seed number per fruit') + ylab('Number of lateral shoots') +
  guides(size = "none") +
  ggtitle("NY_63") + 
  theme(panel.grid.major = element_line(linewidth = 0.3, linetype = 'solid',
                                        colour = "lightgray"), 
        panel.grid.minor = element_line(linewidth = 0.3, linetype = 'solid',
                                        colour = "lightgray"), 
        panel.background = element_rect(fill = "white", colour = "black",
                                        linewidth = 0.3, linetype = "solid"), 
        axis.line = element_line(colour = "black"), 
        plot.title = element_text(size = 30), 
        legend.text = element_text(size = 10), 
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text = element_text(size = 30))

g20p1 <- ggplot(g20, aes(x = seed_per_fruit, y = lateral_shoots)) +
  geom_point(aes(size = 4)) + 
  geom_smooth(method = "lm", formula = y ~ x + I(x^2), se = FALSE, color = "blue") +
  stat_cor(
    aes(label = paste(after_stat(rr.label), after_stat(p.label), sep = "*`,`~")), 
    label.x = 5,
    label.y = 12,
    size = 7.5)

g20p2 <- g20p1 + xlim(0,35) + ylim(0,14) + 
  xlab('Seed number per fruit') + ylab('Number of lateral shoots') +
  guides(size = "none") +
  ggtitle("NY_30") + 
  theme(panel.grid.major = element_line(linewidth = 0.3, linetype = 'solid',
                                        colour = "lightgray"), 
        panel.grid.minor = element_line(linewidth = 0.3, linetype = 'solid',
                                        colour = "lightgray"), 
        panel.background = element_rect(fill = "white", colour = "black",
                                        linewidth = 0.3, linetype = "solid"), 
        axis.line = element_line(colour = "black"), 
        plot.title = element_text(size = 30), 
        legend.text = element_text(size = 10), 
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text = element_text(size = 30))

gm6p1 <- ggplot(gM6, aes(x = seed_per_fruit, y = lateral_shoots)) +
  geom_point(aes(size = 4)) + 
  geom_smooth(method = "lm", formula = y ~ x + I(x^2), se = FALSE, color = "blue") +
  stat_cor(
    aes(label = paste(after_stat(rr.label), after_stat(p.label), sep = "*`,`~")), 
    label.x = 5,
    size = 7.5)

gm6p2 <- gm6p1 + xlim(0,35) + ylim(0,14) +  
  xlab('Seed number per fruit') + ylab('Number of lateral shoots') +
  guides(size = "none") +
  ggtitle("MI_30") +
  theme(panel.grid.major = element_line(linewidth = 0.3, linetype = 'solid',
                                        colour = "lightgray"), 
        panel.grid.minor = element_line(linewidth = 0.3, linetype = 'solid',
                                        colour = "lightgray"), 
        panel.background = element_rect(fill = "white", colour = "black",
                                        linewidth = 0.3, linetype = "solid"), 
        axis.line = element_line(colour = "black"), 
        plot.title = element_text(size = 30), 
        legend.text = element_text(size = 10), 
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text = element_text(size = 30))

g7p1 <- ggplot(g7, aes(x = seed_per_fruit, y = lateral_shoots)) +
  geom_point(aes(size = 4)) + 
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE, color = "blue") +
  stat_cor(
    aes(label = paste(after_stat(rr.label), after_stat(p.label), sep = "*`,`~")), 
    label.x = 5,
    size = 7.5)

g7p2 <- g7p1 + xlim(0,35) + ylim(0,14) + 
  xlab('Seed number per fruit') + ylab('Number of lateral shoots') +
  guides(size = "none") +
  ggtitle("NY_56") + 
  theme(panel.grid.major = element_line(linewidth = 0.3, linetype = 'solid',
                                        colour = "lightgray"), 
        panel.grid.minor = element_line(linewidth = 0.3, linetype = 'solid',
                                        colour = "lightgray"), 
        panel.background = element_rect(fill = "white", colour = "black",
                                        linewidth = 0.3, linetype = "solid"), 
        axis.line = element_line(colour = "black"), 
        plot.title = element_text(size = 30), 
        legend.text = element_text(size = 10), 
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text = element_text(size = 30))

g14p1 <- ggplot(g14, aes(x = seed_per_fruit, y = lateral_shoots)) +
  geom_point(aes(size = 4)) + 
  geom_smooth(method = "lm", formula = y ~ x + I(x^2), se = FALSE, color = "blue") +
  stat_cor(
    aes(label = paste(after_stat(rr.label), after_stat(p.label), sep = "*`,`~")), 
    label.x = 5,
    size = 7.5)


g14p2 <- g14p1 + xlim(0,35) + ylim(0,14) + 
  xlab('Seed number per fruit') + 
  ylab('Number of lateral shoots') +
  guides(size = "none") +
  ggtitle("NY_09") + 
  theme(panel.grid.major = element_line(linewidth = 0.3, linetype = 'solid',
                                        colour = "lightgray"), 
        panel.grid.minor = element_line(linewidth = 0.3, linetype = 'solid',
                                        colour = "lightgray"), 
        panel.background = element_rect(fill = "white", colour = "black",
                                        linewidth = 0.3, linetype = "solid"), 
        axis.line = element_line(colour = "black"), 
        plot.title = element_text(size = 30), 
        legend.text = element_text(size = 10), 
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text = element_text(size = 30))

combo1 <- ggarrange(g2p2, g20p2, gm6p2, g7p2, g14p2, ncol = 2, nrow = 3)

png(filename = 'LS_SN_010726_sign.png', res = 300, width= 3800, height = 4200)
annotate_figure(combo1,
                bottom = text_grob("Number of lateral shoots", 
                                   vjust = 0, size = 30),
                left = text_grob("Seed number per fruit", 
                                 rot = 90, size = 30))
dev.off()

#### Analysis of vegetative traits ####

# leaf number
lm_lf1 <- lm(num_leaves ~ paper_genotype + condition + paper_genotype:condition + (1|round), data = growth_clean)
summary(lm_lf1)
summ(lm_lf1)
aov_lf1 <- Anova(lm_lf1)

lfnp1 <- ggplot(growth_clean, aes(x = factor(condition), y = num_leaves)) + 
  geom_violin(position = position_dodge(0.9)) + 
  geom_boxplot(width = 0.1, position = position_dodge(0.9)) +
  geom_signif(test =  ,comparisons = list(c("16", "20"), c("16","30"), c("20", "30")),
              map_signif_level = TRUE,step_increase = 0.1, textsize = 6)

png(filename = 'leaf_number_temp_010726.png', res = 300, width= 3800, height = 4200)
lfnp1 + xlab('Temperature') + ylab('Number of leaves per individual') +
  theme(panel.grid.major = element_line(linewidth = 0.3, linetype = 'solid',
                                        colour = "lightgray"), 
        panel.grid.minor = element_line(linewidth = 0.3, linetype = 'solid',
                                        colour = "lightgray"), 
        panel.background = element_rect(fill = "white", colour = "black",
                                        linewidth = 0.3, linetype = "solid"), 
        axis.line = element_line(colour = "black"), 
        plot.title = element_text(size = 30), 
        legend.text = element_text(size = 10), 
        axis.title = element_text(size = 30),
        axis.text = element_text(size = 30),
        strip.text = element_text(size=30),
        axis.text.x = element_text(angle = 90,vjust = 0.5))
dev.off()

lfnp2 <- ggplot(growth_clean, aes(x = factor(paper_genotype), y = num_leaves)) + 
  geom_violin(position = position_dodge(0.9)) + 
  geom_boxplot(width = 0.1, position = position_dodge(0.9)) +
  geom_signif(comparisons = list(c("NY_09", "NY_63"), 
                                 c("NY_09", "NY_30"), 
                                 c("NY_09", "NY_56"),
                                 c("NY_09", "MI_30"),
                                 c("NY_63", "NY_30"),
                                 c("NY_63", "NY_56"),
                                 c("NY_63", "MI_30"),
                                 c("NY_30", "NY_56"),
                                 c("NY_30", "MI_30"),
                                 c("NY_56", 'MI_30')),
              map_signif_level = TRUE, step_increase = 0.1, textsize = 6) +
  facet_wrap(~condition)

png(filename = 'leaf_number_by_paper_genotype_010726.png', res = 300, width= 3800, height = 4200)
lfnp2 + ylim(0,125) + xlab('Genotype') + ylab('Number of leaves per individual') + 
  theme(panel.grid.major = element_line(linewidth = 0.3, linetype = 'solid',
                                        colour = "lightgray"), 
        panel.grid.minor = element_line(linewidth = 0.3, linetype = 'solid',
                                        colour = "lightgray"), 
        panel.background = element_rect(fill = "white", colour = "black",
                                        linewidth = 0.3, linetype = "solid"), 
        axis.line = element_line(colour = "black"), 
        plot.title = element_text(size = 30), 
        legend.text = element_text(size = 10), 
        axis.title = element_text(size = 30),
        axis.text = element_text(size = 30),
        strip.text = element_text(size=30),
        axis.text.x = element_text(angle = 90,vjust = 0.5))
dev.off()

# plant width 
lm_pw1 <- lm(plant_width ~ paper_genotype + condition + paper_genotype:condition + (1|round), data = growth_clean)
summary(lm_pw1)
summ(lm_pw1)
aov_lf1 <- Anova(lm_pw1)

pw1 <- ggplot(growth_clean, aes(x = factor(condition), y = plant_width)) + 
  geom_violin(position = position_dodge(0.9)) + 
  geom_boxplot(width = 0.1, position = position_dodge(0.9)) +
  geom_signif(test =  ,comparisons = list(c("16", "20"), c("16","30"), c("20", "30")),
              map_signif_level = TRUE,step_increase = 0.1, textsize = 6)

png(filename = 'plant_width_by_temp_010726.png', res = 300, width= 3800, height = 4200)
pw1 + xlab('Temperature') + ylab('Plant Width (cm) per individual') +
  theme(panel.grid.major = element_line(linewidth = 0.3, linetype = 'solid',
                                        colour = "lightgray"), 
        panel.grid.minor = element_line(linewidth = 0.3, linetype = 'solid',
                                        colour = "lightgray"), 
        panel.background = element_rect(fill = "white", colour = "black",
                                        linewidth = 0.3, linetype = "solid"), 
        axis.line = element_line(colour = "black"), 
        plot.title = element_text(size = 30), 
        legend.text = element_text(size = 10), 
        axis.title = element_text(size = 30),
        axis.text = element_text(size = 30),
        strip.text = element_text(size=30),
        axis.text.x = element_text(angle = 90,vjust = 0.5))
dev.off()

pw2 <- ggplot(growth_clean, aes(x = factor(paper_genotype), y = plant_width)) + 
  geom_violin(position = position_dodge(0.9)) + 
  geom_boxplot(width = 0.1, position = position_dodge(0.9)) +
  geom_signif(comparisons = list(c("NY_09", "NY_63"), 
                                 c("NY_09", "NY_30"), 
                                 c("NY_09", "NY_56"),
                                 c("NY_09", "MI_30"),
                                 c("NY_63", "NY_30"),
                                 c("NY_63", "NY_56"),
                                 c("NY_63", "MI_30"),
                                 c("NY_30", "NY_56"),
                                 c("NY_30", "MI_30"),
                                 c("NY_56", 'MI_30')),
              map_signif_level = TRUE,step_increase = 0.1, textsize = 6) +
  facet_wrap(~condition)

png(filename = 'plant_width_by_paper_genotype_010726.png', res = 300, width= 3800, height = 4200)
pw2 + ylim(0,65) + xlab('Genotype') + ylab('Plant Width (cm) per individual') + 
  theme(panel.grid.major = element_line(linewidth = 0.3, linetype = 'solid',
                                        colour = "lightgray"), 
        panel.grid.minor = element_line(linewidth = 0.3, linetype = 'solid',
                                        colour = "lightgray"), 
        panel.background = element_rect(fill = "white", colour = "black",
                                        linewidth = 0.3, linetype = "solid"), 
        axis.line = element_line(colour = "black"), 
        plot.title = element_text(size = 30), 
        legend.text = element_text(size = 10), 
        axis.title = element_text(size = 30),
        axis.text = element_text(size = 30),
        strip.text = element_text(size=30),
        axis.text.x = element_text(angle = 90,vjust = 0.5))
dev.off()

#### analyzing shape data ####

# solidity

shape1 <- lm(solidity ~ length, shape_angles)
shape2 <- lm(solidity ~ width, shape_angles)
shape3 <- lm(solidity ~ ar, shape_angles)
summary(shape3)
shape4 <- lm(solidity ~ area, shape_angles)

sol_mods <- list(shape1, shape2, shape3, shape4)
sol.names <- c("shape1", "shape2", "shape3", "shape4")
sol_comp <- aictab(cand.set = sol_mods, modnames = sol.names)

col_temp = c("blue", "black", "red")

ar_solp1 <-ggplot(shape_angles, aes(x = ar, y = solidity, col = factor(dataset))) + 
  geom_point() + geom_smooth(method = "lm", formula = y ~ x + I(x^2), se = FALSE, color = "blue")

png(filename = "ar_by_solidity_010726.png", res = 300, width= 3000, height = 3400)
ar_solp1 + 
  xlab("Aspect Ratio") +
  ylab("Solidity") +
  guides(color = guide_legend(title = "Temperature")) + 
  scale_colour_manual(values = col_temp) + 
  theme(panel.grid.major = element_line(linewidth = 0.3, linetype = 'solid',
                                        colour = "lightgray"), 
        panel.grid.minor = element_line(linewidth = 0.3, linetype = 'solid',
                                        colour = "lightgray"), 
        panel.background = element_rect(fill = "white", colour = "black",
                                        linewidth = 0.3, linetype = "solid"), 
        axis.line = element_line(colour = "black"), 
        plot.title = element_text(size = 14), 
        legend.position = "right",
        legend.title = element_text(size = 30),
        legend.text = element_text(size = 30), 
        axis.title = element_text(size = 30), 
        axis.text = element_text(size = 30),
        strip.text = element_text(size=30),
        axis.text.x = element_text(angle = 90,vjust = 0.5))
dev.off()

# asymmetry
shape5 <- lm(asymmetry ~ length, shape_angles)
shape6 <- lm(asymmetry ~ width, shape_angles)
shape7 <- lm(asymmetry ~ ar, shape_angles)
summary(shape7)
shape8 <- lm(asymmetry ~ area, shape_angles)

asy_mods <- list(shape5, shape6, shape7, shape8)
asy_names <- c("shape5", "shape6", "shape7", "shape8")
asy_comp <- aictab(cand.set = asy_mods, modnames = asy_names)

ar_asyp1 <-ggplot(shape_angles, aes(x = ar, y = asymmetry, col = factor(dataset))) + 
  geom_point() + geom_smooth(method = "lm", formula = y ~ x, se = FALSE, color = "blue")

png(filename = "ar_by_asymmmetry_010726.png", res = 300, width= 3000, height = 3400)
ar_asyp1 + 
  xlab("Aspect Ratio") +
  ylab("Asymmetry") +
  guides(color = guide_legend(title = "Temperature")) + 
  scale_colour_manual(values = col_temp) + 
  theme(panel.grid.major = element_line(linewidth = 0.3, linetype = 'solid',
                                        colour = "lightgray"), 
        panel.grid.minor = element_line(linewidth = 0.3, linetype = 'solid',
                                        colour = "lightgray"), 
        panel.background = element_rect(fill = "white", colour = "black",
                                        linewidth = 0.3, linetype = "solid"), 
        axis.line = element_line(colour = "black"), 
        plot.title = element_text(size = 14), 
        legend.position = "right",
        legend.title = element_text(size = 30),
        legend.text = element_text(size = 30), 
        axis.title = element_text(size = 30), 
        axis.text = element_text(size = 30),
        strip.text = element_text(size=30),
        axis.text.x = element_text(angle = 90,vjust = 0.5))
dev.off()

# Circularity

shape9 <- lm(circ ~ length, shape_angles)
shape10 <- lm(circ ~ width, shape_angles)
shape11 <- lm(circ ~ ar, shape_angles)
summary(shape11)
shape12 <- lm(circ ~ area, shape_angles)

circ_mods <- list(shape9, shape10, shape11, shape12)
circ_names <- c("shape9", "shape10", "shape11", "shape12")
circ_comp <- aictab(cand.set = circ_mods, modnames = circ_names)

ar_circp1 <-ggplot(shape_angles, aes(x = ar, y = circ, col = factor(dataset))) + 
  geom_point() + geom_smooth(method = "lm", formula = y ~ x, se = FALSE, color = "blue")

png(filename = "ar_by_circ_010726.png", res = 300, width= 3000, height = 3400)
ar_circp1 + 
  xlab("Aspect Ratio") +
  ylab("Circularity") +
  guides(color = guide_legend(title = "Temperature")) + 
  scale_colour_manual(values = col_temp) + 
  theme(panel.grid.major = element_line(linewidth = 0.3, linetype = 'solid',
                                        colour = "lightgray"), 
        panel.grid.minor = element_line(linewidth = 0.3, linetype = 'solid',
                                        colour = "lightgray"), 
        panel.background = element_rect(fill = "white", colour = "black",
                                        linewidth = 0.3, linetype = "solid"), 
        axis.line = element_line(colour = "black"), 
        plot.title = element_text(size = 14), 
        legend.position = "right",
        legend.title = element_text(size = 30),
        legend.text = element_text(size = 30), 
        axis.title = element_text(size = 30), 
        axis.text = element_text(size = 30),
        strip.text = element_text(size=30),
        axis.text.x = element_text(angle = 90,vjust = 0.5))
dev.off()

#### relationship between shape and reproductive traits ####

seed1 <- lm(seeds_per_fruit ~ length, shape_angles)
summary(seed1)
seed2 <- lm(seeds_per_fruit ~ width, shape_angles)
seed3 <- lm(seeds_per_fruit ~ area, shape_angles)
seed4 <- lm(seeds_per_fruit ~ ar, shape_angles)

seed_mods <- list(seed1, seed2, seed3, seed4)
seed_names <- c("seed1", "seed2", "seed3", "seed4")
seed_comp <- aictab(cand.set = seed_mods, modnames = seed_names)

seed_p1 <- ggplot(shape_angles, aes(x = length, y = seeds_per_fruit, col = factor(dataset))) + 
  geom_point() + geom_smooth(method = "lm", formula = y ~ x, se = FALSE, color = "blue")

png(filename = "seed_number_length_010726.png", res = 300, width= 3000, height = 3400)
seed_p1 + 
  xlab("Silicula Length (cm)") +
  ylab("Number of seeds per fruit") +
  guides(color = guide_legend(title = "Temperature")) + 
  scale_colour_manual(values = col_temp) + 
  theme(panel.grid.major = element_line(linewidth = 0.3, linetype = 'solid',
                                        colour = "lightgray"), 
        panel.grid.minor = element_line(linewidth = 0.3, linetype = 'solid',
                                        colour = "lightgray"), 
        panel.background = element_rect(fill = "white", colour = "black",
                                        linewidth = 0.3, linetype = "solid"), 
        axis.line = element_line(colour = "black"), 
        plot.title = element_text(size = 14), 
        legend.position = "right",
        legend.title = element_text(size = 30),
        legend.text = element_text(size = 30), 
        axis.title = element_text(size = 30), 
        axis.text = element_text(size = 30),
        strip.text = element_text(size=30),
        axis.text.x = element_text(angle = 90,vjust = 0.5))
dev.off()

seed5 <- lm(seeds_per_fruit ~ solidity + ar, shape_angles)
seed6 <- lm(seeds_per_fruit ~ asymmetry + ar, shape_angles)
seed7 <- lm(seeds_per_fruit ~ circ + ar, shape_angles)
summary(seed7)

seed_mods2 <- list(seed5, seed6, seed7)
seed_names2 <- c("seed5", "seed6", "seed7")
seed_comp2 <- aictab(cand.set = seed_mods2, modnames = seed_names2)

seed_p2 <- ggplot(shape_angles, aes(x = ar, y = circ, col = seeds_per_fruit)) + 
  geom_point() + geom_smooth(method = "lm", formula = y ~ x, se = FALSE, color = "blue")

png(filename = "seed_number_circ_ar_010726.png", res = 300, width= 3000, height = 3400)
seed_p2 + 
  xlab("Aspect Ratio") +
  ylab("Circularity") +
  guides(color = guide_legend(title = "Seed Number")) + 
  theme(panel.grid.major = element_line(linewidth = 0.3, linetype = 'solid',
                                        colour = "lightgray"), 
        panel.grid.minor = element_line(linewidth = 0.3, linetype = 'solid',
                                        colour = "lightgray"), 
        panel.background = element_rect(fill = "white", colour = "black",
                                        linewidth = 0.3, linetype = "solid"), 
        axis.line = element_line(colour = "black"), 
        plot.title = element_text(size = 14), 
        legend.position = "right",
        legend.title = element_text(size = 30),
        legend.text = element_text(size = 30), 
        axis.title = element_text(size = 30), 
        axis.text = element_text(size = 30),
        strip.text = element_text(size=30),
        axis.text.x = element_text(angle = 90,vjust = 0.5))
dev.off()

#### shape measurements by paper_genotype ####

mean_data <- shape_angles %>%
  group_by(paper_genotype) %>%
  mutate(ml = mean(length),
         mw = mean(width),
         ma = mean(area), 
         mar = mean(ar),
         mc = mean(circ),
         ms = mean(solidity),
         masy = mean(asymmetry),
         sd_l = sd(length),
         sd_w = sd(width),
         sd_a = sd(area),
         sd_ar = sd(ar),
         sd_c = sd(circ),
         sd_sol = sd(solidity),
         sd_asy = sd(asymmetry)) %>%
  reframe(paper_genotype = paper_genotype,
          ml = ml,
          mw = mw,
          ma = ma,
          mar = mar,
          mc = mc,
          ms = ms, 
          masy = masy,
          sd_l = sd_l,
          sd_w = sd_w, 
          sd_a = sd_a,
          sd_ar = sd_ar,
          as_c = sd_c,
          sd_sol = sd_sol,
          sd_asy = sd_asy) %>%
  distinct(paper_genotype, .keep_all = TRUE)

geno2 <- lm(length ~ paper_genotype + (1|round), shape_angles)
summary(geno2)
geno3 <- lm(width ~ paper_genotype + (1|round), shape_angles)
summary(geno3)
geno4 <- lm(ar ~ paper_genotype + (1|round), shape_angles)
summary(geno4)
geno5 <- lm(area ~ paper_genotype + (1|round), shape_angles)
summary(geno5)

geno_lengthp1 <- ggplot(shape_angles, aes(x = factor(paper_genotype), y = length)) + 
  geom_violin(position = position_dodge(0.9)) + 
  geom_boxplot(width = 0.1, position = position_dodge(0.9))

g1 <- geno_lengthp1 + 
  ylim(0.5, 1.8) +
  xlab("Genotype") +
  ylab("Length (cm)") + 
  theme(panel.grid.major = element_line(linewidth = 0.3, linetype = 'solid',
                                        colour = "lightgray"), 
        panel.grid.minor = element_line(linewidth = 0.3, linetype = 'solid',
                                        colour = "lightgray"), 
        panel.background = element_rect(fill = "white", colour = "black",
                                        linewidth = 0.3, linetype = "solid"), 
        axis.line = element_line(colour = "black"), 
        plot.title = element_text(size = 14), 
        legend.position = "right",
        legend.title = element_text(size = 30),
        legend.text = element_text(size = 30), 
        axis.title = element_text(size = 30), 
        axis.text = element_text(size = 30),
        strip.text = element_text(size=30),
        axis.text.x = element_text(angle = 90,vjust = 0.5))

png(filename = "paper_genotype_length_010726.png", res = 300, width= 3000, height = 3400)
g1
dev.off()

geno_arp1 <- ggplot(shape_angles, aes(x = factor(paper_genotype), y = ar)) + 
  geom_violin(position = position_dodge(0.9)) + 
  geom_boxplot(width = 0.1, position = position_dodge(0.9))

g2 <- geno_arp1 + 
  ylim(0,1.2) +
  xlab("Genotype") +
  ylab("Aspect Ratio") + 
  theme(panel.grid.major = element_line(linewidth = 0.3, linetype = 'solid',
                                        colour = "lightgray"), 
        panel.grid.minor = element_line(linewidth = 0.3, linetype = 'solid',
                                        colour = "lightgray"), 
        panel.background = element_rect(fill = "white", colour = "black",
                                        linewidth = 0.3, linetype = "solid"), 
        axis.line = element_line(colour = "black"), 
        plot.title = element_text(size = 14), 
        legend.position = "right",
        legend.title = element_text(size = 30),
        legend.text = element_text(size = 30), 
        axis.title = element_text(size = 30), 
        axis.text = element_text(size = 30),
        strip.text = element_text(size=30),
        axis.text.x = element_text(angle = 90,vjust = 0.5))

png(filename = "paper_genotype_ar_010726.png", res = 300, width= 3000, height = 3400)
g2
dev.off()

geno_widthp1 <- ggplot(shape_angles, aes(x = factor(paper_genotype), y = width)) + 
  geom_violin(position = position_dodge(0.9)) + 
  geom_boxplot(width = 0.1, position = position_dodge(0.9))

g3 <- geno_widthp1 +
  ylim(0,1.2) +
  xlab("Genotype") +
  ylab("Width (cm)") + 
  theme(panel.grid.major = element_line(linewidth = 0.3, linetype = 'solid',
                                        colour = "lightgray"), 
        panel.grid.minor = element_line(linewidth = 0.3, linetype = 'solid',
                                        colour = "lightgray"), 
        panel.background = element_rect(fill = "white", colour = "black",
                                        linewidth = 0.3, linetype = "solid"), 
        axis.line = element_line(colour = "black"), 
        plot.title = element_text(size = 14), 
        legend.position = "right",
        legend.title = element_text(size = 30),
        legend.text = element_text(size = 30), 
        axis.title = element_text(size = 30), 
        axis.text = element_text(size = 30),
        strip.text = element_text(size=30),
        axis.text.x = element_text(angle = 90,vjust = 0.5))

png(filename = "paper_genotype_width_010726.png", res = 300, width= 3000, height = 3400)
g3
dev.off()

geno_areap1 <- ggplot(shape_angles, aes(x = factor(paper_genotype), y = area)) + 
  geom_violin(position = position_dodge(0.9)) + 
  geom_boxplot(width = 0.1, position = position_dodge(0.9))

png(filename = "paper_genotype_area_010726.png", res = 300, width= 3000, height = 3400)
g4 <- geno_areap1 + 
  ylim(0,1.2) +
  xlab("Genotype") +
  ylab("Area (cm)") + 
  theme(panel.grid.major = element_line(linewidth = 0.3, linetype = 'solid',
                                        colour = "lightgray"), 
        panel.grid.minor = element_line(linewidth = 0.3, linetype = 'solid',
                                        colour = "lightgray"), 
        panel.background = element_rect(fill = "white", colour = "black",
                                        linewidth = 0.3, linetype = "solid"), 
        axis.line = element_line(colour = "black"), 
        plot.title = element_text(size = 14), 
        legend.position = "right",
        legend.title = element_text(size = 30),
        legend.text = element_text(size = 30), 
        axis.title = element_text(size = 30), 
        axis.text = element_text(size = 30),
        strip.text = element_text(size=30),
        axis.text.x = element_text(angle = 90,vjust = 0.5))

png(filename = "paper_genotype_area_010726.png", res = 300, width= 3000, height = 3400)
g4
dev.off()

png(filename = "paper_genotype_all_010726.png", res = 300, width= 5000, height = 5400)
ggarrange(g1, g2, g3, g4, ncol = 2, nrow = 2, hjust = 0, vjust = 0)
dev.off()

#### shape descriptors by paper_genotype ####

lm1 <- lm(solidity ~ paper_genotype + dataset + paper_genotype:factor(dataset), shape_angles)
summary(lm1)

sol_p2 <- ggplot(shape_angles, aes(x = factor(paper_genotype), y = solidity)) + 
  geom_violin(position = position_dodge(0.9)) + 
  geom_boxplot(width = 0.1, position = position_dodge(0.9)) 

png(filename = "paper_genotype_solidity_010726.png", res = 300, width= 3500, height = 3800)
sol_p2 +
  xlab("Genotype") +
  ylab("Solidity") +
  theme(panel.grid.major = element_line(linewidth = 0.3, linetype = 'solid',
                                        colour = "lightgray"), 
        panel.grid.minor = element_line(linewidth = 0.3, linetype = 'solid',
                                        colour = "lightgray"), 
        panel.background = element_rect(fill = "white", colour = "black",
                                        linewidth = 0.3, linetype = "solid"), 
        axis.line = element_line(colour = "black"), 
        plot.title = element_text(size = 14), 
        legend.position = "right",
        legend.title = element_text(size = 30),
        legend.text = element_text(size = 30), 
        axis.title = element_text(size = 30), 
        axis.text = element_text(size = 30),
        strip.text = element_text(size=30),
        axis.text.x = element_text(angle = 90,vjust = 0.5))
dev.off()

lm2 <- lm(asymmetry ~ paper_genotype + dataset + paper_genotype:dataset, shape_angles)
summary(lm2)

asy_p2 <- ggplot(shape_angles, aes(x = factor(paper_genotype), y = asymmetry)) + 
  geom_violin(position = position_dodge(0.9)) + 
  geom_boxplot(width = 0.1, position = position_dodge(0.9))

png(filename = "paper_genotype_asymmetry_010726.png", res = 300, width= 3500, height = 3800)
asy_p2 + 
  xlab("Genotype") +
  ylab("Asymmetry") + 
  guides(fill = "none") +
  theme(panel.grid.major = element_line(linewidth = 0.3, linetype = 'solid',
                                        colour = "lightgray"), 
        panel.grid.minor = element_line(linewidth = 0.3, linetype = 'solid',
                                        colour = "lightgray"), 
        panel.background = element_rect(fill = "white", colour = "black",
                                        linewidth = 0.3, linetype = "solid"), 
        axis.line = element_line(colour = "black"), 
        plot.title = element_text(size = 14), 
        legend.position = "right",
        legend.title = element_text(size = 30),
        legend.text = element_text(size = 30), 
        axis.title = element_text(size = 30), 
        axis.text = element_text(size = 30),
        strip.text = element_text(size=30),
        axis.text.x = element_text(angle = 90,vjust = 0.5))
dev.off()

lm3 <- lm(circ ~ paper_genotype + dataset + paper_genotype:factor(dataset), shape_angles)
summary(lm3)

circ_p2 <- ggplot(shape_angles, aes(x = factor(paper_genotype), y = circ)) + 
  geom_violin(position = position_dodge(0.9)) + 
  geom_boxplot(width = 0.1, position = position_dodge(0.9))

png(filename = "paper_genotype_circ_010726.png", res = 300, width= 3500, height = 3800)
circ_p2 + 
  xlab("Genotype") +
  ylab("Circularity") + 
  guides(fill = "none") +
  theme(panel.grid.major = element_line(linewidth = 0.3, linetype = 'solid',
                                        colour = "lightgray"), 
        panel.grid.minor = element_line(linewidth = 0.3, linetype = 'solid',
                                        colour = "lightgray"), 
        panel.background = element_rect(fill = "white", colour = "black",
                                        linewidth = 0.3, linetype = "solid"), 
        axis.line = element_line(colour = "black"), 
        plot.title = element_text(size = 14), 
        legend.position = "right",
        legend.title = element_text(size = 30),
        legend.text = element_text(size = 30), 
        axis.title = element_text(size = 30), 
        axis.text = element_text(size = 30),
        strip.text = element_text(size=30),
        axis.text.x = element_text(angle = 90,vjust = 0.5))
dev.off()

#### bolting date ####

dap_temp <- ggplot(growth_clean, aes(x = factor(condition), y = relative_dap)) + 
  geom_boxplot() + geom_signif(comparisons = list(c("16", "20"), c("16", "30"), c("20", "30")),
                               map_signif_level = TRUE,step_increase = 0.1, 
                               textsize = 6)

png(filename = 'relative_dap_temp_010726.png', res = 300, width= 3800, height = 4200)
dap_temp + ylim(0,1.25) + xlab('Temperature') + ylab('Relative bolting date (DAP)') + 
  theme(panel.grid.major = element_line(linewidth = 0.3, linetype = 'solid',
                                        colour = "lightgray"), 
        panel.grid.minor = element_line(linewidth = 0.3, linetype = 'solid',
                                        colour = "lightgray"), 
        panel.background = element_rect(fill = "white", colour = "black",
                                        linewidth = 0.3, linetype = "solid"), 
        axis.line = element_line(colour = "black"), 
        plot.title = element_text(size = 30), 
        legend.text = element_text(size = 10), 
        axis.title = element_text(size = 30),
        axis.text = element_text(size = 30),
        strip.text = element_text(size=30),
        axis.text.x = element_text(angle = 90,vjust = 0.5))
dev.off()

dap_genotype <- ggplot(growth_clean, aes(x = factor(paper_genotype), y = relative_dap)) + 
  geom_boxplot() + geom_signif(comparisons = list(c("NY_09", "NY_63"), 
                                                  c("NY_09", "NY_30"), 
                                                  c("NY_09", "NY_56"),
                                                  c("NY_09", "MI_30"),
                                                  c("NY_63", "NY_30"),
                                                  c("NY_63", "NY_56"),
                                                  c("NY_63", "MI_30"),
                                                  c("NY_30", "NY_56"),
                                                  c("NY_30", "MI_30"),
                                                  c("NY_56", 'MI_30')),
                               map_signif_level = TRUE,step_increase = 0.1, textsize = 6)

png(filename = 'relative_dap_genotype__genotype_010726.png', res = 300, width= 3800, height = 4200)
dap_genotype + ylim(0,1.30) + xlab('Genotype') + ylab('Relative bolting date (DAP)') + 
  theme(panel.grid.major = element_line(linewidth = 0.3, linetype = 'solid',
                                        colour = "lightgray"), 
        panel.grid.minor = element_line(linewidth = 0.3, linetype = 'solid',
                                        colour = "lightgray"), 
        panel.background = element_rect(fill = "white", colour = "black",
                                        linewidth = 0.3, linetype = "solid"), 
        axis.line = element_line(colour = "black"), 
        plot.title = element_text(size = 30), 
        legend.text = element_text(size = 10), 
        axis.title = element_text(size = 30),
        axis.text = element_text(size = 30),
        strip.text = element_text(size=30),
        axis.text.x = element_text(angle = 90,vjust = 0.5))
dev.off()



