## Fungal guilds and root colonization ploted by horizons and forests

## Load packages and data ####
library(tidyverse)
library(nlme)
library(emmeans)
library(ggpubr)

load("data/colo3.rda")
load("data/ITS.guild.rdata")

## Colonization ####
## Data manipulation
colo3 <- droplevels(subset(colo3, !horiz %in% 'Ah')) # remove all Ah horizons
levels(colo3$myco)[levels(colo3$myco) == "ECM"] <- "EcM"
levels(colo3$myco)[levels(colo3$myco) == "mixed"] <- "Mixed"
colo3$horiz <- factor(colo3$horiz, levels = c('B', 'Ae', 'H', 'F'))

## EcM colo
colo.ecm.mod3 <- lme(pourc.ECM.colo ~ horiz*myco, random = ~ 1 | block, weights =varPower(), data = colo3)
anova(colo.ecm.mod3)
plot(colo.ecm.mod3)
tuckey.colo.ecm <- emmeans(colo.ecm.mod3, pairwise ~ horiz*myco, adjust = "tukey")
mult.colo.ecm <- CLD(tuckey.colo.ecm,alpha=0.05,Letters=letters, adjust="tukey")
L <- data.frame(horiz = c('L', 'L', 'L'), myco = c('AM','EcM','Mixed'), emmean = c(NA,NA,NA), SE = c(NA,NA,NA)) # add NA for L (no roots)
mult.colo.ecm.df <- rbind(as.data.frame(mult.colo.ecm[,1:4]), L)
mult.colo.ecm.df$myco <- factor(mult.colo.ecm.df$myco, levels = c('AM', 'Mixed', 'EcM'))

# By myco type
colo.ecm.myco.mod3 <- lme(pourc.ECM.colo ~ myco, random = ~ 1 | block, weights =varPower(), data = colo3)
tuckey.colo.myco.ecm <- emmeans(colo.ecm.myco.mod3, pairwise ~ myco, adjust = "tukey")
mult.colo.myco.ecm <- CLD(tuckey.colo.myco.ecm,alpha=0.05,Letters=letters, adjust="tukey")

## AM colo
colo.am.mod1 <- lme(pourc.AM.colo ~ horiz*myco, random = ~ 1 | block, data = colo3)
anova(colo.am.mod1)
plot(colo.am.mod1)
tuckey.colo.am <- emmeans(colo.am.mod1, pairwise ~ horiz*myco, adjust = "tukey")
mult.colo.am <- CLD(tuckey.colo.am,alpha=0.05,Letters=letters, adjust="tukey")
mult.colo.am.df <- rbind(as.data.frame(mult.colo.am[,1:4]), L)
mult.colo.am.df$myco <- factor(mult.colo.am.df$myco, levels = c('AM', 'Mixed', 'EcM'))

# By myco type
colo.am.myco.mod1 <- lme(pourc.AM.colo ~ myco, random = ~ 1 | block, data = colo3)
tuckey.colo.myco.am <- emmeans(colo.am.myco.mod1, pairwise ~ myco, adjust = "tukey")
mult.colo.myco.am <- CLD(tuckey.colo.myco.am,alpha=0.05,Letters=letters, adjust="tukey")

## Sequencing ####
## Data manipulation
ITS.guild$myco <- factor(ITS.guild$myco, levels = c('AM', 'Mixed', 'EcM'))
ITS.guild$horiz <- factor(ITS.guild$horiz, levels = c('B', 'Ae', 'H', 'F', 'L'))

ITS.sapro <- ITS.guild[ITS.guild$Trophic.Mode %in% "Saprotroph", ]
ITS.ecm <- ITS.guild[ITS.guild$Guild %in% "Ectomycorrhizal", ]
ITS.glomero <- ITS.guild[ITS.guild$Phylum %in% "p__Glomeromycota", ] #ITS.am <- ITS.guild[ITS.guild$Guild %in% "Arbuscular Mycorrhizal", ]
ITS.eric <- ITS.guild[ITS.guild$Guild %in% "Ericoid Mycorrhizal", ]
ITS.tot <- ITS.guild[ITS.guild$Kingdom %in% "k__Fungi", ]

## Total fungal reads per horizons and forest types
ITS.tot.sample <- ITS.tot  %>% 
  group_by(Sample) %>% 
  summarise(Abundance = sum(Abundance), ECM_perc = unique(ECM_perc), horiz = unique(horiz), myco = unique(myco), block = unique(block.y)) %>% 
  mutate(rel = Abundance/max(Abundance)*100)

tot.lme <- lme(rel ~ myco*horiz, random = ~ 1 | block,  weights = varPower(), data=ITS.tot.sample)
anova(tot.lme)
plot(tot.lme)
tot.tuckey <- emmeans(tot.lme, pairwise ~ horiz*myco, adjust = "tukey")
tot.mult <- CLD(tot.tuckey,alpha=0.05,Letters=letters, adjust="tukey")

## Sapro
ITS.sapro.sample <- ITS.sapro  %>% 
  group_by(Sample) %>% 
  summarise(Abundance = sum(Abundance), ECM_perc = unique(ECM_perc), horiz = unique(horiz), myco = unique(myco), block = unique(block.y)) %>% 
  mutate(rel = Abundance/max(Abundance)*100)

sapro.lme <- lme(rel ~ myco*horiz, random = ~ 1 | block,  weights = varIdent(form = ~ 1 | horiz), data=ITS.sapro.sample)
anova(sapro.lme)
plot(sapro.lme)
sapro.tuckey <- emmeans(sapro.lme, pairwise ~ horiz*myco, adjust = "tukey")
sapro.mult <- CLD(sapro.tuckey,alpha=0.05,Letters=letters, adjust="tukey")

## EcM
ITS.ecm.sample <- ITS.ecm  %>% 
  group_by(Sample) %>% 
  summarise(Abundance = sum(Abundance), ECM_perc = unique(ECM_perc), horiz = unique(horiz), myco = unique(myco), block = unique(block.y)) %>% 
  mutate(rel = Abundance/max(Abundance)*100)
ecm.lme <- lme(rel ~ myco*horiz, random = ~ 1 | block,  weights = varExp(), data=ITS.ecm.sample)
anova(ecm.lme)
plot(ecm.lme)

ecm.tuckey <- emmeans(ecm.lme, pairwise ~ myco*horiz, adjust = "tukey")
ecm.mult <- CLD(ecm.tuckey,alpha=0.05,Letters=letters, adjust="tukey")

## AM
ITS.am.sample <- ITS.glomero  %>% 
  group_by(Sample) %>% 
  summarise(Abundance = sum(Abundance), ECM_perc = unique(ECM_perc), horiz = unique(horiz), myco = unique(myco), block = unique(block.y)) %>% 
  mutate(rel = Abundance/max(Abundance)*100)
am.lme <- lme(rel ~ myco*horiz, random = ~ 1 | block,  weights = varExp(), data=ITS.am.sample)
anova(am.lme)
plot(am.lme)

am.tuckey <- emmeans(am.lme, pairwise ~ myco*horiz, adjust = "tukey")
am.mult <- CLD(am.tuckey,alpha=0.05,Letters=letters, adjust="tukey")

## Ericoid
ITS.eric.sample <- ITS.eric  %>% 
  group_by(Sample) %>% 
  summarise(Abundance = sum(Abundance), ECM_perc = unique(ECM_perc), horiz = unique(horiz), myco = unique(myco), block = unique(block.y)) %>% 
  mutate(rel = Abundance/max(Abundance)*100)

eric.lme <- lme(rel ~ myco*horiz, random = ~ 1 | block,  weights = varExp(), 
                data=ITS.eric.sample)
anova(eric.lme)
plot(eric.lme)
eric.tuckey <- emmeans(eric.lme, pairwise ~ myco*horiz, adjust = "tukey")
eric.mult <- CLD(eric.tuckey,alpha=0.05,Letters=letters, adjust="tukey")

## Plots ####
theme_depth <- theme_bw() + theme(panel.border = element_blank(),
                                  panel.grid.major = element_blank(),
                                  panel.grid.minor = element_blank(),
                                  axis.line = element_line(colour = "darkgrey"),
                                  text = element_text(size=8))
# Total
tot.rel.plot <-ggplot(tot.mult, aes(x = horiz, y = emmean, color = myco, shape = myco)) +
  geom_line(position=position_dodge(width = .2), aes(linetype=myco, color=myco, group = myco)) +
  geom_point(position=position_dodge(width = .2), size = 3) +
  geom_errorbar(aes(ymin = emmean-SE, ymax = emmean+SE), width = 0.2, position=position_dodge(width = .2)) +
  scale_color_manual(values = c('#1b9e77', '#d95f02', '#7570b3')) +
  coord_flip() +
  scale_y_continuous(position = "right") +
  labs(x = "Depth (cm)", y = "Total fungal\nrelative abundance (%)", linetype = "Forest", color = "Forest", shape = "Forest")+
  theme_depth


# EcM colonization
ecm.colo.rel.plot <- ggplot(mult.colo.ecm.df, aes(x = horiz, y = emmean, color = myco, shape = myco)) +
  geom_line(position=position_dodge(width = .2), aes(linetype=myco, color=myco, group = myco)) +
  geom_point(position=position_dodge(width = .2), size = 3) +
  geom_errorbar(aes(ymin = emmean-SE, ymax = emmean+SE), width = 0.2, position=position_dodge(width = .2)) +
  scale_color_manual(values = c('#1b9e77', '#d95f02', '#7570b3')) +
  coord_flip() +
  scale_y_continuous(position = "right") +
  labs(x = "", y = "EcM colonization\n of roots (%)", linetype = "Forest", color = "Forest", shape = "Forest")+
  theme_depth

# AM colonization
am.colo.rel.plot <-ggplot(mult.colo.am.df, aes(x = horiz, y = emmean, color = myco, shape = myco)) +
  geom_line(position=position_dodge(width = .2), aes(linetype=myco, color=myco, group = myco)) +
  geom_point(position=position_dodge(width = .2), size = 3) +
  geom_errorbar(aes(ymin = emmean-SE, ymax = emmean+SE), width = 0.2, position=position_dodge(width = .2)) +
  scale_color_manual(values = c('#1b9e77', '#d95f02', '#7570b3')) +
  coord_flip() +
  scale_y_continuous(position = "right") +
  labs(x = "Soil profile (horizons)", y = "AM colonization\nof roots (%)", linetype = "Forest", color = "Forest", shape = "Forest") +
  theme_depth

# Sapro
sapro.rel.plot <-ggplot(sapro.mult, aes(x = horiz, y = emmean, color = myco, shape = myco)) +
  geom_line(position=position_dodge(width = .2), aes(linetype=myco, color=myco, group = myco)) +
  geom_point(position=position_dodge(width = .2), size = 3) +
  geom_errorbar(aes(ymin = emmean-SE, ymax = emmean+SE), width = 0.2, position=position_dodge(width = .2)) +
  scale_color_manual(values = c('#1b9e77', '#d95f02', '#7570b3')) +
  coord_flip() +
  scale_y_continuous(position = "right") +
  labs(x = "", y = "Saptrophic fungi\nrelative abundance (%)", linetype = "Forest", color = "Forest", shape = "Forest")+
  theme_depth

#EcM
ecm.rel.plot <-ggplot(ecm.mult, aes(x = horiz, y = emmean, color = myco, shape = myco)) +
  geom_line(position=position_dodge(width = .2), aes(linetype=myco, color=myco, group = myco)) +
  geom_point(position=position_dodge(width = .2), size = 3) +
  geom_errorbar(aes(ymin = emmean-SE, ymax = emmean+SE), width = 0.2, position=position_dodge(width = .2)) +
  scale_color_manual(values = c('#1b9e77', '#d95f02', '#7570b3')) +
  coord_flip() +
  scale_y_continuous(position = "right") +
  labs(x = "", y = "EcM fungi\nrelative abundance (%)", linetype = "Forest", color = "Forest", shape = "Forest")+
  theme_depth

# AM
am.rel.plot <-ggplot(am.mult, aes(x = horiz, y = emmean, color = myco, shape = myco)) +
  geom_line(position=position_dodge(width = .2), aes(linetype=myco, color=myco, group = myco)) +
  geom_point(position=position_dodge(width = .2), size = 3) +
  geom_errorbar(aes(ymin = emmean-SE, ymax = emmean+SE), width = 0.2, position=position_dodge(width = .2)) +
  scale_color_manual(values = c('#1b9e77', '#d95f02', '#7570b3')) +
  coord_flip() +
  scale_y_continuous(position = "right") +
  labs(x = "Soil profile (horizons)", y = "AM fungi\nrelative abundance (%)", linetype = "Forest", color = "Forest", shape = "Forest")+
  theme_depth

eric.rel.plot <-ggplot(eric.mult, aes(x = horiz, y = emmean, color = myco, shape = myco)) +
  geom_line(position=position_dodge(width = .2), aes(linetype=myco, color=myco, group = myco)) +
  geom_point(position=position_dodge(width = .2), size = 3) +
  geom_errorbar(aes(ymin = emmean-SE, ymax = emmean+SE), width = 0.2, position=position_dodge(width = .2)) +
  scale_color_manual(values = c('#1b9e77', '#d95f02', '#7570b3')) +
  coord_flip() +
  scale_y_continuous(position = "right") +
  labs(x = "", y = "ErM fungi\nrelative abundance (%)", linetype = "Forest", color = "Forest", shape = "Forest")+
  theme_depth

guild.rel.plot <- ggarrange(am.colo.rel.plot, ecm.colo.rel.plot, sapro.rel.plot, 
                            am.rel.plot, ecm.rel.plot, eric.rel.plot,
                            nrow = 2, ncol = 3, common.legend = TRUE, legend="right",
                            labels = c('a', 'b', 'c', 'd', 'e', 'f'), font.label = list(size = 8))

