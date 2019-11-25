## Mycorrhizal colonisation analyses / SBL mycorrhizal plots

## load packages and data ####
library(dplyr)
library(ggplot2); theme_set(theme_bw())
library(ggpubr)
library(nlme)

load("1.data/colo3.rda")

colo3$AM.ECM[is.na(colo3$AM.ECM)] <- max(colo3$AM.ECM, na.rm = T) # NA if 0% ECM colonisation
colo3$AM.ECM <- colo3$AM.ECM/max(colo3$AM.ECM)

# remove all Ah horizons
colo3 <- droplevels(subset(colo3, !horiz %in% 'Ah'))

levels(colo3$myco)[levels(colo3$myco) == "ECM"] <- "EcM"
levels(colo3$myco)[levels(colo3$myco) == "mixed"] <- "Mixed"

colo3$horiz <- factor(colo3$horiz, levels = c('F', 'H', 'Ae', 'B'))
colo3$myco <- factor(colo3$myco, levels = c('AM', 'Mixed', 'EcM'))

## ggplot theme for ploting
theme_colo <- theme_bw() + theme(panel.grid.major.x = element_blank() , 
                                 panel.grid.major.y = element_line( size=.1, color="grey"), 
                                 panel.grid.minor.y = element_line( size=.1, color="grey"), 
                                 strip.background = element_rect(fill="white"), 
                                 panel.spacing = unit(0, "mm"))

#### Models ####
## ECM
hist(colo3$pourc.ECM.colo)
colo.ecm.mod1 <- lme(pourc.ECM.colo ~ horiz*myco, random = ~ 1 | block, data = colo3)
colo.ecm.mod2 <- update(colo.ecm.mod1, weights = varPower())
colo.ecm.mod3 <- update(colo.ecm.mod1, weights = varExp())
colo.ecm.mod4 <- update(colo.ecm.mod1, weights = varIdent(form = ~ 1 | horiz))
anova(colo.ecm.mod1, colo.ecm.mod2, colo.ecm.mod3, colo.ecm.mod4) # mod2 mod3  better
anova(colo.ecm.mod3); plot(colo.ecm.mod3)

# Interactions
tuckey.colo.ecm <- emmeans(colo.ecm.mod2, pairwise ~ horiz*myco, adjust = "tukey")
mult.colo.ecm <- CLD(tuckey.colo.ecm,alpha=0.05,Letters=letters, adjust="tukey")

colo.mod.ecm <- ggplot(mult.colo.ecm, aes(x = myco, y = emmean, color = myco)) +
  geom_point(stat = "identity") +
  geom_errorbar(aes(ymin = emmean-SE, ymax = emmean+SE), width = 0.2) +
  facet_wrap(~horiz, ncol = 4)+
  scale_color_manual(values = c('#1b9e77', '#d95f02', '#7570b3'), guide = FALSE) +
  labs(x = 'Forest', y = "Root colonization by EcM fungi (%)") +
  theme_colo
colo.mod.ecm

myco_names <- c(`AM` = "AM forest ", `Mixed` = "Mixed forest", `EcM` = "EcM forest")
colo.mod.horiz.ecm <- ggplot(mult.colo.ecm, aes(x = horiz, y = emmean)) +
  geom_point(stat = "identity") +
  geom_errorbar(aes(ymin = emmean-SE, ymax = emmean+SE), width = 0.2) +
  facet_grid(~myco, labeller = as_labeller(myco_names)) +
  labs(x = 'Horizons', y = "Root colonization by EcM fungi (%)") +
  theme_colo
colo.mod.horiz.ecm

## AM
hist(colo3$pourc.AM.colo)
colo.am.mod1 <- lme(pourc.AM.colo ~ horiz*myco, random = ~ 1 | block, data = colo3)
colo.am.mod2 <- update(colo.am.mod1, weights = varPower())
colo.am.mod3 <- update(colo.am.mod1, weights = varExp())
colo.am.mod4 <- update(colo.am.mod1, weights = varIdent(form = ~ 1 | horiz))
anova(colo.am.mod1, colo.am.mod2, colo.am.mod3, colo.am.mod4)
anova(colo.am.mod1); plot(colo.am.mod1)

## Interactions
tuckey.colo.am <- emmeans(colo.am.mod1, pairwise ~ horiz*myco, adjust = "tukey")
mult.colo.am <- CLD(tuckey.colo.am,alpha=0.05,Letters=letters, adjust="tukey")

colo.mod.am <- ggplot(mult.colo.am, aes(x = myco, y = emmean, color = myco)) +
  geom_point(stat = "identity") +
  geom_errorbar(aes(ymin = emmean-SE, ymax = emmean+SE), width = 0.2) +
  facet_wrap(~horiz, ncol = 4)+
  scale_color_manual(name = "Forest", labels = c(`AM` = "AM", `mixed` = "Mixed", `ECM` = "EcM"), values = c('#1b9e77', '#d95f02', '#7570b3'), guide = FALSE) +
  labs(x = '', y = "Root colonization by AM fungi (%)") +
  theme_colo
colo.mod.am

colo.mod.horiz.am <- ggplot(mult.colo.am, aes(x = horiz, y = emmean)) +
  geom_point(stat = "identity") +
  geom_errorbar(aes(ymin = emmean-SE, ymax = emmean+SE), width = 0.2) +
  facet_grid(~myco, labeller = as_labeller(myco_names)) +
  labs(x = '', y = "Root colonization by AM fungi (%)") +
  theme_colo
colo.mod.horiz.am

colo.mod.plot <- ggarrange(colo.mod.horiz.am, colo.mod.horiz.ecm, labels = c("a", "b"), nrow = 2, ncol = 1)

