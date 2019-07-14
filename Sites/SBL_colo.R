# Mycorrhizal colonisation analyses / SBL mycorrhizal plots

# load packages
library(dplyr)
library(ggplot2); theme_set(theme_bw())
library(ggpubr)
library(nlme)

# Data manipulation
load("data/colo3.rda")

colo3$AM.ECM[is.na(colo3$AM.ECM)] <- max(colo3$AM.ECM, na.rm = T) # NA if 0% ECM colonisation
colo3$AM.ECM <- colo3$AM.ECM/max(colo3$AM.ECM)

# remove all Ah horizons
colo3 <- droplevels(subset(colo3, !horiz %in% 'Ah'))

levels(colo3$myco)[levels(colo3$myco) == "AM"] <- "Maple"
levels(colo3$myco)[levels(colo3$myco) == "ECM"] <- "Beech"
levels(colo3$myco)[levels(colo3$myco) == "mixed"] <- "Mixed"

colo3$horiz <- factor(colo3$horiz, levels = c('F', 'H', 'Ae', 'B'))
colo3$myco <- factor(colo3$myco, levels = c('Maple', 'Mixed', 'Beech'))

#### Exploratory graphs ####
colo.plot <- ggplot(colo3, aes(x = myco, y = AM.ECM)) +
  geom_point(colour = "lightgrey") +
  stat_summary(fun.data = mean_se) +
  facet_wrap(~horiz, nrow = 1)+
  xlab('Forest') +
  ylab('AM:ECM colonisation ratio') +
  theme(axis.text.x = element_text(angle=45, vjust= 1, size = 10, hjust=1))
colo.plot

colo.AM.plot <- ggplot(colo3, aes(x = myco, y = pourc.AM.colo)) +
  geom_boxplot() +
  facet_wrap(~horiz, nrow = 1)+
  ylim(0,100) +
  theme(axis.text.x = element_text(angle=45, vjust= 1, hjust=1), panel.grid.major.x = element_blank()) +
  labs(x='', y='AM colonisation (%)')

colo.ECM.plot <- ggplot(colo3, aes(x = myco, y = pourc.ECM.colo)) +
  geom_boxplot() +
  facet_wrap(~horiz, nrow = 1)+
  ylim(0,100) +
  theme(axis.text.x = element_text(angle=45, vjust= 1, hjust=1), panel.grid.major.x = element_blank()) +
  labs(x='Forest', y='ECM colonisation (%)')

colo.plot <- ggarrange(colo.AM.plot, colo.ECM.plot, labels = c("a", "b"), nrow = 2, ncol = 1)
colo.plot
#ggsave("Sites/colo.plot.pdf", colo.plot, width = 7, height = 5)
#ggsave("Sites/colo.plot.png", colo.plot, width = 7, height = 5)
