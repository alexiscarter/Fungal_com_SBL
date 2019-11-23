## Comparing environmental data among forest and horizons ####

## Load packages and data ####
library(gridExtra)
library(grid)
library(ggplot2); theme_set(theme_bw())
library(emmeans)
library(nlme)
library(ggpubr)
library(dplyr)

load("data/soil3.rda")
load("data/root2.rda")

## Merge dataframes
soil <- left_join(soil3, root2)

## Reorder factors
soil$horiz = factor(soil$horiz, levels = c('L', 'F', 'H', 'Ae', 'B'))
soil$myco = factor(soil$myco, levels = c('AM', 'mixed', 'ECM'))

## Mixed effects linear models ####
pH.lm <- anova(lme(pHCaCl2 ~ horiz * myco, random = ~ 1 | block,  weights = varIdent(form = ~ 1 | horiz), data=soil))
totalC.lm <- anova(lme(totalC ~ horiz * myco, random = ~ 1 | block, weights = varIdent(form = ~ 1 | horiz),data=soil))
totalN.lm <- anova(lme(totalN ~ horiz * myco, random = ~ 1 | block, weights = varExp(), data=soil))
CN.lm <- anova(lme(CN ~ horiz * myco, random = ~ 1 | block, weights = NULL,data=soil))
totalP.lm <- anova(lme(totalP ~ horiz * myco, random = ~ 1 | block, weights = varExp(), data=soil))
orgP.lm <- anova(lme(orgP ~ horiz * myco, random = ~ 1 | block, weights = varExp(), data=soil))
inorgP.lm <- anova(lme(inorgP ~ horiz * myco, random = ~ 1 | block, weights = varIdent(form = ~ 1 | horiz), data=soil))
BrayP.lm <- anova(lme(BrayP ~ horiz * myco, random = ~ 1 | block, weights = varPower(), data=soil))
ECEC.lm <- anova(lme(ECEC ~ horiz * myco, random = ~ 1 | block, weights = NULL, data=soil))
BS.lm <- anova(lme(BS ~ horiz * myco, random = ~ 1 | block, weights = varIdent(form = ~ 1 | horiz), data=soil))

## Checked for variance 
# test <-
# test1 <- update(test, weights = varPower())
# test2 <- update(test, weights = varExp())
# test3 <- update(test, weights = varIdent(form = ~ 1 | horiz))
# anova(test, test1, test2, test3)
# plot(test)
# plot(test1)
# plot(test2)
# plot(test3)
# anova()

# Put results in a table
lme.results <- data.frame(X = c("Intercept", "Horiz", "Myco", "Interaction"),
                          pH.lm = pH.lm$`p-value`,
                          totalC.lm = totalC.lm$`p-value`,
                          totalN.lm = totalN.lm$`p-value`,
                          CN.lm = CN.lm$`p-value`,
                          totalP.lm = totalP.lm$`p-value`,
                          inorgP.lm = inorgP.lm$`p-value`,
                          orgP.lm = orgP.lm$`p-value`,
                          BrayP.lm = BrayP.lm$`p-value`,
                          ECEC.lm = ECEC.lm$`p-value`,
                          BS.lm = BS.lm$`p-value`)
lme.results[,c(2:9)] <- round(lme.results[,c(2:9)], 5)

## Mixed models one by one ####
## pH
pH.lme <- lme(pHCaCl2 ~ horiz * myco, random = ~ 1 | block,  weights = varIdent(form = ~ 1 | horiz), data=soil)
anova(pH.lme)
plot(pH.lme)
pH.tuckey <- emmeans(pH.lme, pairwise ~ horiz * myco, adjust = "tukey")
pH.mult <- CLD(pH.tuckey,alpha=0.05,Letters=letters, adjust="tukey")

ECEC.lme <- lme(ECEC ~ horiz * myco, random = ~ 1 | block, weights = NULL,data=soil)
anova(ECEC.lme)
plot(ECEC.lme)
ECEC.tuckey <- emmeans(ECEC.lme, pairwise ~ myco+horiz, adjust = "tukey")
ECEC.mult <- CLD(ECEC.tuckey,alpha=0.05,Letters=letters, adjust="tukey")

BS.lme <- lme(BS ~ horiz * myco, random = ~ 1 | block, weights = varIdent(form = ~ 1 | horiz),data=soil)
anova(BS.lme)
plot(BS.lme)
BS.tuckey <- emmeans(BS.lme, pairwise ~ myco+horiz, adjust = "tukey")
BS.mult <- CLD(BS.tuckey,alpha=0.05,Letters=letters, adjust="tukey")

CN.lme <- lme(CN ~ horiz * myco, random = ~ 1 | block, weights = NULL,data=soil)
anova(CN.lme)
plot(CN.lme)
CN.tuckey <- emmeans(CN.lme, pairwise ~ myco*horiz, adjust = "tukey")
CN.mult <- CLD(CN.tuckey,alpha=0.05,Letters=letters, adjust="tukey")

totalC.lme <- lme(totalC ~ horiz * myco, random = ~ 1 | block, weights = NULL,data=soil)
anova(totalC.lme)
plot(totalC.lme)
totalC.tuckey <- emmeans(totalC.lme, pairwise ~ myco+horiz, adjust = "tukey")
totalC.mult <- CLD(totalC.tuckey,alpha=0.05,Letters=letters, adjust="tukey")

totalN.lme <- lme(totalN ~ horiz * myco, random = ~ 1 | block, weights = varExp(),data=soil)
anova(totalN.lme)
plot(totalN.lme)
totalN.tuckey <- emmeans(totalN.lme, pairwise ~ myco+horiz, adjust = "tukey")
totalN.mult <- CLD(totalN.tuckey,alpha=0.05,Letters=letters, adjust="tukey")

orgP.lme <- lme(orgP ~ horiz * myco, random = ~ 1 | block, weights = varIdent(form = ~ 1 | horiz), data=soil)
AIC(orgP.lme)
anova(orgP.lme)
plot(orgP.lme)
orgP.tuckey <- emmeans(orgP.lme, pairwise ~ horiz * myco, adjust = "tukey")
orgP.mult <- CLD(orgP.tuckey,alpha=0.05,Letters=letters, adjust="tukey")

totalP.lme <- lme(totalP ~ horiz * myco, random = ~ 1 | block, weights = varIdent(form = ~ 1 | horiz), data=soil)
anova(totalP.lme)
plot(totalP.lme)
totalP.tuckey <- emmeans(totalP.lme, pairwise ~ horiz * myco, adjust = "tukey")
totalP.mult <- CLD(totalP.tuckey,alpha=0.05,Letters=letters, adjust="tukey")

BrayP.lme <- lme(BrayP ~ horiz * myco, random = ~ 1 | block, weights = varPower(), data=soil)
anova(BrayP.lme)
plot(BrayP.lme)
BrayP.tuckey <- emmeans(BrayP.lme, pairwise ~ horiz * myco, adjust = "tukey")
BrayP.mult <- CLD(BrayP.tuckey,alpha=0.05,Letters=letters, adjust="tukey")

## Plots ####
theme_depth <- theme_bw() + theme(panel.border = element_blank(),
                                  panel.grid.major = element_blank(),
                                  panel.grid.minor = element_blank(),
                                  axis.line = element_line(colour = "darkgrey"),
                                  text = element_text(size=8))

pH.plot <-ggplot(pH.mult, aes(x = horiz, y = emmean, color = myco, shape = myco)) +
  geom_line(position=position_dodge(width = .2), aes(linetype=myco, color=myco, group = myco)) +
  geom_point(position=position_dodge(width = .2), size = 2) +
  geom_errorbar(aes(ymin = emmean-SE, ymax = emmean+SE), width = 0.2, position=position_dodge(width = .2)) +
  scale_color_manual(values = c('#1b9e77', '#d95f02', '#7570b3')) +
  coord_flip() +
  scale_y_continuous(position = "right") +
  labs(x="Soil profile (horizons)", y= expression(paste("pH (in CaCl"[2],")")), linetype = "Forest", color = "Forest", shape = "Forest") +
  theme_depth

ECEC.plot <-ggplot(ECEC.mult, aes(x = horiz, y = emmean, color = myco, shape = myco)) +
  geom_line(position=position_dodge(width = .2), aes(linetype=myco, color=myco, group = myco)) +
  geom_point(position=position_dodge(width = .2), size = 2) +
  geom_errorbar(aes(ymin = emmean-SE, ymax = emmean+SE), width = 0.2, position=position_dodge(width = .2)) +
  scale_color_manual(values = c('#1b9e77', '#d95f02', '#7570b3')) +
  coord_flip() +
  scale_y_continuous(position = "right") +
  labs(x="", y= expression(paste("ECEC (kg m"^"-3",")")), linetype = "Forest", color = "Forest", shape = "Forest") +
  theme_depth

BS.plot <-ggplot(BS.mult, aes(x = horiz, y = emmean, color = myco, shape = myco)) +
  geom_line(position=position_dodge(width = .2), aes(linetype=myco, color=myco, group = myco)) +
  geom_point(position=position_dodge(width = .2), size = 2) +
  geom_errorbar(aes(ymin = emmean-SE, ymax = emmean+SE), width = 0.2, position=position_dodge(width = .2)) +
  scale_color_manual(values = c('#1b9e77', '#d95f02', '#7570b3')) +
  coord_flip() +
  scale_y_continuous(position = "right") +
  labs(x="", y= "Base Sat. (%)", linetype = "Forest", color = "Forest", shape = "Forest") +
  theme_depth

CN.plot <-ggplot(CN.mult, aes(x = horiz, y = emmean, color = myco, shape = myco)) +
  geom_line(position=position_dodge(width = .2), aes(linetype=myco, color=myco, group = myco)) +
  geom_point(position=position_dodge(width = .2), size = 2) +
  geom_errorbar(aes(ymin = emmean-SE, ymax = emmean+SE), width = 0.2, position=position_dodge(width = .2)) +
  scale_color_manual(values = c('#1b9e77', '#d95f02', '#7570b3')) +
  coord_flip() +
  scale_y_continuous(position = "right") +
  labs(x="", y= "C:N ratio", linetype = "Forest", color = "Forest", shape = "Forest") +
  theme_depth

totalC.plot <-ggplot(totalC.mult, aes(x = horiz, y = emmean, color = myco, shape = myco)) +
  geom_line(position=position_dodge(width = .2), aes(linetype=myco, color=myco, group = myco)) +
  geom_point(position=position_dodge(width = .2), size = 2) +
  geom_errorbar(aes(ymin = emmean-SE, ymax = emmean+SE), width = 0.2, position=position_dodge(width = .2)) +
  scale_color_manual(values = c('#1b9e77', '#d95f02', '#7570b3')) +
  coord_flip() +
  scale_y_continuous(position = "right") +
  labs(x="Soil profile (horizons)", y= expression(paste("Organic C (kg m"^"-3",")")), linetype = "Forest", color = "Forest", shape = "Forest") +
  theme_depth

totalN.plot <-ggplot(totalN.mult, aes(x = horiz, y = emmean, color = myco, shape = myco)) +
  geom_line(position=position_dodge(width = .2), aes(linetype=myco, color=myco, group = myco)) +
  geom_point(position=position_dodge(width = .2), size = 2) +
  geom_errorbar(aes(ymin = emmean-SE, ymax = emmean+SE), width = 0.2, position=position_dodge(width = .2)) +
  scale_color_manual(values = c('#1b9e77', '#d95f02', '#7570b3')) +
  coord_flip() +
  scale_y_continuous(position = "right") +
  labs(x="", y= expression(paste("Total N (kg m"^"-3",")")), linetype = "Forest", color = "Forest", shape = "Forest") +
  theme_depth

orgP.plot <-ggplot(orgP.mult, aes(x = horiz, y = emmean, color = myco, shape = myco)) +
  geom_line(position=position_dodge(width = .2), aes(linetype=myco, color=myco, group = myco)) +
  geom_point(position=position_dodge(width = .2), size = 2) +
  geom_errorbar(aes(ymin = emmean-SE, ymax = emmean+SE), width = 0.2, position=position_dodge(width = .2)) +
  scale_color_manual(values = c('#1b9e77', '#d95f02', '#7570b3')) +
  coord_flip() +
  scale_y_continuous(position = "right") +
  labs(x="", y= expression(paste("Organic P (kg m"^"-3",")")), linetype = "Forest", color = "Forest", shape = "Forest") +
  theme_depth

totalP.plot <-ggplot(totalP.mult, aes(x = horiz, y = emmean, color = myco, shape = myco)) +
  geom_line(position=position_dodge(width = .2), aes(linetype=myco, color=myco, group = myco)) +
  geom_point(position=position_dodge(width = .2), size = 2) +
  geom_errorbar(aes(ymin = emmean-SE, ymax = emmean+SE), width = 0.2, position=position_dodge(width = .2)) +
  scale_color_manual(values = c('#1b9e77', '#d95f02', '#7570b3')) +
  coord_flip() +
  scale_y_continuous(position = "right") +
  labs(x="Soil profile (horizons)", y= expression(paste("Total P (kg m"^"-3",")")), linetype = "Forest", color = "Forest", shape = "Forest") +
  theme_depth

BrayP.plot <-ggplot(BrayP.mult, aes(x = horiz, y = emmean, color = myco, shape = myco)) +
  geom_line(position=position_dodge(width = .2), aes(linetype=myco, color=myco, group = myco)) +
  geom_point(position=position_dodge(width = .2), size = 2) +
  geom_errorbar(aes(ymin = emmean-SE, ymax = emmean+SE), width = 0.2, position=position_dodge(width = .2)) +
  scale_color_manual(values = c('#1b9e77', '#d95f02', '#7570b3')) +
  coord_flip() +
  scale_y_continuous(position = "right") +
  labs(x="", y= expression(paste("Labile (Bray) P (kg m"^"-3",")")), linetype = "Forest", color = "Forest", shape = "Forest") +
  theme_depth

all.plot <- ggarrange(pH.plot, ECEC.plot, BS.plot, 
                      totalC.plot, totalN.plot, CN.plot, 
                      totalP.plot, orgP.plot, BrayP.plot,
                      labels = c("a", "b", "c", "d", "e", "f", "g", "h", "i"),
                      ncol = 3, nrow = 3, common.legend = TRUE, legend="right", font.label = list(size = 8))
all.plot

