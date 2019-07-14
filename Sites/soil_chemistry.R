#### comparing environmental data among forest and horizons ####

## Packages ####
library(gridExtra)
library(grid)
library(ggplot2); theme_set(theme_bw())
library(emmeans)
library(nlme)
library(ggpubr)
library(dplyr)

## Get data ####
load("data/soil3.rda")
load("data/root2.rda")

# Merge df
soil <- left_join(soil3, root2)

## Reorder factors
soil$horiz = factor(soil$horiz, levels = c('L', 'F', 'H', 'Ae', 'B'))
soil$myco = factor(soil$myco, levels = c('AM', 'mixed', 'ECM'))

# Subset by  forest type
soil.AM <- soil[soil$myco %in% "AM", ]
soil.ECM <- soil[soil$myco %in% "ECM", ]
soil.mixed <- soil[soil$myco %in% "mixed", ]

## Mixed effects linear Model ####
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
#write.csv(lme.results, "Sites/soil.results.lme.csv")

# Just check Tukey HSD
pH.lme <- lme(pHCaCl2 ~ horiz * myco, random = ~ 1 | block,  weights = varIdent(form = ~ 1 | horiz), data=soil)
tuckey <- emmeans(pH.lme, pairwise ~ horiz * myco, adjust = "tukey")
mult <- CLD(tuckey,alpha=0.05,Letters=letters, adjust="tukey")

CN.lme <- lme(CN ~ horiz * myco, random = ~ 1 | block, weights = NULL,data=soil)
tuckey <- emmeans(CN.lme, pairwise ~ myco, adjust = "tukey")
mult <- CLD(tuckey,alpha=0.05,Letters=letters, adjust="tukey")

orgP.lme <- lme(orgP ~ horiz * myco, random = ~ 1 | block, weights = varExp(), data=soil)
tuckey <- emmeans(orgP.lme, pairwise ~ myco, adjust = "tukey")
mult <- CLD(tuckey,alpha=0.05,Letters=letters, adjust="tukey")

inorgP.lme <- lme(inorgP ~ horiz * myco, random = ~ 1 | block, weights = varIdent(form = ~ 1 | horiz), data=soil)
tuckey <- emmeans(inorgP.lme, pairwise ~ myco, adjust = "tukey")
mult <- CLD(tuckey,alpha=0.05,Letters=letters, adjust="tukey")

#### Total C ####
## AM
lme.C.AM <- lme(totalC ~ horiz,  random = ~ 1 | block, weights = varIdent(form = ~ 1 | horiz), data=soil.AM)
tuckey.C.AM <- emmeans(lme.C.AM, pairwise ~ horiz, adjust = "tukey")
mult.C.AM <- as.data.frame(CLD(tuckey.C.AM,alpha=0.05,Letters=letters, adjust="tukey"))
mult.C.AM <- mult.C.AM[order(mult.C.AM$horiz),]
mult.C.AM$chem <- as.factor(rep("totalC", length(unique(soil$horiz))))
mult.C.AM$myco <- as.factor(rep("AM", length(unique(soil$horiz))))

## ECM
lme.C.ECM <- lme(totalC ~ horiz,  random = ~ 1 | block, weights = varIdent(form = ~ 1 | horiz), data=soil.ECM)
tuckey.C.ECM <- emmeans(lme.C.ECM, pairwise ~ horiz, adjust = "tukey")
mult.C.ECM <- as.data.frame(CLD(tuckey.C.ECM,alpha=0.05,Letters=letters, adjust="tukey"))
mult.C.ECM <- mult.C.ECM[order(mult.C.ECM$horiz),]
mult.C.ECM$chem <- as.factor(rep("totalC", length(unique(soil$horiz))))
mult.C.ECM$myco <- as.factor(rep("ECM", length(unique(soil$horiz))))

## Mixed
lme.C.mixed <- lme(totalC ~ horiz,  random = ~ 1 | block, weights = varIdent(form = ~ 1 | horiz), data=soil.mixed)
tuckey.C.mixed <- emmeans(lme.C.mixed, pairwise ~ horiz, adjust = "tukey")
mult.C.mixed <- as.data.frame(CLD(tuckey.C.mixed,alpha=0.05,Letters=letters, adjust="tukey"))
mult.C.mixed <- mult.C.mixed[order(mult.C.mixed$horiz),]
mult.C.mixed$chem <- as.factor(rep("totalC", length(unique(soil$horiz))))
mult.C.mixed$myco <- as.factor(rep("mixed", length(unique(soil$horiz))))

#### Total N ####
## AM
lme.N.AM <- lme(totalN ~ horiz,  random = ~ 1 | block, weights = varPower(), data=soil.AM)
tuckey.N.AM <- emmeans(lme.N.AM, pairwise ~ horiz, adjust = "tukey")
mult.N.AM <- as.data.frame(CLD(tuckey.N.AM,alpha=0.05,Letters=letters, adjust="tukey"))
mult.N.AM <- mult.N.AM[order(mult.N.AM$horiz),]
mult.N.AM$chem <- as.factor(rep("totalN", length(unique(soil$horiz))))
mult.N.AM$myco <- as.factor(rep("AM", length(unique(soil$horiz))))

## ECM
lme.N.ECM <- lme(totalN ~ horiz,  random = ~ 1 | block, weights = NULL, data=soil.ECM)
tuckey.N.ECM <- emmeans(lme.N.ECM, pairwise ~ horiz, adjust = "tukey")
mult.N.ECM <- as.data.frame(CLD(tuckey.N.ECM,alpha=0.05,Letters=letters, adjust="tukey"))
mult.N.ECM <- mult.N.ECM[order(mult.N.ECM$horiz),]
mult.N.ECM$chem <- as.factor(rep("totalN", length(unique(soil$horiz))))
mult.N.ECM$myco <- as.factor(rep("ECM", length(unique(soil$horiz))))

## Mixed
lme.N.mixed <- lme(totalN ~ horiz,  random = ~ 1 | block, weights = varIdent(form = ~ 1 | horiz), data=soil.mixed)
tuckey.N.mixed <- emmeans(lme.N.mixed, pairwise ~ horiz, adjust = "tukey")
mult.N.mixed <- as.data.frame(CLD(tuckey.N.mixed,alpha=0.05,Letters=letters, adjust="tukey"))
mult.N.mixed <- mult.N.mixed[order(mult.N.mixed$horiz),]
mult.N.mixed$chem <- as.factor(rep("totalN", length(unique(soil$horiz))))
mult.N.mixed$myco <- as.factor(rep("mixed", length(unique(soil$horiz))))

#### CN ratio ####
## AM
lme.CN.AM <- lme(totalC/totalN ~ horiz,  random = ~ 1 | block, weights = NULL, data=soil.AM)
tuckey.CN.AM <- emmeans(lme.CN.AM, pairwise ~ horiz, adjust = "tukey")
mult.CN.AM <- as.data.frame(CLD(tuckey.CN.AM,alpha=0.05,Letters=letters, adjust="tukey"))
mult.CN.AM <- mult.CN.AM[order(mult.CN.AM$horiz),]
mult.CN.AM$chem <- as.factor(rep("CN", length(unique(soil$horiz))))
mult.CN.AM$myco <- as.factor(rep("AM", length(unique(soil$horiz))))

## ECM
lme.CN.ECM <- lme(totalC/totalN ~ horiz,  random = ~ 1 | block, weights = NULL, data=soil.ECM)
tuckey.CN.ECM <- emmeans(lme.CN.ECM, pairwise ~ horiz, adjust = "tukey")
mult.CN.ECM <- as.data.frame(CLD(tuckey.CN.ECM,alpha=0.05,Letters=letters, adjust="tukey"))
mult.CN.ECM <- mult.CN.ECM[order(mult.CN.ECM$horiz),]
mult.CN.ECM$chem <- as.factor(rep("CN", length(unique(soil$horiz))))
mult.CN.ECM$myco <- as.factor(rep("ECM", length(unique(soil$horiz))))

## mixed
lme.CN.mixed <- lme(totalC/totalN ~ horiz,  random = ~ 1 | block, weights = NULL, data=soil.mixed)
tuckey.CN.mixed <- emmeans(lme.CN.mixed, pairwise ~ horiz, adjust = "tukey")
mult.CN.mixed <- as.data.frame(CLD(tuckey.CN.mixed,alpha=0.05,Letters=letters, adjust="tukey"))
mult.CN.mixed <- mult.CN.mixed[order(mult.CN.mixed$horiz),]
mult.CN.mixed$chem <- as.factor(rep("CN", length(unique(soil$horiz))))
mult.CN.mixed$myco <- as.factor(rep("mixed", length(unique(soil$horiz))))

#### Organic P ####
## AM
lme.orgP.AM <- lme(orgP ~ horiz,  random = ~ 1 | block, weights = varIdent(form = ~ 1 | horiz), data=soil.AM)
tuckey.orgP.AM <- emmeans(lme.orgP.AM, pairwise ~ horiz, adjust = "tukey")
mult.orgP.AM <- as.data.frame(CLD(tuckey.orgP.AM,alpha=0.05,Letters=letters, adjust="tukey"))
mult.orgP.AM <- mult.orgP.AM[order(mult.orgP.AM$horiz),]
mult.orgP.AM$chem <- as.factor(rep("orgP", length(unique(soil$horiz))))
mult.orgP.AM$myco <- as.factor(rep("AM", length(unique(soil$horiz))))

## ECM
lme.orgP.ECM <- lme(orgP ~ horiz,  random = ~ 1 | block, weights = varIdent(form = ~ 1 | horiz), data=soil.ECM)
tuckey.orgP.ECM <- emmeans(lme.orgP.ECM, pairwise ~ horiz, adjust = "tukey")
mult.orgP.ECM <- as.data.frame(CLD(tuckey.orgP.ECM,alpha=0.05,Letters=letters, adjust="tukey"))
mult.orgP.ECM <- mult.orgP.ECM[order(mult.orgP.ECM$horiz),]
mult.orgP.ECM$chem <- as.factor(rep("orgP", length(unique(soil$horiz))))
mult.orgP.ECM$myco <- as.factor(rep("ECM", length(unique(soil$horiz))))

## Mixed
lme.orgP.mixed <- lme(orgP ~ horiz,  random = ~ 1 | block, weights = varIdent(form = ~ 1 | horiz), data=soil.mixed) # weights = NULL lower AIC and BIC but not equal variance
tuckey.orgP.mixed <- emmeans(lme.orgP.mixed, pairwise ~ horiz, adjust = "tukey")
mult.orgP.mixed <- as.data.frame(CLD(tuckey.orgP.mixed,alpha=0.05,Letters=letters, adjust="tukey"))
mult.orgP.mixed <- mult.orgP.mixed[order(mult.orgP.mixed$horiz),]
mult.orgP.mixed$chem <- as.factor(rep("orgP", length(unique(soil$horiz))))
mult.orgP.mixed$myco <- as.factor(rep("mixed", length(unique(soil$horiz))))

#### Inorganic P ####
## AM
lme.inorgP.AM <- lme(inorgP ~ horiz,  random = ~ 1 | block, weights = varIdent(form = ~ 1 | horiz), data=soil.AM)
tuckey.inorgP.AM <- emmeans(lme.inorgP.AM, pairwise ~ horiz, adjust = "tukey")
mult.inorgP.AM <- as.data.frame(CLD(tuckey.inorgP.AM,alpha=0.05,Letters=letters, adjust="tukey"))
mult.inorgP.AM <- mult.inorgP.AM[order(mult.inorgP.AM$horiz),]
mult.inorgP.AM$chem <- as.factor(rep("inorgP", length(unique(soil$horiz))))
mult.inorgP.AM$myco <- as.factor(rep("AM", length(unique(soil$horiz))))

## ECM
lme.inorgP.ECM <- lme(inorgP ~ horiz,  random = ~ 1 | block, weights = varIdent(form = ~ 1 | horiz), data=soil.ECM)
tuckey.inorgP.ECM <- emmeans(lme.inorgP.ECM, pairwise ~ horiz, adjust = "tukey")
mult.inorgP.ECM <- as.data.frame(CLD(tuckey.inorgP.ECM,alpha=0.05,Letters=letters, adjust="tukey"))
mult.inorgP.ECM <- mult.inorgP.ECM[order(mult.inorgP.ECM$horiz),]
mult.inorgP.ECM$chem <- as.factor(rep("inorgP", length(unique(soil$horiz))))
mult.inorgP.ECM$myco <- as.factor(rep("ECM", length(unique(soil$horiz))))

## Mixed
lme.inorgP.mixed <- lme(inorgP ~ horiz,  random = ~ 1 | block, weights = varIdent(form = ~ 1 | horiz), data=soil.mixed)
tuckey.inorgP.mixed <- emmeans(lme.inorgP.mixed, pairwise ~ horiz, adjust = "tukey")
mult.inorgP.mixed <- as.data.frame(CLD(tuckey.inorgP.mixed,alpha=0.05,Letters=letters, adjust="tukey"))
mult.inorgP.mixed <- mult.inorgP.mixed[order(mult.inorgP.mixed$horiz),]
mult.inorgP.mixed$chem <- as.factor(rep("inorgP", length(unique(soil$horiz))))
mult.inorgP.mixed$myco <- as.factor(rep("mixed", length(unique(soil$horiz))))

#### Bray P ####
## AM
lme.BrayP.AM <- lme(BrayP ~ horiz,  random = ~ 1 | block, weights = varIdent(form = ~ 1 | horiz), data=soil.AM)
tuckey.BrayP.AM <- emmeans(lme.BrayP.AM, pairwise ~ horiz, adjust = "tukey")
mult.BrayP.AM <- as.data.frame(CLD(tuckey.BrayP.AM,alpha=0.05,Letters=letters, adjust="tukey"))
mult.BrayP.AM <- mult.BrayP.AM[order(mult.BrayP.AM$horiz),]
mult.BrayP.AM$chem <- as.factor(rep("BrayP", length(unique(soil$horiz))))
mult.BrayP.AM$myco <- as.factor(rep("AM", length(unique(soil$horiz))))

## ECM
lme.BrayP.ECM <- lme(BrayP ~ horiz,  random = ~ 1 | block, weights = varPower(), data=soil.ECM)
tuckey.BrayP.ECM <- emmeans(lme.BrayP.ECM, pairwise ~ horiz, adjust = "tukey")
mult.BrayP.ECM <- as.data.frame(CLD(tuckey.BrayP.ECM,alpha=0.05,Letters=letters, adjust="tukey"))
mult.BrayP.ECM <- mult.BrayP.ECM[order(mult.BrayP.ECM$horiz),]
mult.BrayP.ECM$chem <- as.factor(rep("BrayP", length(unique(soil$horiz))))
mult.BrayP.ECM$myco <- as.factor(rep("ECM", length(unique(soil$horiz))))

## Mixed
lme.BrayP.mixed <- lme(BrayP ~ horiz,  random = ~ 1 | block, weights = varPower(), data=soil.mixed)
tuckey.BrayP.mixed <- emmeans(lme.BrayP.mixed, pairwise ~ horiz, adjust = "tukey")
mult.BrayP.mixed <- as.data.frame(CLD(tuckey.BrayP.mixed,alpha=0.05,Letters=letters, adjust="tukey"))
mult.BrayP.mixed <- mult.BrayP.mixed[order(mult.BrayP.mixed$horiz),]
mult.BrayP.mixed$chem <- as.factor(rep("BrayP", length(unique(soil$horiz))))
mult.BrayP.mixed$myco <- as.factor(rep("mixed", length(unique(soil$horiz))))

#### pH CaCl2 ####
#### AM
lme.pHCaCl2.AM <- lme(pHCaCl2 ~ horiz,  random = ~ 1 | block, weights = NULL, data=soil.AM)
tuckey.pHCaCl2.AM <- emmeans(lme.pHCaCl2.AM, pairwise ~ horiz, adjust = "tukey")
mult.pHCaCl2.AM <- as.data.frame(CLD(tuckey.pHCaCl2.AM,alpha=0.05,Letters=letters, adjust="tukey"))
mult.pHCaCl2.AM <- mult.pHCaCl2.AM[order(mult.pHCaCl2.AM$horiz),]
mult.pHCaCl2.AM$chem <- as.factor(rep("pHCaCl2", length(unique(soil$horiz))))
mult.pHCaCl2.AM$myco <- as.factor(rep("AM", length(unique(soil$horiz))))

## ECM
lme.pHCaCl2.ECM <- lme(pHCaCl2 ~ horiz,  random = ~ 1 | block, weights = NULL, data=soil.ECM)
tuckey.pHCaCl2.ECM <- emmeans(lme.pHCaCl2.ECM, pairwise ~ horiz, adjust = "tukey")
mult.pHCaCl2.ECM <- as.data.frame(CLD(tuckey.pHCaCl2.ECM,alpha=0.05,Letters=letters, adjust="tukey"))
mult.pHCaCl2.ECM <- mult.pHCaCl2.ECM[order(mult.pHCaCl2.ECM$horiz),]
mult.pHCaCl2.ECM$chem <- as.factor(rep("pHCaCl2", length(unique(soil$horiz))))
mult.pHCaCl2.ECM$myco <- as.factor(rep("ECM", length(unique(soil$horiz))))

## Mixed
lme.pHCaCl2.mixed <- lme(pHCaCl2 ~ horiz,  random = ~ 1 | block, weights = NULL, data=soil.mixed)
tuckey.pHCaCl2.mixed <- emmeans(lme.pHCaCl2.mixed, pairwise ~ horiz, adjust = "tukey")
mult.pHCaCl2.mixed <- as.data.frame(CLD(tuckey.pHCaCl2.mixed,alpha=0.05,Letters=letters, adjust="tukey"))
mult.pHCaCl2.mixed <- mult.pHCaCl2.mixed[order(mult.pHCaCl2.mixed$horiz),]
mult.pHCaCl2.mixed$chem <- as.factor(rep("pHCaCl2", length(unique(soil$horiz))))
mult.pHCaCl2.mixed$myco <- as.factor(rep("mixed", length(unique(soil$horiz))))

#### ECEC ####
## AM
lme.ECEC.AM <- lme(ECEC ~ horiz,  random = ~ 1 | block, weights = varIdent(form = ~ 1 | horiz), data=soil.AM)
tuckey.ECEC.AM <- emmeans(lme.ECEC.AM, pairwise ~ horiz, adjust = "tukey")
mult.ECEC.AM <- as.data.frame(CLD(tuckey.ECEC.AM,alpha=0.05,Letters=letters, adjust="tukey"))
mult.ECEC.AM <- mult.ECEC.AM[order(mult.ECEC.AM$horiz),]
mult.ECEC.AM$chem <- as.factor(rep("ECEC", length(unique(soil$horiz))))
mult.ECEC.AM$myco <- as.factor(rep("AM", length(unique(soil$horiz))))

## ECM
lme.ECEC.ECM <- lme(ECEC ~ horiz,  random = ~ 1 | block, weights = varIdent(form = ~ 1 | horiz), data=soil.ECM)
tuckey.ECEC.ECM <- emmeans(lme.ECEC.ECM, pairwise ~ horiz, adjust = "tukey")
mult.ECEC.ECM <- as.data.frame(CLD(tuckey.ECEC.ECM,alpha=0.05,Letters=letters, adjust="tukey"))
mult.ECEC.ECM <- mult.ECEC.ECM[order(mult.ECEC.ECM$horiz),]
mult.ECEC.ECM$chem <- as.factor(rep("ECEC", length(unique(soil$horiz))))
mult.ECEC.ECM$myco <- as.factor(rep("ECM", length(unique(soil$horiz))))

## Mixed
lme.ECEC.mixed <- lme(ECEC ~ horiz,  random = ~ 1 | block, weights = varIdent(form = ~ 1 | horiz), data=soil.mixed)
tuckey.ECEC.mixed <- emmeans(lme.ECEC.mixed, pairwise ~ horiz, adjust = "tukey")
mult.ECEC.mixed <- as.data.frame(CLD(tuckey.ECEC.mixed,alpha=0.05,Letters=letters, adjust="tukey"))
mult.ECEC.mixed <- mult.ECEC.mixed[order(mult.ECEC.mixed$horiz),]
mult.ECEC.mixed$chem <- as.factor(rep("ECEC", length(unique(soil$horiz))))
mult.ECEC.mixed$myco <- as.factor(rep("mixed", length(unique(soil$horiz))))


#### Merge all data ####
mult.chem <- rbind(mult.C.AM, mult.C.ECM, mult.C.mixed,
                   mult.N.AM, mult.N.ECM, mult.N.mixed,
                   mult.pHCaCl2.AM, mult.pHCaCl2.ECM, mult.pHCaCl2.mixed,
                   mult.orgP.AM, mult.orgP.ECM, mult.orgP.mixed,
                   mult.inorgP.AM, mult.inorgP.ECM, mult.inorgP.mixed,
                   mult.BrayP.AM, mult.BrayP.ECM, mult.BrayP.mixed,
                   mult.ECEC.AM, mult.ECEC.ECM, mult.ECEC.mixed,
                   mult.CN.AM, mult.CN.ECM, mult.CN.mixed)

## Reorder factors
mult.chem$horiz = factor(mult.chem$horiz, levels = c('B', 'Ae', 'H', 'F', 'L'))
mult.chem$myco = factor(mult.chem$myco, levels = c('AM', 'mixed', 'ECM'))

#### Plot one by one and arrange ####
## Set plot theme
mytheme <- theme_bw() + theme(text = element_text(size=13), 
                              strip.text.y = element_blank(), 
                              strip.background = element_rect(fill = NA, colour = "black"),
                              panel.spacing = unit(1, "lines"), 
                              panel.grid.major.y = element_blank())
theme_set(mytheme)
myco_names <- c(`AM` = "AM forest", `mixed` = "Mixed forest", `ECM` = "EcM forest")

pH.plot <- ggplot(mult.chem[mult.chem$chem == "pHCaCl2",], aes(x = horiz, y = emmean, fill = horiz)) + 
  geom_bar(stat = "identity", show.legend = FALSE, colour="black", width =.8) +
  geom_errorbar(aes(ymin = emmean, ymax = emmean+SE), width = 0.2) +
  scale_fill_manual(values = c('darkorange2', 'grey60', 'grey1', 'sienna4', 'darkgreen')) +
  #geom_text(aes(y = emmean, label = .group, vjust = 0.4, hjust	= 1.1)) + 
  facet_grid(chem~myco, scales = "fixed",  labeller = as_labeller(myco_names)) +  
  labs(x="", y= expression(paste("pH (in CaCl"[2],")"))) + 
  coord_flip()
  
totalC.plot <- ggplot(mult.chem[mult.chem$chem == "totalC",], aes(x = horiz, y = emmean, fill = horiz)) + 
  geom_bar(stat = "identity", show.legend = FALSE, colour="black", width =.8) +
  geom_errorbar(aes(ymin = emmean, ymax = emmean+SE), width = 0.2) +
  scale_fill_manual(values = c('darkorange2', 'grey60', 'grey1', 'sienna4', 'darkgreen')) +
  #geom_text(aes(y = emmean, label = .group, vjust = 0.4, hjust	= 1.1)) + 
  facet_grid(chem~myco, scales = "fixed", labeller = as_labeller(myco_names)) +  
  labs(x="", y= expression(paste("Total C (kg m"^"-3",")"))) + 
  coord_flip()

totalN.plot <- ggplot(mult.chem[mult.chem$chem == "totalN",], aes(x = horiz, y = emmean, fill = horiz)) + 
  geom_bar(stat = "identity", show.legend = FALSE, colour="black", width =.8) +
  geom_errorbar(aes(ymin = emmean, ymax = emmean+SE), width = 0.2) +
  scale_fill_manual(values = c('darkorange2', 'grey60', 'grey1', 'sienna4', 'darkgreen')) +
  #geom_text(aes(y = emmean, label = .group, vjust = 0.4, hjust	= 1.1)) + 
  facet_grid(chem~myco, scales = "fixed", labeller = as_labeller(myco_names)) +  
  labs(x="", y= expression(paste("Total N (kg m"^"-3",")"))) + 
  coord_flip()

CN.plot <- ggplot(mult.chem[mult.chem$chem == "CN",], aes(x = horiz, y = emmean, fill = horiz)) + 
  geom_bar(stat = "identity", show.legend = FALSE, colour="black", width =.8) +
  geom_errorbar(aes(ymin = emmean, ymax = emmean+SE), width = 0.2) +
  scale_fill_manual(values = c('darkorange2', 'grey60', 'grey1', 'sienna4', 'darkgreen')) +
  #geom_text(aes(y = emmean, label = .group, vjust = 0.4, hjust	= 1.1)) + 
  facet_grid(chem~myco, scales = "fixed", labeller = as_labeller(myco_names)) +  
  labs(x="", y= expression(paste("C:N ratio"))) + 
  coord_flip()

orgP.plot <- ggplot(mult.chem[mult.chem$chem == "orgP",], aes(x = horiz, y = emmean, fill = horiz)) + 
  geom_bar(stat = "identity", show.legend = FALSE, colour="black", width =.8) +
  geom_errorbar(aes(ymin = emmean, ymax = emmean+SE), width = 0.2) +
  scale_fill_manual(values = c('darkorange2', 'grey60', 'grey1', 'sienna4', 'darkgreen')) +
  #geom_text(aes(y = emmean, label = .group, vjust = 0.4, hjust	= 1.1)) + 
  facet_grid(chem~myco, scales = "fixed", labeller = as_labeller(myco_names)) +  
  labs(x="", y= expression(paste("Organic P (kg m"^"-3",")"))) + 
  coord_flip()

inorgP.plot <- ggplot(mult.chem[mult.chem$chem == "inorgP",], aes(x = horiz, y = emmean, fill = horiz)) + 
  geom_bar(stat = "identity", show.legend = FALSE, colour="black", width =.8) +
  geom_errorbar(aes(ymin = emmean, ymax = emmean+SE), width = 0.2) +
  scale_fill_manual(values = c('darkorange2', 'grey60', 'grey1', 'sienna4', 'darkgreen')) +
  #geom_text(aes(y = emmean, label = .group, vjust = 0.4, hjust	= 1.1)) + 
  facet_grid(chem~myco, scales = "fixed", labeller = as_labeller(myco_names)) +  
  labs(x="", y= expression(paste("Inorganic P (kg m"^"-3",")"))) + 
  coord_flip()

BrayP.plot <- ggplot(mult.chem[mult.chem$chem == "BrayP",], aes(x = horiz, y = emmean, fill = horiz)) + 
  geom_bar(stat = "identity", colour="black", width =.8, show.legend = FALSE) +
  geom_errorbar(aes(ymin = emmean, ymax = emmean+SE), width = 0.2) +
  scale_fill_manual(values = c('darkorange2', 'grey60', 'grey1', 'sienna4', 'darkgreen')) +
  #geom_text(aes(y = emmean, label = .group, vjust = 0.4, hjust	= 1.1)) + 
  facet_grid(chem~myco, scales = "fixed", labeller = as_labeller(myco_names)) +  
  labs(x="", y= expression(paste("Labile P (kg m"^"-3",")"))) + 
  coord_flip()

ECEC.plot <- ggplot(mult.chem[mult.chem$chem == "ECEC",], aes(x = horiz, y = emmean, fill = horiz)) + 
  geom_errorbar(aes(ymin = emmean, ymax = emmean+SE), width = 0.2) +
  geom_bar(stat = "identity", colour="black", width =.8, show.legend = FALSE) +
  scale_fill_manual(values = c('darkorange2', 'grey60', 'grey1', 'sienna4', 'darkgreen')) +
  #geom_text(aes(y = emmean, label = .group, vjust = 0.4, hjust	= 1.1)) + 
  facet_grid(chem~myco, scales = "fixed", labeller = as_labeller(myco_names)) +  
  labs(x="", y= expression(paste("ECEC (kg m"^"-3",")"))) + 
  coord_flip()

all.plot <- ggarrange(pH.plot, ECEC.plot, CN.plot, orgP.plot, inorgP.plot, BrayP.plot,
                      labels = c("a", "b", "c", "d", "e", "f"),
                      ncol = 1, nrow = 6)
all.plot
#ggsave("Sites/soil.chemistry.vol.SE.arrange.pdf", all.plot, width = 8, height = 10)
#ggsave("Sites/soil.chemistry.vol.SE.arrange.png", all.plot, width = 8, height = 10)

## Thickness ####
soil$horiz <- factor(soil$horiz, levels = c('B', 'Ae', 'H', 'F', 'L'))
ggplot(soil, aes(x = myco, y = thick.avg, fill = horiz)) +
  geom_bar(stat = "summary", fun.y = "mean")  +
  scale_fill_manual('Horizons',  labels= levels(soil$horiz), values = rev(c('darkgreen', 'sienna4', 'grey1', 'grey60', 'darkorange2'))) +
  scale_y_reverse() +
  guides(fill = guide_legend(reverse = TRUE))+
  labs(x = 'Forest', y = 'Thickness (cm)')
