# ITS Paired
library(phyloseq)
library(dplyr)
library(ggplot2)
library(vegan)
library(ddpcr)
library(RVAideMemoire)
library(stringi)
library(gridExtra)
library(labdsv)
library(emmeans)
library(ggedit)


# Package version
packageVersion("phyloseq")
packageVersion("dplyr")
packageVersion("ggplot2"); theme_set(theme_bw())
packageVersion("vegan")
packageVersion("ddpcr") # use for the 'quiet()' function
packageVersion("RVAideMemoire")
packageVersion("stringi") 
packageVersion("gridExtra") 
packageVersion("labdsv") 

#### Filtering and transformation ####
load("data/ps.paired.ITS.29.05.minboot50.rda")
ntaxa(ps)

# Remove ESV with no "real" reads
ps = filter_taxa(ps, function(x) sum(x) > 0, TRUE)

# Subset of the dataset with only fungi
ps.fungi = subset_taxa(ps, Kingdom=="k__Fungi")

# Remove singletons and doubletons
ps.fungi.nosd = filter_taxa(ps.fungi, function(x) sum(x) > 2, TRUE)
sum(colSums(otu_table(ps.fungi.nosd))) # 2711899 sequences
ntaxa(ps.fungi.nosd) # 2521 ASV

## Shifted log transformation
ps.nosd.log <- transform_sample_counts(ps.fungi.nosd, function(x) log(x+1))

## Agglomerate taxa
## At the species level
ps.nosd.log.sp = tax_glom(ps.nosd.log, "Species", NArm = FALSE) # NArm = FALSE because it keeps ESV not taxonomically assigned at the species level

## Create a column with the exact sequence of nucleotides 
tax_table(ps.nosd.log.sp) <- cbind(tax_table(ps.nosd.log.sp), ASV=taxa_names(ps.nosd.log.sp))

## Label changed with ESV prefix
taxa_names(ps.nosd.log.sp) <- paste0("ESV", seq(ntaxa(ps.nosd.log.sp)))
ntaxa(ps.nosd.log.sp) # 781 taxa

# Extract phyloseq S4 object in 1 df
mdf = psmelt(ps.nosd.log.sp)


## Exploration ####
## Sequencing depth per samples
depth.ps <- plot_bar(ps, x = "ID", y = "Abundance", fill = "Kingdom") +
  geom_bar(stat="identity", color = "black") +
  labs(x = "Sample ID", y = "Sequence abundance (untransformed)", title = "Distribution sequencing depth per sample")
depth.ps
# ggsave("ITS/depth.ps.pdf", width = 8, height = 8)

## reads per kingdom 
depth.kingdom <- plot_bar(ps, fill = "Kingdom", x = "myco", y = "Abundance") +
  facet_grid(~horiz) +
  geom_bar(stat="identity", color = "black") +
  theme(legend.text = element_text(size=8))
# ggsave("ITS/depth.kingdom.pdf", width = 6, height = 6)

## Check distribution of the lengths of the sequences after filtering
mdf$seq_length <- stri_length(mdf$ASV)

seq_length <- ggplot(mdf, aes(y = Abundance, x = seq_length, color = Phylum)) +
  geom_point(alpha = .4) +
  labs(x="Sequence length", y =  "Sequence abundance (Shifted log)")
seq_length
# ggsave("ITS/seq_length.pdf", width = 6, height = 6)


## plot ESV-abundance curve after filtering
# calculate the mean and choose taxa
clusterAgg = aggregate(Abundance ~ OTU + Phylum,data=mdf, sum)
# filtering and picking the number to display
clusterAgg = clusterAgg[order(-clusterAgg$Abundance),][1:50,]
ggplot(clusterAgg,aes(x=reorder(OTU,-Abundance),y=Abundance)) +
  geom_point(aes(color=Phylum),size=3) + 
  labs(y =  "Sequence abundance (Shifted log)")
  theme(axis.ticks = element_blank(), axis.text.x = element_blank())

# Values for most abundant species (Shifted log)
clusterAggSpecies = aggregate(Abundance ~ OTU + Genus + Species,data=mdf,sum) # sum or mean
most.abund.50 <- clusterAggSpecies[order(-clusterAggSpecies$Abundance),][1:50,]
#write.csv(most.abund.50, "ITS/most.abund.50.csv")

#plot_bar(ps.fungi, fill = "Phylum") + facet_wrap(horiz~myco, ncol =3, scales = "free_x")

# Alhpa diversity
rich <- plot_richness(ps.fungi.nosd, x="myco", measures= "Observed") + 
  geom_boxplot() +
  facet_wrap(~horiz, nrow = 1) +
  scale_x_discrete(label = c(`AM` = "Maple", `mixed` = "Mixed", `ECM` = "Beech")) +
  labs(x="Forest", y= "Observed richness (number of ASV)") +
  theme(axis.text.x = element_text(angle=45, vjust= 1, hjust=1), panel.grid.major.x = element_blank())
rich %>% remove_geom('point')
#ggsave("ITS/rihness.fungi.nosd.pdf", width = 7, height = 4)
#ggsave("ITS/rihness.fungi.nosd.png", width = 7, height = 4)


#### Ordination ####
ps.fungi.ord <- ps.nosd.log.sp
myco_names <- c(`AM` = "AM plots", `mixed` = "Mixed plots", `ECM` = "EcM plots")

## nMDS with Bray-Curtis dissimilarity
quiet(ord.nmds <- ordinate(ps.fungi.ord, method="NMDS", k = 2, try = 4000, distance = "bray"))

plot_nmds <- plot_ordination(ps.fungi.ord, ord.nmds, type="samples", color = "horiz", shape = "horiz") +
  geom_point(size = 8, alpha = .7) +
  facet_wrap(~myco, labeller = as_labeller(myco_names)) +
  scale_color_manual(values = c('darkgreen', 'sienna4', 'grey1', 'grey60', 'darkorange2')) +
  scale_shape_manual(values = c(19, 15, 17, 18, 20)) +
  labs(x="Axis 1", y= "Axis 2", color = "Horizon", shape = "Horizon") +
  coord_fixed() +
  labs(caption = paste("Stress =", round(ord.nmds$stress, 2))) +
  theme(text = element_text(size=20), strip.text = element_text(size=20)) +
  xlim(c(-1.6,1.6)) + ylim(c(-1.6,1.6))
plot_nmds
#ggsave("ITS/nmds.bray.png", width = 12, height = 12)
#ggsave("ITS/nmds.bray.pdf", width = 12, height = 12)

## nMDS with Sorensen dissimilarity
quiet(ord.nmds.bin <- ordinate(ps.fungi.ord, method="NMDS", k = 2, try = 4000, distance = "bray", binary =TRUE))

plot_nmds_bin <- plot_ordination(ps.fungi.ord, ord.nmds.bin, type="samples", color = "horiz", shape = "horiz") +
  geom_point(size = 8, alpha = .7) +
  facet_wrap(~myco, labeller = as_labeller(myco_names)) +
  scale_color_manual(values = c('darkgreen', 'sienna4', 'grey1', 'grey60', 'darkorange2')) +
  scale_shape_manual(values = c(19, 15, 17, 18, 20)) +
  labs(x="Axis 1", y= "Axis 2", color = "Horizon", shape = "Horizon") +
  coord_fixed() +
  labs(caption = paste("Stress =", round(ord.nmds.bin$stress, 2))) +
  theme(text = element_text(size=20), strip.text = element_text(size=20)) +
  xlim(c(-1.6,1.6)) + ylim(c(-1.6,1.6))
plot_nmds_bin
#ggsave("ITS/nmds.bin.png", width = 12, height = 12)
#ggsave("ITS/nmds.bin.pdf", width = 12, height = 12)


#### Permanova ####
# Data
ps.perm <- ps.nosd.log.sp
metadata <- as(sample_data(ps.fungi), "data.frame")

## horizon and myco effect with interaction
perm.horiz.myco <- adonis(phyloseq::distance(ps.perm, method="bray") ~ horiz * myco, strata = metadata$block.y,
       data = metadata,
       permutations = 99999)
perm.horiz.myco
## Test for the homogeneity of dispersion for Bray data
(betadisp.myco <- permutest(betadisper(phyloseq::distance(ps.perm, method="bray"), metadata$myco), permutations = 99999))
(betadisp.horiz <- permutest(betadisper(phyloseq::distance(ps.perm, method="bray"), metadata$horiz), permutations = 99999))

## Sorensen
perm.horiz.myco.bin <- adonis(phyloseq::distance(ps.perm, method = "bray", binary =TRUE) ~ horiz * myco, strata = metadata$block.y, 
       data = metadata,
       permutations = 99999)
perm.horiz.myco.bin
## Homogeneity of dispersion for binary data
(betadisp.myco.bin <- permutest(betadisper(phyloseq::distance(ps.perm, method = "bray", binary =TRUE), metadata$myco), permutations = 99999)) 
(betadisp.horiz.bin <- permutest(betadisper(phyloseq::distance(ps.perm, method = "bray", binary =TRUE), metadata$horiz), permutations = 99999)) 

## Pairwise comparison for Bray data
pairwise.perm.manova(phyloseq::distance(ps.perm, method = "bray"), metadata$myco, nperm=99999,
                     p.method = "fdr") # Benjamini-Hochberg adjustment method

pairwise.perm.manova(phyloseq::distance(ps.perm, method = "bray"), metadata$horiz, nperm=99999,
                     p.method = "fdr")

## Pairwise comparison for Sorensen data
pairwise.perm.manova(phyloseq::distance(ps.perm, method = "bray", binary =TRUE), metadata$myco, nperm=99999,
                     p.method = "fdr")

pairwise.perm.manova(phyloseq::distance(ps.perm, method = "bray", binary =TRUE), metadata$horiz, nperm=99999,
                     p.method = "fdr")

#### Canonical Analysis (CCA) ####
## CCA
ord.cca <- ordinate(ps.fungi.ord, formula = ps.fungi.ord ~ CN+inorgP+orgP+BrayP+pHCaCl2+ECEC, method="CCA")
#ord.cca <- ordinate(ps.fungi.ord, formula = ps.fungi.ord ~ CN+inorgP+orgP+BrayP+pHCaCl2+Al+Ca+Fe+K+Mn+Na, method="CCA")
evals <- ord.cca$CCA$eig

# Get stat
rs <- RsquareAdj(ord.cca); rs
pv <- anova.cca(ord.cca); pv

# Check variance inflation factors 
vif.cca(ord.cca)

## Plot CCA with arrows for environmental data
# Add the environmental variables as arrows
arrowmat = vegan::scores(ord.cca, display = "bp")
# Add labels, make a data.frame
arrowdf <- data.frame(arrowmat)
rownames(arrowdf) <- c("CN","inorgP","orgP","BrayP","pHCaCl2","ECEC")
arrowdf <- data.frame(labels = c("CN","Pi","Po","BrayP","pHCaCl2","ECEC"), arrowmat)

CA.log.paired.ITS.sp <- plot_ordination(ps.fungi.ord, vegan::scores(ord.cca, scaling=1), type="sites", color = "horiz") + 
  coord_fixed(evals[2] / evals[1]) +
  geom_point(size = 8, alpha = 1) +
  geom_segment(aes(xend = CCA1, yend = CCA2, x = 0, y = 0), size = 0.8, data = arrowdf, color = "black", arrow = arrow(length = unit(0.025, "npc"))) + 
  geom_text(aes(x = CCA1*1.3, y = CCA2*1.3, color = NULL, label = labels), size = 6, data = arrowdf) +
  facet_wrap(~myco, labeller = as_labeller(myco_names)) +
  theme(text = element_text(size=20), strip.text = element_text(size=20)) +
  scale_color_manual(values = c('darkgreen', 'sienna4', 'grey1', 'grey60', 'darkorange2')) +
  labs(x="CC axis 1", y= "CC axis 2", color = "Horizon") +
  #labs(caption = paste("Adjusted R2 =", round(rs$adj.r.squared*100, 2), "%")) +
  xlim(c(-1.4,1.4)) + ylim(c(-1.4,1.4))
CA.log.paired.ITS.sp

#### RDA and Variance paritioning ####
# matrices Y
ps.vp <- ps.nosd.log.sp
metadata <- sample_data(ps.vp)
rownames(metadata) <- metadata$ID

# Data
env.chem <- as.matrix(metadata[,c(chem <- c("CN", "inorgP", "orgP", "BrayP", "pHCaCl2", "ECEC"))])
horiz <- as(metadata[,"horiz"], "data.frame")
myco <- as(metadata[,"sp"], "data.frame")
block <- as(metadata[,"block.y"], "data.frame"); block$block.y <- as.factor(block$block.y)

## Matrix X
seqtab.dbray <- phyloseq::distance(ps.vp, method="bray")

# Converting Categorical Columns into Multiple Binary Columns and trimming latent variables
env.horiz <- as.matrix(model.matrix(~horiz -1, data=horiz))
env.myco <- as.matrix(model.matrix(~sp -1, data=myco))
env.block <- as.matrix(model.matrix(~block.y -1, data=block))

# deleting redundant variables
env.horiz <- env.horiz[,1:4]
env.myco <- env.myco[,1:2]
env.block <- env.block[,1:4]

# Appying mutliple partial RDA to partitionate the variance with varpart
var.part.all <- varpart(seqtab.dbray, env.chem, env.horiz, env.myco) # env.block
plot(var.part.all, bg = 1:3, digits = 2, Xnames = c("Soil Chemistry","Horizon","Mycorrhizal Type"), cex = 1.5) # cutoff = 0.001
var.part.all$part$fract$Adj.R.square

# with capscale
rda.all <- capscale(formula = seqtab.dbray ~ env.chem + env.horiz + env.myco, distance = NULL) # Inertia is variance, proportion is r.squared (in vegan)
RsquareAdj(rda.all) # same than var.part, 0.3333739
anova(rda.all) #0.001

## Use function ‘capscale’ (db-RDA) to test significance of fractions of interest.
## Semi-Partial correlation
prda.X1 <- capscale(formula = seqtab.dbray ~ env.chem + Condition(env.horiz + env.myco))
r.prda.X1 <- RsquareAdj(prda.X1)
a.prda.X1 <- anova(prda.X1)

prda.X2 <- capscale(formula = seqtab.dbray ~ env.block + Condition(env.chem + env.horiz + env.myco))
r.prda.X2 <- RsquareAdj(prda.X2)
a.prda.X2 <- anova(prda.X2)

prda.X3 <- capscale(formula = seqtab.dbray ~ env.horiz + Condition(env.chem + env.myco))
r.prda.X3 <- RsquareAdj(prda.X3)
a.prda.X3 <- anova(prda.X3)

prda.X4 <- capscale(formula = seqtab.dbray ~ env.myco + Condition(env.horiz + env.horiz))
r.prda.X4 <- RsquareAdj(prda.X4)
a.prda.X4 <- anova(prda.X4)

# Partial correlation
rda.X1 <- capscale(formula = seqtab.dbray ~ env.chem)
r.rda.X1 <- RsquareAdj(rda.X1)
a.rda.X1 <- anova(rda.X1)

rda.X2 <- capscale(formula = seqtab.dbray ~ env.block)
r.rda.X2 <- RsquareAdj(rda.X2)
a.rda.X2 <- anova(rda.X2)

rda.X3 <- capscale(formula = seqtab.dbray ~ env.horiz)
r.rda.X3 <- RsquareAdj(rda.X3)
a.rda.X3 <- anova(rda.X3)

rda.X4 <- capscale(formula = seqtab.dbray ~ env.myco)
r.rda.X4 <- RsquareAdj(rda.X4)
a.rda.X4 <- anova(rda.X4)

# Put results in a table
vp.results <- data.frame(X = c("Chemistry", "Block", "Horizon", "Stand Type"),
           RDA = c(r.rda.X1$adj.r.squared, r.rda.X2$adj.r.squared, r.rda.X3$adj.r.squared, r.rda.X4$adj.r.squared),
           Pr = c(a.rda.X1$`Pr(>F)`[1], a.rda.X2$`Pr(>F)`[1], a.rda.X3$`Pr(>F)`[1], a.rda.X4$`Pr(>F)`[1]),
           pRDA = c(r.prda.X1$adj.r.squared, r.prda.X2$adj.r.squared, r.prda.X3$adj.r.squared, r.prda.X4$adj.r.squared),
           Pr = c(a.prda.X1$`Pr(>F)`[1], a.prda.X2$`Pr(>F)`[1], a.prda.X3$`Pr(>F)`[1], a.prda.X4$`Pr(>F)`[1]))
vp.results[,c(2,4)] <- vp.results[,c(2,4)]*100 # r2adj as pourcentage
vp.results[,c(2,4)] <- round(vp.results[,c(2,4)], 1)
vp.results
#write.csv(vp.results, "ITS/vp.results.csv")

## Plot RDA for soil chemistry only
# Check variance inflation factors 
vif.cca(rda.X1)
evals <- rda.X1$CCA$eig/sum(rda.X1$CCA$eig)
# Add the environmental variables as arrows
arrowmat = vegan::scores(rda.X1, display = "bp", scaling = 1)
# Add labels, make a data.frame
arrowdf <- data.frame(arrowmat)
rownames(arrowdf) <- c("CN","inorgP","orgP","BrayP","pHCaCl2","ECEC")
arrowdf <- data.frame(labels = c("C:N","Pi","Po","labileP","pH","cations"), arrowmat)

plot_ordination(ps.nosd.log.sp, vegan::scores(rda.X1, scaling = 1), type="sites", color = "horiz") + 
  coord_fixed() +
  geom_point(size = 7, alpha = 1) +
  scale_color_manual(name = "Horizon", values = c('darkgreen', 'sienna4', 'grey1', 'grey60', 'darkorange2')) +
  geom_segment(aes(xend = CAP1, yend = CAP2, x = 0, y = 0), size = .7, data = arrowdf, color = "black", arrow = arrow(length = unit(0.05, "npc"))) +
  geom_text(aes(x = CAP1*1.35, y = CAP2*1.35, color = NULL, label = labels), size = 6, data = arrowdf) +
  #facet_wrap(~myco, labeller = as_labeller(myco_names)) +
  theme(text = element_text(size=18), strip.text = element_text(size=18)) +
  labs(x = paste("Constrained axis 1 (",round(evals[1]*100, 0),"%)"), y = paste("Constrained axis 2 (",round(evals[2]*100, 0),"%)"))
# ggsave("1.ITS/rda.chem.pdf", width = 8, height = 8)

## Plot RDA for myco data
vif.cca(rda.X4)
evals <- rda.X4$CCA$eig/sum(rda.X1$CCA$eig)
# Add the environmental variables as arrows
arrowmat = vegan::scores(rda.X4, display = "bp", scaling = 1)
# Add labels, make a data.frame
arrowdf <- data.frame(arrowmat)
rownames(arrowdf) <- c("AM","EcM")
arrowdf <- data.frame(labels = c("AM","EcM"), arrowmat)

plot_ordination(ps.nosd.log.sp, vegan::scores(rda.X4, scaling = 1), type="sites",  shape = "myco", color = "myco") + 
  geom_point(size = 6, alpha = .7) +
  scale_shape_discrete(name = "Plots", labels = c(`AM` = "AM", `mixed` = "Mixed", `ECM` = "EcM")) +
  scale_color_discrete(name = "Plots", labels = c(`AM` = "AM", `mixed` = "Mixed", `ECM` = "EcM")) +
  coord_fixed(sqrt(evals[2] / evals[1])) +
  theme(text = element_text(size=18)) +
  labs(color = "Forest", shape = "Forest", x = paste("Constrained axis 1 (",round(evals[1]*100, 0),"%)"), y = paste("Constrained axis 2 (",round(evals[2]*100, 0),"%)"))
#ggsave("ITS/rda.myco.png", width = 8, height = 8)
#ggsave("ITS/rda.myco.pdf", width = 8, height = 8)

#### Funguild ####
guild <- read.table("ITS/ps.nosd.log.sp.29.05.dynamic_all_02.02.2019.minboot50.guilds.txt", header = TRUE, sep = "\t") # /!\ Has to be created with the same mdf

# Merge guild info with abundance data
ITS.guild  <- left_join(mdf, guild, by = c("OTU" = "OTU_ID"))

# Simplify and relevel dataframe
ITS.guild$Confidence.Ranking = factor(ITS.guild$Confidence.Ranking, levels = c('Highly Probable', 'Probable', 'Possible', '-'))
ITS.guild$Trophic.Mode <- gsub("Pathogen-Saprotroph-Symbiotroph" , "Others", ITS.guild$Trophic.Mode)
ITS.guild$Trophic.Mode <- gsub("Pathotroph-Saprotroph-Symbiotroph" , "Others", ITS.guild$Trophic.Mode)
ITS.guild$Trophic.Mode <- gsub("Pathotroph-Saprotroph" , "Patho_Sapro", ITS.guild$Trophic.Mode)
ITS.guild$Trophic.Mode <- gsub("Pathotroph-Symbiotroph" , "Others", ITS.guild$Trophic.Mode)
ITS.guild$Trophic.Mode <- gsub("Saprotroph-Symbiotroph" , "Sapro_Symbio", ITS.guild$Trophic.Mode)
ITS.guild$Trophic.Mode <- gsub("-" , "Others", ITS.guild$Trophic.Mode)
ITS.guild$Trophic.Mode <- factor(ITS.guild$Trophic.Mode, levels = c('Others', 'Pathotroph', 'Patho_Sapro', 'Saprotroph', 'Sapro_Symbio', 'Symbiotroph'))
ITS.guild$myco <- gsub("AM", "Maple", ITS.guild$myco)
ITS.guild$myco <- gsub("ECM", "Beech", ITS.guild$myco)
ITS.guild$myco <- gsub("mixed", "Mixed", ITS.guild$myco)
ITS.guild$myco = factor(ITS.guild$myco, levels = c('Maple', 'Mixed', 'Beech'))

## Plot trophic mode abundance in samples
myco_names <- c(`Maple` = "AM plots", `Mixed` = "Mixed plots", `Beech` = "EcM plots", "L" = "L", "F" = "F", "H" = "H", "Ae" = "Ae", "B" = "B")

## Plot with subset of AM, EcM and Sapro only
ITS.guild$myco <- gsub("AM", "Maple", ITS.guild$myco)
ITS.guild$myco <- gsub("ECM", "Beech", ITS.guild$myco)
ITS.guild$myco = factor(ITS.guild$myco, levels = c('Maple', 'Mixed', 'Beech'))

ITS.sapro <- ITS.guild[ITS.guild$Trophic.Mode == "Saprotroph", ]
ITS.ecm <- ITS.guild[ITS.guild$Guild == "Ectomycorrhizal", ]
ITS.am <- ITS.guild[ITS.guild$Guild == "Arbuscular Mycorrhizal", ]

ggplot(ITS.sapro, aes(x = Trophic.Mode, y=Abundance, fill = Confidence.Ranking)) +
  geom_bar(stat = "summary", fun.y = "mean") +
  scale_fill_grey() +
  facet_grid(horiz~myco, labeller = as_labeller(myco_names)) +
  theme(axis.text.x = element_text(size = 0), text = element_text(size=20), plot.title = element_text(size=20)) +
  labs(x = "Saprotroph", y = "Sequence abundance (Shifted log)", title = "")

## LME on guilds and tukey adjustment for SE
# Sapro
ITS.sapro.sample <- ITS.sapro  %>% 
  group_by(Sample) %>% 
  summarise(Abundance = sum(Abundance), ECM_perc = unique(ECM_perc), horiz = unique(horiz), myco = unique(myco), block = unique(block.y))

sapro.lme <- lme(Abundance ~ myco*horiz, random = ~ 1 | block,  weights = varIdent(form = ~ 1 | horiz), data=ITS.sapro.sample)
sapro.tuckey <- emmeans(sapro.lme, pairwise ~ horiz*myco, adjust = "tukey")
sapro.mult <- CLD(sapro.tuckey,alpha=0.05,Letters=letters, adjust="tukey")

sapro.plot <- ggplot(sapro.mult, aes(x = myco, y = emmean, fill = horiz)) + 
  geom_bar(stat = "identity", width = .7, colour="black", show.legend = FALSE) +
  scale_fill_manual(values = c('darkgreen', 'sienna4', 'grey1', 'grey60', 'darkorange2')) +
  geom_errorbar(aes(ymin = emmean, ymax = emmean+SE), width = 0.2) +
  facet_wrap(~horiz,  nrow = 5, scales = "fixed", labeller = as_labeller(myco_names)) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1), text = element_text(size=20), strip.text = element_blank(), strip.background = element_blank(), plot.title = element_text(hjust = 0.5, size=20), panel.spacing = unit(1, "lines"), panel.grid.major.x = element_blank()) +
  labs(x = "Forest", y = "", title = "Saprotrophic")

# ECM
ITS.ecm.sample <- ITS.ecm  %>% 
  group_by(Sample) %>% 
  summarise(Abundance = sum(Abundance), ECM_perc = unique(ECM_perc), horiz = unique(horiz), myco = unique(myco), block = unique(block.y))

ecm.lme <- lme(Abundance ~ myco*horiz, random = ~ 1 | block,  weights = varIdent(form = ~ 1 | horiz), data=ITS.ecm.sample)
ecm.tuckey <- emmeans(ecm.lme, pairwise ~ myco*horiz, adjust = "tukey")
ecm.mult <- CLD(ecm.tuckey,alpha=0.05,Letters=letters, adjust="tukey")
ecm.plot <- ggplot(ecm.mult, aes(x = myco, y = emmean, fill = horiz)) + 
  geom_bar(stat = "identity", width = .7, colour="black", show.legend = FALSE) +
  scale_fill_manual(values = c('darkgreen', 'sienna4', 'grey1', 'grey60', 'darkorange2')) +
  geom_errorbar(aes(ymin = emmean, ymax = emmean+SE), width = 0.2) +
  facet_wrap(~horiz,  ncol = 1, scales = "fixed", strip.position = "right", labeller = as_labeller(myco_names)) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1), text = element_text(size=20), strip.text.y = element_text(angle = 0), strip.background.y = element_blank(), plot.title = element_text(hjust = 0.5, size=20), panel.spacing = unit(1, "lines"), panel.grid.major.x = element_blank()) +
  labs(x = "", y = "", title = "EcM")

# AM
ITS.am.sample <- ITS.am  %>% 
  group_by(Sample) %>% 
  summarise(Abundance = sum(Abundance), ECM_perc = unique(ECM_perc), horiz = unique(horiz), myco = unique(myco), block = unique(block.y))

am.lme <- lme(Abundance ~ myco*horiz, random = ~ 1 | block,  weights = varIdent(form = ~ 1 | horiz), data=ITS.am.sample)
am.tuckey <- emmeans(am.lme, pairwise ~ myco*horiz, adjust = "tukey")
am.mult <- CLD(am.tuckey,alpha=0.05,Letters=letters, adjust="tukey")
am.plot <- ggplot(am.mult, aes(x = myco, y = emmean, fill = horiz)) + 
  geom_bar(stat = "identity", width = .7, colour="black", show.legend = FALSE) +
  scale_fill_manual(values = c('darkgreen', 'sienna4', 'grey1', 'grey60', 'darkorange2')) +
  geom_errorbar(aes(ymin = emmean, ymax = emmean+SE), width = 0.2) +
  facet_wrap(~horiz,  nrow = 5, scales = "fixed", labeller = as_labeller(myco_names)) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1), text = element_text(size=20), strip.text = element_blank(), strip.background = element_blank(), plot.title = element_text(hjust = 0.5, size=20), panel.spacing = unit(1, "lines"), panel.grid.major.x = element_blank()) +
  labs(x = "", y = "Sequence mean abundance (Shifted log)", title = "AM")

guild.plot <- grid.arrange(am.plot, sapro.plot, ecm.plot, nrow = 1, widths = c(1,1,1.1))
guild.plot
# ggsave("ITS/guild.am.ecm.sapro.SE.pdf", guild.plot, width = 8, height = 10)
# ggsave("ITS/guild.am.ecm.sapro.SE.png", guild.plot, width = 8, height = 10)


### Indicator species #### 
iva <- indval(otu_table(ps.nosd.log.sp), sample_data(ps.nosd.log.sp)$myco_horiz)

# Table of the significant indicator species
gr <- iva$maxcls[iva$pval <= 0.05]
iv <- iva$indcls[iva$pval <= 0.05]
pv <- iva$pval[iva$pval <= 0.05]
fr <- apply(otu_table(ps.nosd.log.sp) > 0, 2, sum)[iva$pval <= 0.05]
fidg <- data.frame(group=gr, indval=iv, pvalue=pv, freq=fr)
fidg <- fidg[order(fidg$group, -fidg$indval),]
fidg$OTU <- rownames(fidg)

groups <- data.frame(group = as.integer(rep(1:15)),
                     myco_horiz=levels(sample_data(ps.nosd.log.sp)$myco_horiz))

tax <- as.data.frame(as(tax_table(ps.nosd.log.sp), "matrix"))
tax$OTU <- rownames(tax)

indic.species <- fidg %>%
  left_join(groups, by = 'group') %>%
  left_join(tax, by = 'OTU')

indic.max.species <- indic.species %>% 
  group_by(group) %>% 
  top_n(1, indval)
indic.max.species[,c(2:13)]
#write.csv(indic.max.species, "ITS/indic.max.species.csv")


