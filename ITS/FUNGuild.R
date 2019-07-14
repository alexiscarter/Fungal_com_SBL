## Extract taxonomy information
tax <- as.data.frame(as(tax_table(ps.nosd.log.sp), "matrix"))
tax$OTU <- rownames(tax)
colnames(tax)
fg <- as.matrix(data.frame(OTU_ID = tax$OTU,
                           taxonomy = paste(tax$Phylum, tax$Class, tax$Order, tax$Family, tax$Genus, tax$Species, sep=";")))
## save matrix
# write.csv(fg, file = "1.ITS/ps.nosd.log.sp.29.05.dynamic_all_02.02.2019.minboot50.csv", row.names = FALSE)
# get rid of the "

## Submit matrix to http://www.stbates.org/guilds/app.php
## Download OTU file
## Load it
guild <- read.table("ITS/ps.nosd.log.sp.13.05.dynamic_all_02.02.2019.minboot50.guilds.txt", header = TRUE, sep = "\t") # /!\ Has to be created with the same mdf

