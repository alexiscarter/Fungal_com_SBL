#### Map of the locations of the sites ####

#library
library(ggplot2)
library(ggrepel)
library(ggmap)
library(ggsn)

# Need API key to use get_map, not free access anymore

# Geographic coordinates x and y from the plots data frame
load("data/plots.rda")
plots <- plots[1:15,]

plots$plot <- gsub("M_FS_1/2", "Mixed", plots$plot)
plots$plot <- gsub("M_AS_FS", "Maple", plots$plot)
plots$plot <- gsub("M_FG_FS", "Beech", plots$plot)

# correct longitudes
plots$long <- plots$long * -1

# get the map
boite <- make_bbox(lon=-73.995, lat=45.989)
SBL <- get_map(location=boite, maptype="satellite", source="google", zoom=14)
sbl <- ggmap(SBL)

## For the deocmposition experiment
#plots <- plots[plots$myco %in% c('AM', 'ECM'),]

# add points
sbl2 <- sbl + 
  geom_point(data = plots, aes(long, lat), size = 2, color = 'lightgrey') +
  xlab("Longitude (°)") +
  ylab("Latitude (°)") +
  geom_label_repel(data = plots, aes(long, lat, fill = block, label = plot),
    fontface = 'bold', color = 'lightgrey',
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.5, "lines"),
    segment.color = 'lightgrey',
    show.legend = FALSE, size = 4) +
  #scale_shape_discrete((name="Mycorrhizal\nType"), labels=c("AM", "EcM", "Mixed")) +
  #scale_fill_discrete(name="Block") +
  scalebar(sbl2, dist = .5, dist_unit = "km",transform = FALSE, model = "WGS84")

