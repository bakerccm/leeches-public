# R code to draw Fig 1a (continent-scale map to show location of study site)

# see e.g.https://datacarpentry.org/r-raster-vector-geospatial/06-vector-open-shapefile-in-r/
# and https://www.ecologi.st/post/spatial-data-in-r-2-a-practical-example/

# Note. download DEM data from USGS earth explorer:
# https://earthexplorer.usgs.gov

library("here")
library("sf")
library("gplots")
library("ggrepel")
library("readxl")
library("viridis")
library("raster")
library("fasterize")
library("rasterVis")
library("rgl")
library("magick")
library("tidyverse")
library("egg")

# get data

    load(file=here("rdata", "Ailaoshan_gis.rdata"))

# get world map

    library("rnaturalearth")
    library("rnaturalearthdata")
    library("rgeos")
    # library("rnaturalearthhires") # load if you want to use scale = "large" below

    world <- ne_countries(scale = "medium", returnclass = "sf")

# draw map

    ggplot() +
        geom_sf(data = world, fill="grey", color=NA) +
        geom_sf(data = ailaoshan.polygons[ailaoshan.polygons$Polygon_ID == 9,], fill="red", color="red") +
        geom_sf(data = world, fill=NA, color="black") +
        # geographic limits of map:
        coord_sf(xlim = c(88.5, 109.5), ylim = c(10, 31), expand = FALSE)

    ggsave(here("figures", "Fig1a_Ailaoshan_map_continent.pdf"), device = "pdf", width=4, height=4)
