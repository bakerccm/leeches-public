# R code to process GIS files mostly in preparation for making figures of study site

# see e.g.https://datacarpentry.org/r-raster-vector-geospatial/06-vector-open-shapefile-in-r/
# and https://www.ecologi.st/post/spatial-data-in-r-2-a-practical-example/

# Note. download DEM data from USGS earth explorer:
# https://earthexplorer.usgs.gov

library("here")
library("sf")
library("raster")
library("rgdal")
library("fasterize")
library("rasterVis")
library("tidyverse")

# get data

    # Ailaoshan polygons
        ailaoshan.polygons <- st_read(here("data", "gis", "ailaoshan_polygons", "study_area.shp"))

    # SRTM data
        ailaoshan.srtm <- list()
        ailaoshan.srtm[[1]] <- raster(here("data", "gis", "srtm_data", "n23_e100_1arc_v3.tif"))
        ailaoshan.srtm[[2]] <- raster(here("data", "gis", "srtm_data", "n23_e101_1arc_v3.tif"))
        ailaoshan.srtm[[3]] <- raster(here("data", "gis", "srtm_data", "n24_e100_1arc_v3.tif"))
        ailaoshan.srtm[[4]] <- raster(here("data", "gis", "srtm_data", "n24_e101_1arc_v3.tif"))
        ailaoshan.srtm[[5]] <- raster(here("data", "gis", "srtm_data", "n25_e100_1arc_v3.tif"))
        ailaoshan.srtm[[6]] <- raster(here("data", "gis", "srtm_data", "n25_e101_1arc_v3.tif"))

# add Polygon_ID to ailaoshan.polygons

    # OBJECTID in ailaoshan.polygons is incremented by 1 compared to Polygon_ID but the datasets otherwise appear to be the same:
        # load(file="../Ailaoshan_environmental.rdata")
        # env.data$Polygon_ID == as.numeric(ailaoshan.polygons$OBJECTID) - 1
        # round(env.data$shape_perimeter,2) == round(as.numeric(ailaoshan.polygons$Shape_Leng),2)
        # round(env.data$shape_area_ha,2) == round(as.numeric(ailaoshan.polygons$Shape_Area),2)
        # etc

    ailaoshan.polygons <- data.frame(ailaoshan.polygons, Polygon_ID = as.numeric(ailaoshan.polygons$OBJECTID) - 1) %>%
        select(OBJECTID, Polygon_ID, everything()) %>% st_sf()

# examine Ailaoshan polygons

    # st_geometry_type(ailaoshan.polygons)
    # st_crs(ailaoshan.polygons)
    # extent(ailaoshan.polygons)

    # ailaoshan.polygons %>% ggplot() + geom_sf(size = 0.1, color = "black", fill = "burlywood2") + ggtitle("Ailaoshan") + coord_sf()

# process SRTM data

    # merge srtm tiles into one big raster
        ailaoshan.srtm <- do.call(merge, ailaoshan.srtm)

    # crop to a little outside Ailaoshan nature reserve
        ailaoshan.srtm.crop <- crop(ailaoshan.srtm, extent(100.6, 101.6, 23.9, 25.1))

    # reproject and crop SRTM data}
    # R typically uses the PROJ.4 conventions for cartographic projections (or coordinate reference systems - CRS).
    # Check out http://proj4js.org or http://spatialreference.org/ or google for the “proj4string” for various coordinate reference systems.

    # note projection does not match that for Ailaoshan polygons:

        st_crs(ailaoshan.polygons)
        proj4string(ailaoshan.srtm.crop)

        extent(ailaoshan.polygons)
        extent(ailaoshan.srtm.crop)

    # so reproject data

        ailaoshan.polygons.crs <- st_crs(ailaoshan.polygons)$proj4string
        ailaoshan.srtm.reproj <- projectRaster(ailaoshan.srtm.crop, crs = CRS(ailaoshan.polygons.crs))

    # now projection should match Ailaoshan polygons:

        st_crs(ailaoshan.polygons)
        proj4string(ailaoshan.srtm.reproj)

        extent(ailaoshan.polygons)
        extent(ailaoshan.srtm.reproj)

    #### rasterize ####
        ailaoshan.polygons.rast <- fasterize(ailaoshan.polygons, ailaoshan.srtm.reproj)
        ailaoshan.polygons.rast <- ailaoshan.polygons.rast * ailaoshan.srtm.reproj + 1

# save data to file

    save(ailaoshan.polygons, ailaoshan.polygons.crs, ailaoshan.polygons.rast, ailaoshan.srtm, ailaoshan.srtm.crop, ailaoshan.srtm.reproj,
        file=here("rdata","Ailaoshan_gis.rdata"))
