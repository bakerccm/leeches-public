# R code to draw Fig 1b (3D map of Ailaoshan study site)

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
library("ggspatial")

# get data

    load(file=here("rdata", "Ailaoshan_gis.rdata"))

# make 3D plots

    # set window size
        #par3d(windowRect = c(200, 000, 1200, 1000))

    # user matrix for setting viewpoint
        # um <- par3d()$userMatrix
        um <- matrix(data = c(0.97, 0.23, 0.01, 0,
            -0.13, 0.51, 0.85, 0,
            0.19, -0.83, 0.52, 0,
            0, 0, 0, 1), nrow=4, ncol=4, byrow=TRUE)
        view3d(userMatrix = um)

    # red on colour terrain
        plot3D(ailaoshan.srtm.reproj, drape=ailaoshan.srtm.reproj, zfac=3, adjust=FALSE)
        plot3D(ailaoshan.polygons.rast, drape=ailaoshan.polygons.rast, zfac=3, adjust=FALSE, col="red")
        rgl.snapshot(here("figures", "Fig1b_Ailaoshan_site_colour.png"), fmt = "png")
    
    # red on grey
        rgl.clear()
        plot3D(ailaoshan.srtm.reproj, drape=ailaoshan.srtm.reproj, zfac=3, adjust=FALSE, col=gray(seq(0.8, 1, 0.02)))
        plot3D(ailaoshan.polygons.rast, drape=ailaoshan.polygons.rast, zfac=3, adjust=FALSE, col="red")
        rgl.snapshot(here("figures", "Fig1b_Ailaoshan_site_grey.png"), fmt = "png")
    
    # make movie
        #movie3d(spin3d(axis = c(0, 0, 1), rpm = 10), duration = 6, fps = 25, movie="ailaoshan_site", dir = "figures")
