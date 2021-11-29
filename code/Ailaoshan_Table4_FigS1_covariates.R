# R code to generate Table 4 and Fig S1 to show environmental covariates

library("here")
library("sf")
library("raster")
library("sp")
library("tidyverse")

# get data

    load(file=here("rdata","Ailaoshan_environmental.rdata"))

    load(file=here("rdata","Ailaoshan_gis.rdata"))

# prepare data

    selected.env.data <- env.data %>%
        # get covariates of interest
            select(Polygon_ID, elevation_median, elevation_median, tpi_median, distance_to_road_median, distance_to_stream_median, distance_to_nature_reserve_boundary) %>%
        # rename covariates
            rename(elevation = elevation_median, TPI = tpi_median, road = distance_to_road_median, stream = distance_to_stream_median,
                reserve = distance_to_nature_reserve_boundary) %>%
        # filter remove Polygon_IDs that are just a ranger ID (these don't have environmental data anyway)
            filter(substr(Polygon_ID,1,1) != "R")

# attach to spatial polygons for mapping

    ailaoshan.polygons.env <- ailaoshan.polygons %>% mutate(Polygon_ID = as.character(Polygon_ID)) %>% left_join(selected.env.data, by="Polygon_ID")

# draw maps showing env covariates

    ailaoshan.polygons.env %>%
        ggplot() + geom_sf(aes(fill = elevation), size = 0.1, color = "black") +
            theme(legend.position = c(0,0), legend.justification = c(0,0), legend.background = element_rect(fill = NA)) +
            coord_sf() + labs(fill = "elevation (m)")
    ggsave(here("figures","FigS1a_elevation.pdf"), width=4, height=4)

    ailaoshan.polygons.env %>%
        ggplot() + geom_sf(aes(fill = TPI), size = 0.1, color = "black") +
            theme(legend.position = c(0,0), legend.justification = c(0,0), legend.background = element_rect(fill = NA)) +
            coord_sf()
    ggsave(here("figures","FigS1c_TPI.pdf"), width=4, height=4)

    ailaoshan.polygons.env %>%
        ggplot() + geom_sf(aes(fill = road), size = 0.1, color = "black") +
            theme(legend.position = c(0,0), legend.justification = c(0,0), legend.background = element_rect(fill = NA)) +
            coord_sf() + labs(fill = "road (m)")
    ggsave(here("figures","FigS1e_road.pdf"), width=4, height=4)

    ailaoshan.polygons.env %>%
        ggplot() + geom_sf(aes(fill = stream), size = 0.1, color = "black") +
            theme(legend.position = c(0,0), legend.justification = c(0,0), legend.background = element_rect(fill = NA)) +
            coord_sf() + labs(fill = "stream (m)")
    ggsave(here("figures","FigS1g_stream.pdf"), width=4, height=4)

    ailaoshan.polygons.env %>%
        ggplot() + geom_sf(aes(fill = reserve), size = 0.1, color = "black") +
            theme(legend.position = c(0,0), legend.justification = c(0,0), legend.background = element_rect(fill = NA)) +
            coord_sf() + labs(fill = "reserve (m)")
    ggsave(here("figures","FigS1i_reserve.pdf"), width=4, height=4)

# histograms for environmental covariates

    selected.env.data %>% ggplot(aes(x = elevation)) + geom_histogram() +
        labs(x = "median elevation (m)")
    ggsave(here("figures","FigS1b_elevation.pdf"), width=3.2, height=4)

    selected.env.data %>% ggplot(aes(x = TPI)) + geom_histogram() +
        labs(x = "median TPI") +
        scale_y_continuous(breaks = seq(0,50,10), labels = c(seq(0,40,10), "")) # omit top tick mark label to avoid overlap with subfigure label
    ggsave(here("figures","FigS1d_TPI.pdf"), width=3.2, height=4)

    selected.env.data %>% ggplot(aes(x = road)) + geom_histogram() +
        labs(x = "median distance to road (m)")
    ggsave(here("figures","FigS1f_road.pdf"), width=3.2, height=4)

    selected.env.data %>% ggplot(aes(x = stream)) + geom_histogram() +
        labs(x = "median distance to stream (m)")
    ggsave(here("figures","FigS1h_stream.pdf"), width=3.2, height=4)

    selected.env.data %>% ggplot(aes(x = reserve)) + geom_histogram() +
        labs(x = "centroid distance to reserve boundary (m)")
    ggsave(here("figures","FigS1j_reserve.pdf"), width=3.2, height=4)

# summary table for environmental covariates
# Table 1 in manuscript

    env.summary <- selected.env.data %>%
        pivot_longer(-Polygon_ID, names_to = "covariate") %>%
        # summarize covariates
            group_by(covariate) %>%
            summarize(mean = mean(value), sd = sd(value), min = min(value), max = max(value), .groups = "drop") %>%
        # reorder rows for publication
            mutate(covariate_order = recode(covariate, elevation="A", TPI="B", road="C", stream="D", reserve="E")) %>%
            arrange(covariate_order) %>% select(-covariate_order)

    write_csv(env.summary, path = here("tables", "Table4_environmental_covariates_summary.csv"))
