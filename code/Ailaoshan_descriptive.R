# R code for descriptive analyses of Ailaoshan data

library("here")
library("tidyverse")
library("vegan")
library("segmented")

# read in raw data and MCMC output

    load(file=here("rdata", "Ailaoshan_OTU_table.rdata"))

# some summary counts

    leech %>% group_by(dataset) %>%
        summarize(
            number.PolygonIDs = length(unique(Polygon_ID)),
            number.LabIDs = length(unique(Lab_ID)),
            number.OTUs = length(unique(OTU)),
            .groups = "drop")

# Ranger_IDs and Polygon_IDs

    # most Ranger_IDs are associated with one polygon, but some are associated with >1.

        leech %>% group_by(dataset, Ranger_ID) %>%
            summarize(number.Polygon_IDs = length(unique(Polygon_ID)), .groups = "drop") %>%
            select(dataset, number.Polygon_IDs) %>% table()

    # Note that some Ranger_IDs that are associated with 1 Polygon_ID were actually associated with
    # no polygons originally (this field was NA) and the Polygon_ID is just the Ranger_ID.

        leech %>% filter(is.na(Polygon_ID_original)) %>%
            group_by(dataset, Ranger_ID) %>%
            summarize(number.Polygon_IDs = length(unique(Polygon_ID)), .groups = "drop") %>%
            select(dataset, number.Polygon_IDs) %>% table()

# leeches and Lab_IDs

    # number of leeches in each dataset
    # note these are just the ones that made it through to the OTU table
    # LSU = 23774, SSU = 26775

    leech %>%
        distinct(dataset, Lab_ID, leech_qty) %>%
        group_by(dataset) %>%
        summarize(total.leeches = sum(leech_qty), .groups = "drop")

    # some Lab_IDs have zero species observed OTU richness (because they only matched to humans -- the one
    # Lab_ID that had no OTUs in the LSU dataset after excluding humans was taken out altogether)

        leech %>% group_by(dataset, Lab_ID) %>%
            summarize(observed.richness = sum(reads > 0), .groups = "drop_last") %>%
            select(dataset, observed.richness) %>% table()

        leech %>% group_by(dataset, Lab_ID) %>%
            summarize(observed.richness = sum(reads > 0), .groups = "drop_last") %>%
            filter(observed.richness > 0) %>% tally() # LSU: 553, SSU: 675 Lab_IDs with at least one non-human OTU

    # number of leeches per Lab_ID

        leech %>%
            distinct(dataset, Lab_ID, leech_qty) %>%
            group_by(dataset) %>%
            summarize(mean.leech_qty = mean(leech_qty), # LSU: 36.40735, SSU 36.18243
                sd.leech_qty = sd(leech_qty), # LSU: 18.91909, SSU 18.42684
                median.leech_qty = median(leech_qty), # LSU 35, SSU 34
                min.leech_qty = min(leech_qty),
                max.leech_qty = max(leech_qty),
                .groups = "drop")

        leech %>%
            distinct(dataset, Lab_ID, leech_qty) %>%
            ggplot(aes(x=leech_qty)) + geom_histogram(bins = 10) + facet_wrap(~dataset) # positive skew in both datasets

# Ranger_IDs and Lab_IDs

    # number of rangers

        leech %>% group_by(dataset) %>%
            summarize(number.RangerIDs = length(unique(Ranger_ID)), .groups = "drop") # LSU: 125 rangers, SSU: 126 rangers

    # how many Lab_IDs were associated with each Ranger_ID?
    # (note this differs from calculations above that examine Polygon_IDs per Ranger_ID)

        leech %>% group_by(dataset, Ranger_ID) %>%
            summarize(number.Lab_IDs = length(unique(Lab_ID)), .groups = "drop") %>%
            select(dataset, number.Lab_IDs) %>% table()

# polygons and samples

    # number of polygons

        leech %>% group_by(dataset) %>%
            summarize(number.PolygonIDs = length(unique(Polygon_ID)), .groups = "drop") # LSU: 126 polygons, SSU: 127 polygons

    # how many polygons were actually real polygons?
    # imputed Polygon_IDs are those derived from Ranger_IDs
    # real Polygon_IDs are actual polygons on the ground
    # LSU: 35 imputed, 91 real; SSU: 34 imputed, 93 real

        leech %>%
            mutate(polygon.type = ifelse(substr(Polygon_ID,1,1) == "R", "imputed", "real")) %>%
            distinct(dataset, Polygon_ID, polygon.type) %>%
            group_by(dataset, polygon.type) %>% tally()

    # how many Lab_IDs were associated with each Polygon_ID?

        leech %>% group_by(dataset, Polygon_ID) %>%
            summarize(number.Lab_IDs = length(unique(Lab_ID)), .groups = "drop") %>%
            select(dataset, number.Lab_IDs) %>% table()

    # mean/median number of Lab_IDs associated with each Polygon_ID

        LabIDs.per.PolygonID <- leech %>% group_by(dataset, Polygon_ID) %>%
            summarize(number.Lab_IDs = length(unique(Lab_ID)), .groups = "drop_last")

        LabIDs.per.PolygonID %>%
            summarize(mean.number.Lab_IDs = mean(number.Lab_IDs), # LSU: 5.182540, SSU: 5.826772
                median.number.Lab_IDs = median(number.Lab_IDs), .groups = "drop") # LSU: 3, SSU: 4

        LabIDs.per.PolygonID %>%
            ggplot(aes(x=number.Lab_IDs)) + stat_count() + facet_wrap("dataset") +
                labs(x="LabIDs per PolygonID", y="number of PolygonIDs")

        rm(LabIDs.per.PolygonID)

# species

    # how many species did we observe in total
    # LSU: 59 OTUs excluding humans, SSU: 72 OTUs excluding humans

        leech %>%
            group_by(dataset, OTU) %>% filter(sum(reads) > 0) %>%
            group_by(dataset) %>% distinct(dataset, OTU) %>% tally()

    # how many OTUs were identified to species level?
    # 58 identified to species level

        leech %>%
            select(consensus.short, consensus.species) %>%
            distinct() %>%
            mutate(identified.to.species = !grepl("\\d", consensus.species)) %>%
            group_by(identified.to.species) %>% tally()

    # how many species were detected in each Lab_ID?

        richness.per.LabID <- leech %>% group_by(dataset, Lab_ID) %>%
            summarize(richness = sum(reads>0), .groups = "drop")

        richness.per.LabID %>% ggplot(aes(x=richness)) + stat_count() + facet_wrap("dataset") +
                labs(x="OTUs per LabID", y="number of LabIDs")

        richness.per.LabID %>% group_by(dataset) %>%
            summarize(mean.richness = mean(richness), # LSU: 1.448698, SSU: 1.956757
                median.richness = median(richness), .groups = "drop") # LSU: 1, SSU: 2

        rm(richness.per.LabID)

    # how many species were detected in each Polygon_ID?

        richness.per.PolygonID <- leech %>%
            group_by(dataset, Polygon_ID, OTU) %>%
            summarize(detected = (sum(reads)>0), .groups = "drop_last") %>%
            summarize(richness = sum(detected), .groups = "drop")

        richness.per.PolygonID %>% ggplot(aes(x=richness)) + stat_count() + facet_wrap("dataset") +
                labs(x="OTUs per PolygonID", y="number of PolygonIDs")

        richness.per.PolygonID %>% group_by(dataset) %>%
            summarize(mean.richness = mean(richness), # LSU: 3.349206, SSU: 5.464567
                median.richness = median(richness), .groups = "drop") # LSU: 3, SSU: 4

        rm(richness.per.PolygonID)

    # how many Lab_IDs or Polygon_IDs was each species detected in?

        # median number of Lab_IDs per OTU
        leech %>%
            group_by(dataset, OTU) %>% summarize(num.LabIDs = sum(reads>0)) %>%
            group_by(dataset) %>% summarize(median.num.LabIDs = median(num.LabIDs))

        # median number of Polygon_IDs per OTU
        leech %>%
            group_by(dataset, Polygon_ID, OTU) %>% summarize(detected = (sum(reads)>0)) %>%
            group_by(dataset, OTU) %>% summarize(num.PolygonIDs = sum(detected)) %>%
            group_by(dataset) %>% summarize(median.num.PolygonIDs = median(num.PolygonIDs))

# observed species and Polygon_ID characteristics

    temp <- leech %>% filter(!is.na(latitude)) %>%
        group_by(dataset, Polygon_ID, shape_area_ha, shape_perimeter, elevation_mean,
            distance_to_nature_reserve_boundary, tpi_mean, distance_to_stream_mean, distance_to_road_median, OTU) %>%
        summarize(detected = (sum(reads) > 0), .groups = "drop_last") %>%
        summarize(richness = sum(detected), .groups = "drop")

    # size of polygons
        # perimeter ... no effect
            temp %>% ggplot(aes(x = shape_perimeter, y = richness)) + geom_point() +
                geom_smooth(aes(x = shape_perimeter, y=richness), method = "loess") + facet_wrap("dataset")
        # area ... nope
            temp %>% ggplot(aes(x = shape_perimeter, y = richness)) + geom_point() +
                geom_smooth(aes(x = shape_perimeter, y = richness), method = "loess") + facet_wrap("dataset")

    # elevation and distance to park edge are positively related to observed species richness
        # greater species richness at higher elevations
            temp %>% ggplot(aes(x = elevation_mean, y = richness)) + geom_point() +
                geom_smooth(aes(x = elevation_mean, y = richness), method = "loess") + facet_wrap("dataset")
        # positive correlation with distance to nature reserve boundary
            temp %>% ggplot(aes(x = distance_to_nature_reserve_boundary, y = richness)) + geom_point() +
                geom_smooth(aes(x = distance_to_nature_reserve_boundary, y = richness), method = "loess") + facet_wrap("dataset")
        # distance_to_nature_reserve_boundary and elevation_mean are positively correlated
            leech %>% filter(!is.na(latitude)) %>%
                distinct(Lab_ID, Polygon_ID, distance_to_nature_reserve_boundary, elevation_mean) %>%
                ggplot(aes(x = distance_to_nature_reserve_boundary, y = elevation_mean)) + geom_point()

    # not much to see in relation to TPI
        temp %>% ggplot(aes(x = tpi_mean, y = richness)) + geom_point() +
                geom_smooth(aes(x=tpi_mean, y=richness), method = "loess") + facet_wrap("dataset")
    # not much to see with distance to stream
        temp %>% ggplot(aes(x = distance_to_stream_mean, y = richness)) + geom_point() +
                geom_smooth(aes(x=distance_to_stream_mean, y=richness), method = "loess") + facet_wrap("dataset")
    # mild positive correlation with distance to road
        temp %>% ggplot(aes(x = distance_to_road_median, y = richness)) + geom_point() +
                geom_smooth(aes(x=distance_to_road_median, y=richness), method = "loess") + facet_wrap("dataset")

    rm(temp)

# detections vs leech quantity

    # plot

        leech %>% group_by(dataset, Lab_ID, leech_qty) %>% summarize(richness = sum(reads>0)) %>%
            ggplot(aes(x=leech_qty, y= richness)) + geom_jitter(height=0.3, width=0) + facet_wrap(~dataset) + labs(y="observed richness")

    # segmented regressions to illustrate drop off in slope

        # prepare data
        lsu.data <- leech %>% filter(dataset =="LSU") %>%
            group_by(Lab_ID, leech_qty) %>% summarize(richness = sum(reads>0))
        ssu.data <- leech %>% filter(dataset =="SSU") %>%
            group_by(Lab_ID, leech_qty) %>% summarize(richness = sum(reads>0))

        # regular linear models
        lsu.lm <- lm(richness ~ leech_qty, data= lsu.data)
        ssu.lm <- lm(richness ~ leech_qty, data= ssu.data)

        # segmented regression models
        lsu.seg <- segmented(lsu.lm, seg.Z = ~leech_qty, psi=50)
        ssu.seg <- segmented(ssu.lm, seg.Z = ~leech_qty, psi=50)

        # plots with regressions overlaid

        plot(richness ~ leech_qty, data = lsu.data)
        plot(lsu.seg, add=T, conf.level = 0.95, col= "red")

        plot(richness ~ leech_qty, data = ssu.data)
        plot(ssu.seg, add=T, conf.level = 0.95, col= "red")

        # test whether slope decreases
        davies.test(lsu.lm, seg.Z = ~leech_qty, k=10, alternative = "less")
        davies.test(ssu.lm, seg.Z = ~leech_qty, k=10, alternative = "less")
