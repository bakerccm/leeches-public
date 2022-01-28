# R code for Fig 3 (occupancy maps and scatterplots) and Figs S3a,b (occupancy histograms)

library("here")
library("sf")
library("raster")
library("sp")
library("tidyverse")
library("ggExtra") # for marginal histograms
library("cowplot")

# get data

    load(file=here("rdata","Ailaoshan_OTU_table.rdata"))

    load(file=here("rdata","Ailaoshan_environmental.rdata"))

    load(file=here("rdata","Ailaoshan_gis.rdata"))

    # model summaries generated with Ailaoshan_model_summary.R
    model.summary <- list(
        LSU = readRDS(file = here("rds", "Ailaoshan_model_summary_200_final_LSU.rds")),
        SSU = readRDS(file = here("rds", "Ailaoshan_model_summary_200_final_SSU.rds"))
    )
    Nsite.output <- lapply(model.summary, FUN = function (X) as.data.frame(X$Nsite.output))
    rm(model.summary)

# colorblind-friendly plot colors from Wong 2011 (https://www.nature.com/articles/nmeth.1618)

    plot.colors <- c("estimated" = "black", "observed" = "grey")

# add observed richness data to ailaoshan.polygons sf object

    observed.OTU.richness <- leech %>%
        group_by(dataset, Polygon_ID, OTU) %>%
        summarize(OTU.observed = ifelse(sum(reads) > 0, 1, 0), .groups = "drop_last") %>%
        summarize(richness = sum(OTU.observed), .groups = "drop") %>%
        pivot_wider(id_cols = Polygon_ID, values_from = richness, names_prefix = "observed.richness.", names_from = dataset)

    ailaoshan.polygons <- ailaoshan.polygons %>%
        mutate(Polygon_ID = as.character(Polygon_ID)) %>%
        left_join(observed.OTU.richness, by="Polygon_ID")

# add estimated richness data to ailaoshan.polygons sf object

    estimated.OTU.richness <- bind_rows(LSU = Nsite.output$LSU, SSU = Nsite.output$SSU, .id = "dataset") %>%
        select(dataset, Polygon_ID, mean) %>%
        pivot_wider(id_cols = Polygon_ID, values_from = mean, names_prefix = "estimated.richness.", names_from = dataset)

    ailaoshan.polygons <- ailaoshan.polygons %>%
        left_join(estimated.OTU.richness, by="Polygon_ID")

# prepare data for plotting

    OTU.richness <- inner_join(observed.OTU.richness, estimated.OTU.richness, by = "Polygon_ID") %>%
        pivot_longer(-Polygon_ID, names_to = c("source", NA, "dataset"), names_sep = "\\.", values_to = "richness") %>%
        filter(complete.cases(.))

    OTU.richness.medians <- OTU.richness %>% group_by(`source`, dataset) %>%
        summarize(median.richness = median(richness), .groups = "drop")

# how many observed species per replicate (i.e. LabID)?

    leech %>% filter(!is.na(leech_qty)) %>%
        group_by(dataset, Lab_ID) %>%
        summarize(observed.richness = sum(ifelse(reads > 0, 1, 0)), .groups = "drop_last") %>%
        summarize(median.observed.richness = median(observed.richness), .groups = "drop")

# how many replicates per polygonID?

    leech %>% filter(!is.na(leech_qty)) %>%
        group_by(dataset, Polygon_ID, OTU) %>%
        summarize(num.replicates = n(), .groups = "drop_last") %>%
        summarize(num.replicates = mean(num.replicates), .groups = "drop_last") %>%
        summarize(num.replicates = round(median(num.replicates), 1), .groups = "drop")

# histograms for Figs S3a,b

    LSU.medians <- OTU.richness.medians %>% filter(dataset == "LSU")
    SSU.medians <- OTU.richness.medians %>% filter(dataset == "SSU")

    # plot Fig S3a: LSU richness
        OTU.richness %>% filter(dataset == "LSU") %>%
            ggplot(aes(x = richness, color = source)) +
                geom_histogram(fill = "white", alpha = 0.5, position = "identity", binwidth = 2) + coord_cartesian(xlim=c(0,46), ylim=c(0,55)) +
                geom_vline(data = LSU.medians, aes(xintercept = median.richness, color = source), linetype = "dashed") +
                labs(x = "LSU species richness") + scale_color_manual(values = plot.colors) +
                theme(legend.title = element_blank(), legend.position = c(0.842, 0.71), legend.margin=margin(t = 0, r = 0.2, b = 0.2, l = 0.2, unit = "cm")) +
                geom_label(data = LSU.medians, aes(x = median.richness + 0.1, y = 54,
                    label = paste0("median = ", round(median.richness, 0))), show.legend = FALSE, size = 3)
        ggsave(here("figures","FigS3a_LSU_richhist.pdf"), width=4, height=3, useDingbats = FALSE)

    # plot Fig S3b: SSU richness
        OTU.richness %>% filter(dataset == "SSU") %>%
            ggplot(aes(x = richness, color = source)) +
                geom_histogram(fill = "white", alpha = 0.5, position = "identity", binwidth = 2) + coord_cartesian(xlim=c(0,46), ylim=c(0,55)) +
                geom_vline(data = SSU.medians, aes(xintercept = median.richness, color = source), linetype = "dashed") +
                labs(x = "SSU species richness") + scale_color_manual(values = plot.colors) +
                theme(legend.title = element_blank(), legend.position = c(0.842, 0.71), legend.margin=margin(t = 0, r = 0.2, b = 0.2, l = 0.2, unit = "cm")) +
                geom_label(data = SSU.medians, aes(x = median.richness + 0.1, y = 54,
                    label = paste0("median = ", round(median.richness, 0))), show.legend = FALSE, size = 3)
        ggsave(here("figures","FigS3b_SSU_richhist.pdf"), width=4, height=3, useDingbats = FALSE)

    rm(LSU.medians, SSU.medians)

# draw estimated richness maps for Figs 3a,b

    # Fig 3a: LSU observed richness
        ailaoshan.polygons %>%
            ggplot() + geom_sf(aes(fill=observed.richness.LSU), size = 0.1, color = "black") +
                theme(legend.position=c(0,0), legend.justification=c(0,0), legend.background = element_rect(fill=NA)) +
                coord_sf() + labs(fill = "LSU\nobserved\nspecies\nrichness") + scale_fill_continuous(na.value=NA)
        ggsave(here("figures","Fig3a_LSU_observed_richness.pdf"), width=4, height=4)

    # Fig 3b: SSU observed richness
        ailaoshan.polygons %>%
            ggplot() + geom_sf(aes(fill=observed.richness.SSU), size = 0.1, color = "black") +
                theme(legend.position=c(0,0), legend.justification=c(0,0), legend.background = element_rect(fill=NA)) +
                coord_sf() + labs(fill = "SSU\nobserved\nspecies\nrichness") + scale_fill_continuous(na.value=NA)
        ggsave(here("figures","Fig3b_SSU_observed_richness.pdf"), width=4, height=4)

# plot estimated richness maps for Figs 3c,d

    # Fig 3c: LSU estimated richness
        ailaoshan.polygons %>%
            ggplot() + geom_sf(aes(fill=estimated.richness.LSU), size = 0.1, color = "black") +
                theme(legend.position=c(0,0), legend.justification=c(0,0), legend.background = element_rect(fill=NA)) +
                coord_sf() + labs(fill = "LSU\nestimated\nspecies\nrichness") + scale_fill_continuous(na.value=NA)
        ggsave(here("figures","Fig3c_LSU_estimated_richness.pdf"), width=4, height=4)

    # Fig 3d: SSU estimated richness
        ailaoshan.polygons %>%
            ggplot() + geom_sf(aes(fill=estimated.richness.SSU), size = 0.1, color = "black") +
                theme(legend.position=c(0,0), legend.justification=c(0,0), legend.background = element_rect(fill=NA)) +
                coord_sf() + labs(fill = "SSU\nestimated\nspecies\nrichness") + scale_fill_continuous(na.value=NA)
        ggsave(here("figures","Fig3d_SSU_estimated_richness.pdf"), width=4, height=4)

# estimated species richness (Nsite) vs env covariates for Figs 3e,f

    # LSU scatter plot

        LSU.scatter <- Nsite.output$LSU %>% left_join(env.data, by = "Polygon_ID") %>%
            select(mean, elevation_median, distance_to_nature_reserve_boundary) %>%
            pivot_longer(cols =c(elevation_median, distance_to_nature_reserve_boundary), values_to = "covariate_value", names_to = "covariate") %>%
            mutate(covariate = factor(covariate, levels = c("elevation_median","distance_to_nature_reserve_boundary"))) %>%
            rename(`estimated species richness` = mean) %>%
            ggplot(aes(x = covariate_value, y = `estimated species richness`)) + geom_point() +labs(x="")+
            facet_wrap(~covariate, scales = "free_x") + theme(strip.background = element_blank(), strip.text.x = element_blank()) +
            ylim(20,48)

        LSU.scatter
        ggsave(here("figures", "Fig3e_LSU_scatter.pdf"), width=4, height=3, useDingbats = FALSE)

    # LSU histogram

        LSU.hist <- Nsite.output$LSU %>% left_join(env.data, by = "Polygon_ID") %>%
            select(mean, elevation_median) %>%
            rename(`estimated species richness` = mean) %>%
            ggplot(aes(x = elevation_median, y = `estimated species richness`)) + geom_point() +
            ylim(20,48)

        # need to wrap the plot command with ggsave (rather than running sequenctially) for it to work when running from command line
        ggsave(
            here("figures", "Fig3e_LSU_histogram.pdf"),
            LSU.hist %>% ggMarginal(type = "histogram", margins = "y", binwidth = 1, fill = "cornsilk3"),
            width=4, height=3, useDingbats = FALSE
        )

    # SSU scatter plot

        SSU.scatter <- Nsite.output$SSU %>% left_join(env.data, by = "Polygon_ID") %>%
            select(mean, elevation_median) %>%
            rename(`estimated species richness` = mean) %>%
            ggplot(aes(x = elevation_median, y = `estimated species richness`)) + geom_point() +
            labs(x = "elevation (m)") + ylim(20,48)

        SSU.scatter
        ggsave(here("figures", "Fig3f_SSU_scatter.pdf"), width=2, height=3, useDingbats = FALSE)

    # SSU histogram

        # need to wrap the plot command with ggsave (rather than running sequenctially) for it to work when running from command line
        ggsave(
            here("figures", "Fig3f_SSU_histogram.pdf"),
            SSU.scatter  %>% ggMarginal(type = "histogram", margins = "y", binwidth = 1, fill = "cornsilk3"),
            width=4, height=3, useDingbats = FALSE
        )
