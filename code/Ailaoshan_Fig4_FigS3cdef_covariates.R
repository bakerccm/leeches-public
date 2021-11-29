# R code to draw Fig 4

library("here")
library("tidyverse")
library("scales")
library("cowplot") # to layout individual plots in a grid
# see https://wilkelab.org/cowplot/articles/aligning_plots.html for options

########################################################################################
# get predictions generated previously

    load(here("rdata", "Ailaoshan_occupancy_covariates_LSU_community_predictions.rdata"))
    load(here("rdata", "Ailaoshan_occupancy_covariates_LSU_individual_predictions.rdata"))
    load(here("rdata", "Ailaoshan_occupancy_covariates_SSU_community_predictions.rdata"))
    load(here("rdata", "Ailaoshan_occupancy_covariates_SSU_individual_predictions.rdata"))

# relabel (non-avian) reptiles as squamates
    LSU.community.predictions <- LSU.community.predictions %>% lapply(function (X) {
            X %>% mutate(group = ifelse(as.character(group) == "amphibians/reptiles", "amphibians/squamates", as.character(group))) %>%
                mutate(group = factor(group, levels = c("mammals/birds","amphibians/squamates")))
    })
    LSU.individual.predictions <- LSU.individual.predictions %>%
        lapply(function (X) {
            X %>% mutate(group.name  = ifelse(group.name == "amphibians/reptiles", "amphibians/squamates", group.name)) %>%
                mutate(consensus.class  = ifelse(consensus.class == "reptiles", "squamates", consensus.class))
    })
    SSU.community.predictions <- SSU.community.predictions %>% lapply(function (X) {
        X %>% mutate(group = ifelse(as.character(group) == "amphibians/reptiles", "amphibians/squamates", as.character(group))) %>%
            mutate(group = factor(group, levels = c("mammals/birds","amphibians/squamates")))
    })
    SSU.individual.predictions <- SSU.individual.predictions %>%
        lapply(function (X) {
            X %>% mutate(group.name  = ifelse(group.name == "amphibians/reptiles", "amphibians/squamates", group.name)) %>%
                mutate(consensus.class  = ifelse(consensus.class == "reptiles", "squamates", consensus.class))
    })

########################################################################################
# generate individual plots

plot.list <- list()

    # set colors for taxonomic groups in community plots to avoid confusion with taxonomic classes
        group.colors <- c("mammals/birds" = "#00BFC4", "amphibians/squamates" = "#F8766D")

        # ggplot colors
        # hue_pal()(2) # "#F8766D" "#00BFC4"
        # show_col(hue_pal()(2))

# LSU community

    plot.list$LSU.community.elevation <- LSU.community.predictions$occupancy.elevation %>%
        ggplot(aes(x=`elevation (m)`)) +
            geom_ribbon(aes(ymin = `2.5%`, ymax = `97.5%`, fill = group), alpha = 0.4) +
            geom_line(aes(y = mean, col = group)) + labs(y = "mean occupancy") +
            theme(legend.title = element_blank()) + scale_fill_manual(values = group.colors) + scale_color_manual(values = group.colors)

    plot.list$LSU.community.reserve <- LSU.community.predictions$occupancy.reserve %>%
        ggplot(aes(x=`distance to reserve edge (m)`)) +
            geom_ribbon(aes(ymin = `2.5%`, ymax = `97.5%`, fill = group), alpha = 0.4) +
            geom_line(aes(y = mean, col = group)) + labs(y = "mean occupancy") +
            theme(legend.title = element_blank()) + scale_fill_manual(values = group.colors) + scale_color_manual(values = group.colors)

    plot.list$LSU.community.detection <- LSU.community.predictions$detection %>%
        ggplot(aes(x=`number of leeches`)) +
            geom_ribbon(aes(ymin = `2.5%`, ymax = `97.5%`, fill = group), alpha = 0.4) +
            geom_line(aes(y = mean, col = group)) + labs(y = "mean detection") +
            theme(legend.title = element_blank()) + scale_fill_manual(values = group.colors) + scale_color_manual(values = group.colors)

# SSU community

    plot.list$SSU.community.elevation <- SSU.community.predictions$occupancy.elevation %>%
        ggplot(aes(x=`elevation (m)`)) +
            geom_ribbon(aes(ymin = `2.5%`, ymax = `97.5%`, fill = group), alpha = 0.4) +
            geom_line(aes(y = mean, col = group)) + labs(y = "mean occupancy") +
            theme(legend.title = element_blank()) + scale_fill_manual(values = group.colors) + scale_color_manual(values = group.colors)

    plot.list$SSU.community.detection <- SSU.community.predictions$detection %>%
        ggplot(aes(x=`number of leeches`)) +
            geom_ribbon(aes(ymin = `2.5%`, ymax = `97.5%`, fill = group), alpha = 0.4) +
            geom_line(aes(y = mean, col = group)) + labs(y = "mean detection") +
            theme(legend.title = element_blank()) + scale_fill_manual(values = group.colors) + scale_color_manual(values = group.colors)

# LSU individual

    plot.list$LSU.individual.elevation <- LSU.individual.predictions$elev %>%
        ggplot(aes(x=elev.unscaled, y=`estimated occupancy`, group = OTU, color = consensus.class)) +
        geom_line() + labs(x = "elevation (m)", y = "occupancy", color = "taxonomic group") + theme(legend.title = element_blank())

    plot.list$LSU.individual.reserve <- LSU.individual.predictions$reserve %>%
        ggplot(aes(x=reserve.unscaled, y=`estimated occupancy`, group = OTU, color = consensus.class)) +
        geom_line() + labs(x = "distance to reserve edge (m)", y = "occupancy", color = "taxonomic group") + theme(legend.title = element_blank())

    plot.list$LSU.individual.detection <- LSU.individual.predictions$numleeches %>%
        ggplot(aes(x=leeches, y=`estimated detection`, group = OTU, color = consensus.class)) +
        geom_line() + labs(x = "number of leeches", y = "detection", color = "taxonomic group") + theme(legend.title = element_blank()) +
        scale_y_continuous(breaks = seq(0, 0.9, 0.3))

# SSU individual

    plot.list$SSU.individual.elevation <- SSU.individual.predictions$elev %>%
        ggplot(aes(x=elev.unscaled, y=`estimated occupancy`, group = OTU, color = consensus.class)) +
        geom_line() + labs(x = "elevation (m)", y = "occupancy", color = "taxonomic group") + theme(legend.title = element_blank())

    plot.list$SSU.individual.detection <- SSU.individual.predictions$numleeches %>%
        ggplot(aes(x=leeches, y=`estimated detection`, group = OTU, color = consensus.class)) +
        geom_line() + labs(x = "number of leeches", y = "detection", color = "taxonomic group") + theme(legend.title = element_blank())

########################################################################################
# arrange occupancy plots in grid for Fig 4

    # see https://wilkelab.org/cowplot/articles/shared_legends.html for info on shared legends

    # plots without legends
        occupancy.plot.grid <- plot_grid(
            plot.list$LSU.community.elevation + theme(legend.position="none"),
            plot.list$LSU.individual.elevation + theme(legend.position="none"),
            plot.list$LSU.community.reserve + theme(legend.position="none"),
            plot.list$LSU.individual.reserve + theme(legend.position="none"),
            plot.list$SSU.community.elevation + theme(legend.position="none"),
            plot.list$SSU.individual.elevation + theme(legend.position="none"),
            ncol=2, align = 'hv',
            hjust=-1, vjust = -1)

    # extract legends
        occupancy.community.legend <- get_legend(
          # create some space to the left of the legend
          plot.list$LSU.community.elevation + guides(color = guide_legend(nrow = 1)) + theme(legend.position = "bottom")
        )
        occupancy.individual.legend <- get_legend(
          # create some space to the left of the legend
          plot.list$LSU.individual.elevation + guides(color = guide_legend(nrow = 1)) + theme(legend.position = "bottom")
        )

    # arrange legends side by side
        occupancy.legend.grid <- plot_grid(occupancy.community.legend, occupancy.individual.legend, ncol = 2)

    # add the legend below the plots
    # give it 20% of the height of one plot (via rel_heights)
        plot_grid(occupancy.plot.grid, occupancy.legend.grid, ncol = 1, rel_heights = c(3, .2))

    # save to file
        ggsave(filename = here("figures", "Fig4_occupancy_covariates.pdf"), width = 8, height = 9, useDingbats = FALSE)

########################################################################################
# arrange occupancy plots in grid for Fig S3cdef

    # see https://wilkelab.org/cowplot/articles/shared_legends.html for info on shared legends

    # plots without legends
        detection.plot.grid <- plot_grid(
            plot.list$LSU.community.detection + theme(legend.position="none"),
            plot.list$LSU.individual.detection + theme(legend.position="none"),
            plot.list$SSU.community.detection + theme(legend.position="none"),
            plot.list$SSU.individual.detection + theme(legend.position="none"),
            ncol=2, align = 'hv',
            hjust=-1, vjust = -1)

    # extract legends
        detection.community.legend <- get_legend(
            # create some space to the left of the legend
            plot.list$LSU.community.detection + guides(color = guide_legend(nrow = 1)) + theme(legend.position = "bottom")
        )
        detection.individual.legend <- get_legend(
            # create some space to the left of the legend
            plot.list$LSU.individual.detection + guides(color = guide_legend(nrow = 1)) + theme(legend.position = "bottom")
        )

    # arrange legends side by side
        detection.legend.grid <- plot_grid(detection.community.legend, detection.individual.legend, ncol = 2)

    # add the legend below the plots
    # give it 20% of the height of one plot (via rel_heights)
        plot_grid(detection.plot.grid, detection.legend.grid, ncol = 1, rel_heights = c(2, .2))

    # save to file
    ggsave(filename = here("figures", "FigS3cdef_detection_covariates.pdf"), width = 8, height = 6.2, useDingbats = FALSE)
