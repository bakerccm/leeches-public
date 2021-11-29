# R code to draw Fig 2c

library("here")
library("tidyverse")
library("ggrepel")

# get data

    load(file=here("rdata","Ailaoshan_OTU_table.rdata"))

    # model summaries generated with Ailaoshan_model_summary.R
    model.summary <- list(
        LSU = readRDS(file = here("rds", "Ailaoshan_model_summary_200_final_LSU.rds")),
        SSU = readRDS(file = here("rds", "Ailaoshan_model_summary_200_final_SSU.rds"))
    )
    estocc.output <- lapply(model.summary, FUN = function (X) as.data.frame(X$estocc.output))
    beta0.output <- lapply(model.summary, FUN = function (X) as.data.frame(X$beta0.output))
    gamma0.output <- lapply(model.summary, FUN = function (X) as.data.frame(X$gamma0.output))
    rm(model.summary)

# prepare occupancy and detection estimates

    occupancy.estimates <- bind_rows(LSU = estocc.output$LSU, SSU = estocc.output$SSU, .id = "dataset") %>%
        select(dataset, OTU, mean, `2.5%`, `97.5%`) %>%
        rename(occupancy_mean = mean, `occupancy_2.5%` = `2.5%`, `occupancy_97.5%` = `97.5%`)

    # N.B. Using estocc.output (i.e. actual fraction of sites occupied, z_k / n) gives similar but slightly different
    # results than using beta0.output to measure occupancy probability.
    # occupancy.estimates <- bind_rows(LSU = beta0.output$LSU, SSU = beta0.output$SSU, .id = "dataset") %>%
    #     select(dataset, OTU, mean, `2.5%`, `97.5%`) %>%
    #     rename(occupancy_mean = mean, `occupancy_2.5%` = `2.5%`, `occupancy_97.5%` = `97.5%`) %>%
    #     mutate(occupancy_mean = plogis(occupancy_mean), `occupancy_2.5%` = plogis(`occupancy_2.5%`), `occupancy_97.5%` = plogis(`occupancy_97.5%`))

    detection.estimates <- bind_rows(LSU = gamma0.output$LSU, SSU = gamma0.output$SSU, .id = "dataset") %>%
        select(dataset, OTU, prob_mean, `prob_2.5%`, `prob_97.5%`) %>%
        rename(detection_prob_mean = prob_mean, `detection_prob_2.5%` = `prob_2.5%`, `detection_prob_97.5%` = `prob_97.5%`)

    model.estimates <- full_join(occupancy.estimates, detection.estimates, by = c("dataset", "OTU"))

# prepare data for graphing

    # join observed occupancy with model estimates
    OTU.occupancy.summary <- leech %>%
        # calculate observed presence/absence (i.e. detection = true/false) by Polygon_ID
            group_by(dataset, OTU, consensus.short, consensus.class, consensus.order, consensus.family, consensus.genus, consensus.species, Polygon_ID) %>%
            summarize(occupied = ifelse(sum(reads) > 0 , 1, 0), .groups = "drop_last") %>%
        # calculate fraction of sites occupied
            summarize(occupancy_observed = sum(occupied) / n(), .groups = "drop") %>%
        # attach model estimates for occupancy and detection
            left_join(model.estimates, by = c("dataset", "OTU")) %>%
        # add rank for observed occupancy, estimated occupancy and estimated detection probability within each dataset
            group_by(dataset) %>%
            mutate(occupancy_observed_rank = rank(desc(occupancy_observed))) %>% # observed occupancy
            mutate(occupancy_mean_rank = rank(desc(occupancy_mean))) %>% # estimated occupancy probability
            mutate(detection_prob_mean_rank = rank(desc(detection_prob_mean))) %>% # estimated detection probability
            ungroup() %>%
        # relabel (non-avian) reptiles as squamates
            mutate(consensus.class = ifelse(consensus.class == "Reptiles", "Squamates", consensus.class))
        
# draw plots

    OTU.occupancy.summary.plot <- OTU.occupancy.summary %>%
        mutate(consensus.class = tolower(consensus.class)) %>%
        ggplot(aes(x=occupancy_mean, y=detection_prob_mean, col=consensus.class, label=consensus.short)) + geom_point(alpha=0.5) +
            labs(x="fraction of patrol areas occupied", y = "detection probability per 100 leeches") +
            theme(legend.title = element_blank()) +
            facet_wrap("dataset") + xlim(0,1)

    # unlabelled
    OTU.occupancy.summary.plot
    ggsave(here("figures","Fig2c_occupancy_detection.pdf"), width=7, height=4)

    # with text labels
    # use this to work out how to label points on plot above in illustrator
    OTU.occupancy.summary.plot + geom_text(size = 1, show.legend = FALSE)
    ggsave(here("figures","Fig2c_occupancy_detection_labels.pdf"), width=7, height=4)

    # with error bars
    OTU.occupancy.summary.plot +
        geom_errorbar(aes(ymin = `detection_prob_2.5%`, ymax = `detection_prob_97.5%`)) +
        geom_errorbarh(aes(xmin = `occupancy_2.5%`, xmax = `occupancy_97.5%`))

    rm(OTU.occupancy.summary.plot)

# median species-wise occupancy and detection estimates

    OTU.occupancy.summary %>%
        group_by(dataset) %>%
        summarize(median.occupancy_mean = median(occupancy_mean) %>% round(.,2), median.detection_prob_mean = median(detection_prob_mean) %>% round(.,2), .groups = "drop")
