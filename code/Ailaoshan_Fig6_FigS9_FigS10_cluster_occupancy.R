# R code to produce tables of occupancy for each species in each cluster

library("here")
library("vegan")
library("tidyverse")
library("cowplot") # to layout individual plots in a grid
# see https://wilkelab.org/cowplot/articles/aligning_plots.html for options

########################################################################################
# get data

    load(file=here("rdata","Ailaoshan_OTU_table.rdata"))

    # occupancy for each species in each cluster
    cluster.occ <- list(
        LSU = readRDS(here("rds", "Ailaoshan_final_LSU_cluster_occupancy.rds")),
        SSU = readRDS(here("rds", "Ailaoshan_final_SSU_cluster_occupancy.rds"))
    )

    # model summaries generated with Ailaoshan_model_summary.R
    model.summary <- list(
        LSU = readRDS(file = here("rds", "Ailaoshan_model_summary_200_final_LSU.rds")),
        SSU = readRDS(file = here("rds", "Ailaoshan_model_summary_200_final_SSU.rds"))
    )
    estocc.output <- lapply(model.summary, FUN = function (X) as.data.frame(X$estocc.output))
    gamma0.output <- lapply(model.summary, FUN = function (X) as.data.frame(X$gamma0.output))
    rm(model.summary)

########################################################################################
# colorblind-friendly plot colors from Wong 2011 (https://www.nature.com/articles/nmeth.1618)

    plot.colors <- c("amphibians" = rgb(0, 158, 115, maxColorValue = 255), # bluish green
                     "birds" = rgb(86, 180, 233, maxColorValue = 255), # sky blue
                     "mammals" = rgb(230, 159, 0, maxColorValue = 255), # orange
                     "squamates" = rgb(0, 0, 0, maxColorValue = 255)) # black

########################################################################################
# add taxonomic info to cluster occupancy results

    cluster.occ <- lapply(cluster.occ, function (X) {
        X %>% left_join(leech %>% select(OTU, consensus, consensus.short, consensus.class, consensus.order, consensus.family, consensus.genus, consensus.species, Chinese_common_name, domestic, AdultBodyMass_g, AdultBodyMass_g_source) %>% distinct(), by = "OTU") %>%
            # relabel (non-avian) reptiles as squamates
                mutate(consensus.class = ifelse(consensus.class == "Reptiles", "Squamates", consensus.class)) %>%
            # remove "augmented" species -- theseare NA for consensus.class
                filter(!is.na(consensus.class))
    })

########################################################################################
# top species in each cluster

    cluster.occ$LSU %>% filter(cluster == 1) %>% arrange(desc(mean))

    cluster.occ$LSU %>% filter(cluster == 2) %>% arrange(desc(mean))

    cluster.occ$LSU %>% filter(cluster == 3) %>% arrange(desc(mean))

########################################################################################
# graphs of all species ordered by occupancy in each cluster

cluster.occ$LSU %>% filter(cluster == 1) %>% arrange(desc(mean)) %>%
    ggplot(aes(x=reorder(consensus.short,mean), y=mean,ymin = `2.5%`, ymax=`97.5%`, color = consensus.class)) + geom_point() + geom_errorbar() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.45), axis.title.x = element_blank())

cluster.occ$LSU %>% filter(cluster == 2) %>% arrange(desc(mean)) %>%
    ggplot(aes(x=reorder(consensus.short,mean), y=mean,ymin = `2.5%`, ymax=`97.5%`, color = consensus.class)) + geom_point() + geom_errorbar() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.45), axis.title.x = element_blank())

cluster.occ$LSU %>% filter(cluster == 3) %>% arrange(desc(mean)) %>%
    ggplot(aes(x=reorder(consensus.short,mean), y=mean,ymin = `2.5%`, ymax=`97.5%`, color = consensus.class)) + geom_point() + geom_errorbar() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.45), axis.title.x = element_blank())

########################################################################################
# top species in cluster 1
cluster.occ$LSU %>%
    select(OTU, consensus.short, cluster, mean,`2.5%`,`97.5%`) %>%
    pivot_wider(id_cols = c(OTU, consensus.short), names_from = cluster, values_from = c(mean,`2.5%`,`97.5%`), names_glue = "cluster{cluster}_{.value}") %>%
    select(OTU, consensus.short, starts_with("cluster1_"), starts_with("cluster2_"), starts_with("cluster3_")) %>%
    arrange(desc(cluster1_mean))

# top species in cluster 3
cluster.occ$LSU %>%
    select(OTU, consensus.short, cluster, mean,`2.5%`,`97.5%`) %>%
    pivot_wider(id_cols = c(OTU, consensus.short), names_from = cluster, values_from = c(mean,`2.5%`,`97.5%`), names_glue = "cluster{cluster}_{.value}") %>%
    select(OTU, consensus.short, starts_with("cluster1_"), starts_with("cluster2_"), starts_with("cluster3_")) %>%
    arrange(desc(cluster3_mean))

########################################################################################
# mammals over 10kg
LSU.10kg.mammals <- cluster.occ$LSU %>%
    filter(consensus.class == "Mammals", AdultBodyMass_g > 10000)
SSU.10kg.mammals <- cluster.occ$SSU %>%
    filter(consensus.class == "Mammals", AdultBodyMass_g > 10000)

# facet_wrap
    LSU.10kg.mammals %>%
        ggplot(aes(x = cluster, col = domestic)) +
        geom_point(aes(y=mean)) +
        geom_errorbar(aes(ymin = `2.5%`, ymax = `97.5%`), width = 0.2) +
        facet_wrap(~ consensus.short)
    SSU.10kg.mammals %>%
        ggplot(aes(x = cluster, col = domestic)) +
        geom_point(aes(y=mean)) +
        geom_errorbar(aes(ymin = `2.5%`, ymax = `97.5%`), width = 0.2) +
        facet_wrap(~ consensus.short)

# by cluster
    LSU.10kg.mammals %>%
        ggplot(aes(x = consensus.short, col = domestic)) +
        geom_point(aes(y=mean)) + geom_errorbar(aes(ymin = `2.5%`, ymax = `97.5%`), width = 0.2) +
        facet_wrap(~ cluster, ncol = 1)
    SSU.10kg.mammals %>%
        ggplot(aes(x = consensus.short, col = domestic)) +
        geom_point(aes(y=mean)) + geom_errorbar(aes(ymin = `2.5%`, ymax = `97.5%`), width = 0.2) +
        facet_wrap(~ cluster, ncol = 1)

########################################################################################
# max cluster occupancy for all species
    cluster.occ$LSU %>%
        group_by(OTU) %>% mutate(max.occupancy = max(mean)) %>% ungroup() %>%
        select(consensus.short, max.occupancy) %>% distinct() %>%
        ggplot(aes(x = reorder(consensus.short, max.occupancy), y= max.occupancy)) + geom_point() +
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.45), axis.title.x = element_blank())

# cluster occupancy for species with max occupancy > x% -- arranged by cluster
    cluster.occ$LSU %>%
        group_by(OTU) %>% mutate(max.occupancy = max(mean)) %>% ungroup() %>%
        filter(max.occupancy >= 0.75) %>%
        ggplot(aes(x = consensus.short, col = domestic)) +
            geom_point(aes(y=mean)) + geom_errorbar(aes(ymin = `2.5%`, ymax = `97.5%`), width = 0.2) +
            facet_wrap(~ cluster, ncol = 1) +
            theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.45), axis.title.x = element_blank())

# cluster occupancy for mammals with overall occupancy > x% -- arranged by cluster
    cluster.occ$LSU %>%
        left_join(estocc.output$LSU %>% select(OTU, mean) %>% rename(estocc.mean = mean), by = "OTU") %>%
        filter(consensus.class == "Mammals", estocc.mean >= 0.30) %>%
        ggplot(aes(x = reorder(consensus.short, -estocc.mean), col = domestic)) +
            geom_point(aes(y=mean)) + geom_errorbar(aes(ymin = `2.5%`, ymax = `97.5%`), width = 0.2) +
            facet_wrap(~ cluster, ncol = 1) +
            theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.45), axis.title.x = element_blank())

    # cluster occupancy for species with overall occupancy > x% and detection >x% -- arranged by cluster -- trying to get more reliable estimates and more abundant species
    cluster.occ$LSU %>%
        left_join(estocc.output$LSU %>% select(OTU, mean) %>% rename(estocc.mean = mean), by = "OTU") %>%
        left_join(gamma0.output$LSU %>% select(OTU, prob_mean) %>% rename(gamma0.prob_mean = prob_mean), by = "OTU") %>%
        filter(estocc.mean >= 0.3, gamma0.prob_mean >= 0.1) %>%
        ggplot(aes(x = reorder(consensus.short, -mean), col = domestic)) +
        geom_point(aes(y=mean)) + geom_errorbar(aes(ymin = `2.5%`, ymax = `97.5%`), width = 0.2) +
        facet_wrap(~ cluster, ncol = 1) +
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.45), axis.title.x = element_blank())

########################################################################################
# cluster occupancy for species with max occupancy > x% and detection > x% -- arranged by cluster -- ordered by occupancy in cluster 3

    # prepare data
    clusterocc.plot.data <- sapply(names(cluster.occ), simplify = FALSE, function (X) {
        cluster.occ[[X]] %>%
            # get detection estimates
                left_join(estocc.output[[X]] %>% select(OTU, mean) %>% rename(estocc.mean = mean), by = "OTU") %>%
                left_join(gamma0.output[[X]] %>% select(OTU, prob_mean) %>% rename(gamma0.prob_mean = prob_mean), by = "OTU") %>%
            # calculate max occupancy out of the three clusters
                group_by(OTU) %>% mutate(max.occupancy = max(mean)) %>% ungroup() %>%
            # add dummy reptile rows to ensure reptiles facet is wide enough in supplementary figure
                bind_rows(tibble(consensus.short = paste0("dummy", 1:4), consensus.class = rep("squamates",4), mean = c(0,0,1,1), cluster = "3")) %>%
            # order OTUs by occupancy in cluster 3
                mutate(consensus.short = as.factor(consensus.short)) %>%
                group_by(OTU) %>% mutate(mean.cluster3 = mean[cluster==3]) %>% ungroup() %>%
                mutate(consensus.short = reorder(consensus.short, mean.cluster3)) %>%
            # make consensus.class lower case and relabel clusters
                mutate(consensus.class = tolower(consensus.class)) %>%
                mutate(cluster = recode(cluster, "1" = "high", "2" = "intermediate", "3" = "low")) %>%
                mutate(cluster = factor(cluster, levels = c("low", "intermediate", "high")))
    })

    ########################################################################################
    ## main text plots ##

    # individual plots
    clusterocc.plot.list <- lapply(clusterocc.plot.data, function (X) {
        X %>%
            # get only common taxa with reasonable detection probs
                filter(estocc.mean >= 0.4, gamma0.prob_mean >= 0.1) %>%
                #filter(max.occupancy >= 0.67, gamma0.prob_mean >= 0.1) %>%
            ggplot(aes(x = cluster, col = consensus.class)) +
                geom_point(aes(y = mean), size = 1.7) + geom_errorbar(aes(ymin = `2.5%`, ymax = `97.5%`), width = 0.25) +
                geom_linerange(aes(ymin = `25%`, ymax = `75%`), size = 1.2, show.legend = FALSE) +
                facet_wrap(~ reorder(consensus.short, mean.cluster3), nrow = 4) +
                theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.45), axis.title.x = element_blank()) +
                theme(legend.title = element_blank()) +
                labs(y = "estimated occupancy by elevation cluster") +
                scale_color_manual(breaks = c("amphibians", "birds", "mammals"), values = plot.colors)
    })

    # arrange plots in grid

    # see https://wilkelab.org/cowplot/articles/shared_legends.html for info on shared legends

    # plots without legends
    clusterocc.plot.grid <- plot_grid(
        clusterocc.plot.list$LSU + theme(legend.position="none"),
        clusterocc.plot.list$SSU + theme(legend.position="none"),
        ncol=2, align = 'hv', axis = "bt",
        hjust = -1, rel_widths = c(5.5,10))

    # extract legend from one of the plots
    clusterocc.plot.legend <- get_legend(
        clusterocc.plot.list$SSU + guides(color = guide_legend(nrow = 1)) + theme(legend.title = element_blank(), legend.position = "bottom")
    )

    # add the legend below the plots
    # give it 20% of the height of one plot (via rel_heights)
    plot_grid(clusterocc.plot.grid, clusterocc.plot.legend, ncol = 1, rel_heights = c(1, 0.05))

    # save plots to file
    ggsave(filename = here("figures", "Fig6ab_cluster_occupancy.pdf"), width = 10, height = 6, useDingbats = FALSE)

    ########################################################################################
    ## supplementary figures ##

    suppl.clusterocc.plot.list <- lapply(clusterocc.plot.data, function (X) {
        X %>% mutate(cluster = factor(cluster, levels = c("high", "intermediate", "low"))) %>%
            ggplot(aes(x = consensus.short, col = consensus.class)) +
            geom_point(aes(y=mean)) + geom_errorbar(aes(ymin = `2.5%`, ymax = `97.5%`), width = 0.5) +
            geom_linerange(aes(ymin = `25%`, ymax = `75%`), size = 1, show.legend = FALSE) +
            # facet_wrap(~ cluster, ncol = 1) +
            facet_grid(rows = vars(cluster), cols = vars(consensus.class), scales = "free_x", space = "free") +
            theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.45), axis.title.x = element_blank()) +
            theme(legend.title = element_blank()) +
            labs(y = "estimated occupancy within cluster")
    })

    # LSU plots

        # basic plot
        suppl.clusterocc.plot.list$LSU + scale_color_manual(values = plot.colors)
        ggsave(here("figures","FigS9_LSU_cluster_occupancy.pdf"), width = 12, height = 6, useDingbats = FALSE)

        # add domestic to help with annotating figure
        suppl.clusterocc.plot.list$LSU + aes(col = domestic)
        ggsave(here("figures","FigS9_LSU_cluster_occupancy_domestic.pdf"), width = 12, height = 6, useDingbats = FALSE)

        # add mammal >10kg to help with annotating figure
        suppl.clusterocc.plot.list$LSU + aes(col = (AdultBodyMass_g > 10000))
        ggsave(here("figures","FigS9_LSU_cluster_occupancy_10kg.pdf"), width = 12, height = 6, useDingbats = FALSE)

    # SSU plots

        # basic plot
        suppl.clusterocc.plot.list$SSU + scale_color_manual(values = plot.colors)
        ggsave(here("figures","FigS10_SSU_cluster_occupancy.pdf"), width = 12, height = 6, useDingbats = FALSE)

        # add domestic to help with annotating figure
        suppl.clusterocc.plot.list$SSU + aes(col = domestic)
        ggsave(here("figures","FigS10_SSU_cluster_occupancy_domestic.pdf"), width = 12, height = 6, useDingbats = FALSE)

        # add mammal >10kg to help with annotating figure
        suppl.clusterocc.plot.list$SSU + aes(col = (AdultBodyMass_g > 10000))
        ggsave(here("figures","FigS10_SSU_cluster_occupancy_10kg.pdf"), width = 12, height = 6, useDingbats = FALSE)

########################################################################################
