# R code to draw Fig 2ab

library("here")
library("tidyverse")

# filenames

    input.otu.table <- here("rdata","Ailaoshan_OTU_table.rdata")
    input.Ntotal.summary <- here("rds", "Ailaoshan_model_summary_Ntotal.rds")

# get data

    load(file = input.otu.table)

    Ntotal.summary <- readRDS(input.Ntotal.summary)

# colorblind-friendly plot colors from Wong 2011 (https://www.nature.com/articles/nmeth.1618)

    plot.colors <- c("amphibians" = rgb(0, 158, 115, maxColorValue = 255), # bluish green
                     "birds" = rgb(86, 180, 233, maxColorValue = 255), # sky blue
                     "mammals" = rgb(230, 159, 0, maxColorValue = 255), # orange
                     "squamates" = rgb(0, 0, 0, maxColorValue = 255)) # black

# Fig 2a: plot breakdown of detected OTUs by taxonomic class

    leech %>%
        distinct(dataset, consensus.class, consensus.short) %>%
        group_by(dataset, consensus.class) %>%
        mutate(consensus.class = tolower(consensus.class)) %>%
        # relabel (non-avian) reptiles as squamates
            mutate(consensus.class = ifelse(consensus.class == "reptiles", "squamates", consensus.class)) %>%
        ggplot(aes(x = consensus.class, fill=consensus.class)) + stat_count(show.legend = FALSE) + facet_wrap("dataset") +
        labs(y = "number of species") + scale_fill_manual(values = plot.colors) +
        theme(axis.text.x = element_text(angle=45, hjust=1), axis.title.x = element_blank())

    ggsave(filename = here("figures","Fig2a_taxonomic_breakdown.pdf"), width = 4, height=3)

# Fig 2b: draw Ntotal results using ggplot

    # count OTUs actually detected
        detected.OTUs <- leech %>%
            distinct(dataset, consensus.class, consensus.short) %>%
            group_by(dataset) %>% summarize(count = n())

    Ntotal.summary %>%
        mutate(M = as.factor(M)) %>%
        ggplot(aes(x=M, y= mean)) + facet_wrap(~dataset) + expand_limits(y = 0) +
        geom_point(aes(y=mean)) +
        geom_linerange(aes(ymin = `25%`, ymax = `75%`), size = 1, show.legend = FALSE) +
        geom_errorbar(aes(ymin = `2.5%`, ymax = `97.5%`), size = 0.5, width = 0.1) +
        geom_hline(data = detected.OTUs, aes(yintercept = count), linetype = "dashed") +
        labs(x = "supercommunity size (M)", y = "estimated species richness") + guides(linetype = 'none')

    ggsave(here("figures","Fig2b_estimated_richness.pdf"), width = 3.2, height = 3.0)
