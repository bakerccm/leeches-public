# R code to draw Fig 5 and Fig S8

library("here")
library("vegan")
library("sf")
library("raster")
library("sp")
library("tidyverse")
library("scales") # to get hue_pal()
library("dendextend") # for drawing trees
library("rcompanion") # for Cramer's V

########################################################################################
# averaged jaccard distances from occupancy model
jaccard.site <- list(
    LSU = readRDS(here("rds", "Ailaoshan_final_LSU_jaccard_site.rds")),
    SSU = readRDS(here("rds", "Ailaoshan_final_SSU_jaccard_site.rds"))
)

# data generated with Ailaoshan_OTU_table.Rmd
load(file=here("rdata","Ailaoshan_OTU_table.rdata"))

# Ailaoshan processed GIS data
load(file=here("rdata","Ailaoshan_gis.rdata"))

# environmental data
load(here("rdata","Ailaoshan_environmental.rdata"))

########################################################################################
# good match in distances across datasets

    mantel(as.dist(1-jaccard.site$LSU), as.dist(1-jaccard.site$SSU)) # Mantel statistic r: 0.9447, p= 0.001

# make trees for sites

    jaccard.sites.tree <- lapply(jaccard.site, function (X) hclust(as.dist(1-X), method="ward.D2"))

# split sites trees into groups

    no.sites.groups <- 3

    jaccard.sites.cut <- lapply(jaccard.sites.tree, function (X) {
        data.frame(jaccard.cut = factor(cutree(X, k = no.sites.groups))) %>%
            rownames_to_column("Polygon_ID")
    })

    # if necessary, swap groups here to ensure group numbers match up between LSU and SSU clusters
    # swap groups in SSU so 1 = high elevation and 3 = low elevation, as in LSU
    jaccard.sites.cut$SSU <- jaccard.sites.cut$SSU %>%
        mutate(jaccard.cut = recode(as.character(jaccard.cut), "2" = "3",  "3" = "2")) %>%
        mutate(jaccard.cut = factor(jaccard.cut))
    jaccard.sites.cut.df <- full_join(jaccard.sites.cut$LSU %>% rename(jaccard.cut.LSU = jaccard.cut), jaccard.sites.cut$SSU %>% rename(jaccard.cut.SSU = jaccard.cut), by = "Polygon_ID")

    # some of the group 1 sites in the LSU dataset are group 2 in the SSU dataset
    (jaccard.sites.cut.table <- jaccard.sites.cut.df %>% select(jaccard.cut.LSU, jaccard.cut.SSU) %>% table())

    # test how well clustering matches across the two datasets using Cramer's V and chi-squared test

    cramerV(jaccard.sites.cut.table, ci = TRUE)
    # Cramer.V lower.ci upper.ci
    # 1   0.7922   0.7337   0.8544

    chisq.test(jaccard.sites.cut.table) # X-squared = 215.9, df = 4, p-value < 2.2e-16

# draw tree plots

    # rotating and coloring branches here to match other figures
    ## double check this is correct if the dataset changes! ##

    tree.filesnames <- list(LSU = here("figures", "FigS8a_tree_sites_LSU.pdf"), SSU = here("figures", "FigS8b_tree_sites_SSU.pdf"))
    jaccard.sites.tree <- lapply(jaccard.sites.tree, as.dendrogram)
    # rotate clusters 2 and 3 for LSU tree only
        jaccard.sites.tree$LSU <- all_couple_rotations_at_k(jaccard.sites.tree$LSU, k = 3)[[2]]
    # draw plots
        lapply(names(jaccard.sites.tree), function(X) {
            jaccard.sites.tree[[X]] %>%
                set("branches_k_color", hue_pal()(3), k=3) %>%
                set("labels_cex", 0.3) %>% set("branches_lwd", 0.5) %>%
                as.ggdend() %>% ggplot(horiz = TRUE)
            ggsave(tree.filesnames[[X]], width=4, height=12, useDingbats = FALSE)
        })

# ordinations for sites

    set.seed(1)
    nmds.sites <- lapply(jaccard.site, function (X) metaMDS(as.dist(1-X), k=2, trymax = 200))

    lapply(nmds.sites, function (X) X$stress) # LSU: 0.08808181; SSU: 0.1058787

    sites.env <- sapply(names(jaccard.site), simplify = FALSE, function (X) {
        tibble(Polygon_ID = rownames(jaccard.site[[X]])) %>%
            left_join(env.data %>% select(Polygon_ID, elevation_median, distance_to_nature_reserve_boundary), by = "Polygon_ID")
    })

    # draw NMDS plots
    nmds.filesnames <- list(LSU = here("figures", "Fig5a_nmds_sites_LSU.pdf"), SSU = here("figures", "Fig5b_nmds_sites_SSU.pdf"))
    # note: hue_pal()(3) gives the default ggplot colours so we can match the colours in other figures
    lapply(names(nmds.sites), function(X) {
        pdf(nmds.filesnames[[X]], width=4, height=4, useDingbats = FALSE)
        par(mar = c(4,4,0,0) + 0.1)
        plot(nmds.sites[[X]], display = "sites", type="n", xlab = "", ylab = "") # set up plot area
        points(nmds.sites[[X]], col = hue_pal()(3)[jaccard.sites.cut[[X]]$jaccard.cut], pch=16)
        ordisurf(nmds.sites[[X]], sites.env[[X]]$elevation_median, add=TRUE, col="red", labcex = 0) # set labels manually
        ordisurf(nmds.sites[[X]], sites.env[[X]]$distance_to_nature_reserve_boundary, add=TRUE, col= "blue", labcex = 0) # set labels manually
        mtext(text = "NMDS1", side = 1, line = 2.2)
        mtext(text = "NMDS2", side = 2, line = 2.2)
        dev.off()
    })

    # draw additional NMDS plots with labelled contours
    nmds.labels.filesnames <- list(LSU = here("figures", "Fig5a_nmds_sites_LSU_labels.pdf"), SSU = here("figures", "Fig5b_nmds_sites_SSU_labels.pdf"))
    # note: hue_pal()(3) gives the default ggplot colours so we can match the colours in other figures
    lapply(names(nmds.sites), function(X) {
        pdf(nmds.labels.filesnames[[X]], width=4, height=4, useDingbats = FALSE)
        par(mar = c(4,4,0,0) + 0.1)
        plot(nmds.sites[[X]], display = "sites", type="n", xlab = "", ylab = "") # set up plot area
        points(nmds.sites[[X]], col = hue_pal()(3)[jaccard.sites.cut[[X]]$jaccard.cut], pch=16)
        ordisurf(nmds.sites[[X]], sites.env[[X]]$elevation_median, add=TRUE, col="red")
        ordisurf(nmds.sites[[X]], sites.env[[X]]$distance_to_nature_reserve_boundary, add=TRUE, col= "blue")
        mtext(text = "NMDS1", side = 1, line = 2.2)
        mtext(text = "NMDS2", side = 2, line = 2.2)
        legend(x="topright", legend = c("high", "intermediate", "low"), col = hue_pal()(3), pch = 16)
        legend(x="bottomleft", legend = c("elevation","distance to reserve edge"), col = c("red","blue"), lwd = 1, seg.len = 1, cex = 0.9, y.intersp=0.8, x.intersp=0.8, text.width=0.525)
        dev.off()
    })

# plot maps

    ailaoshan.polygons.jaccard <- ailaoshan.polygons %>%
        mutate(Polygon_ID = as.character(Polygon_ID)) %>%
        left_join(jaccard.sites.cut.df, by = "Polygon_ID") %>%
        mutate(jaccard.cut.LSU = recode(jaccard.cut.LSU, "1" = "high", "2" = "intermediate", "3" = "low")) %>%
        mutate(jaccard.cut.SSU = recode(jaccard.cut.SSU, "1" = "high", "2" = "intermediate", "3" = "low"))

    # draw maps
    ailaoshan.polygons.jaccard %>% ggplot() + geom_sf(aes(fill=jaccard.cut.LSU), size = 0.1, color = 'black') +
        theme(legend.position=c(0,0), legend.justification=c(0,0), legend.background = element_rect(fill=NA)) + # scale_fill_brewer(direction = -1) +
        coord_sf() + labs(fill = "LSU site clusters")
    ggsave(here("figures", "Fig5c_map_sites_LSU.pdf"), width=4, height=4)
    ailaoshan.polygons.jaccard %>% ggplot() + geom_sf(aes(fill=jaccard.cut.SSU), size = 0.1, color = 'black') +
        theme(legend.position=c(0,0), legend.justification=c(0,0), legend.background = element_rect(fill=NA)) +
        coord_sf() + labs(fill = "SSU site clusters")
    ggsave(here("figures", "Fig5d_map_sites_SSU.pdf"), width=4, height=4)

########################################################################################

# export site cluster membership
# jaccard.cut: 1 = high elevation sites; 2 = intermediate sites; 3 = low elevation sites
    saveRDS(jaccard.sites.cut$LSU, file = here("rds", "Ailaoshan_jaccard_sites_cut_LSU.rds"))
    saveRDS(jaccard.sites.cut$SSU, file = here("rds", "Ailaoshan_jaccard_sites_cut_SSU.rds"))

########################################################################################
