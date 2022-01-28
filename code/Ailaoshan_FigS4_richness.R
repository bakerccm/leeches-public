# R code to draw Fig S4

library("here")
library("tidyverse")
library("ggeffects")

# get data

    load(file=here("rdata","Ailaoshan_OTU_table.rdata"))

    # model summaries generated with Ailaoshan_model_summary.R
    model.summary <- list(
        LSU = readRDS(file = here("rds", "Ailaoshan_model_summary_200_final_LSU.rds")),
        SSU = readRDS(file = here("rds", "Ailaoshan_model_summary_200_final_SSU.rds"))
    )
    Nsite.output <- lapply(model.summary, FUN = function (X) as.data.frame(X$Nsite.output))
    rm(model.summary)

# calculate observed richness per Lab_ID
# note: uses all OTUs, not just those in common between the two datasets
# (though this gives qualitatively similar results)

    LabID.obsv.richness <- leech %>% select(dataset, Lab_ID, consensus.short, reads) %>%
        # calculate observed richness per Lab_ID
            mutate(occupied = ifelse(reads > 0, 1, 0)) %>%
            group_by(dataset, Lab_ID) %>%
            summarize(observed.richness = as.integer(sum(occupied, na.rm = TRUE)), .groups = "drop") %>%
        # convert to wide format
            pivot_wider(id_cols = Lab_ID, names_from = dataset, values_from = observed.richness, names_prefix = "observed.richness.") %>%
            filter(complete.cases(.)) %>%
        # add back leech_qty (shouldn't matter whether we use LSU or SSU as they are the same samples)
            left_join(leech %>% select(Lab_ID, leech_qty) %>% distinct(), by = "Lab_ID")

# FigS4a : observed richness per Lab_ID

    LabID.obsv.richness %>% filter(complete.cases(.)) %>%
        # ggplot(aes(x = observed.richness.LSU, y = observed.richness.SSU, colour = leech_qty)) +
        ggplot(aes(x = observed.richness.LSU, y = observed.richness.SSU)) +
        geom_jitter(width = 0.2, height=0.2, alpha=0.4) +
        geom_abline(slope=1, intercept=0, linetype="dotted") +
        # scale_color_gradient(low = "blue", high = "red") + labs(col="# leeches   ") +
        theme(aspect.ratio = 1) + labs(x="LSU observed species richness", y= "SSU observed species richness")
    ggsave(here("figures","FigS4a_obsvrich_LabID.pdf"), width=4, height=4)

    # test statistical significance of correlation
        temp <- LabID.obsv.richness %>%
            filter(complete.cases(.)) %>%
            select(observed.richness.LSU, observed.richness.SSU)
        cor.test(x = temp$observed.richness.LSU, y = temp$observed.richness.SSU, alternative = "greater")   # cor = 0.6490133; t = 21.173, df = 616, p-value < 2.2e-16
        rm(temp)

# FigS4b : observed richness per Lab_ID scatterplots with GLMs

    # does observed richness increase with number of leeches per replicate?

    # estimate GLM regressions
        # z-value in summary is "... z ratio if the dispersion is known (or fixed by the family)."
        # p-value in summary is "two-tailed p-value corresponding to the t or z ratio based on a Student t or Normal reference distribution."
        lsu.glm <- glm(observed.richness.LSU ~ log(leech_qty), data = LabID.obsv.richness, family=poisson())
        lsu.glm %>% summary() # z = 6.913, Pr(>|z|) = 4.76e-12 (--> halve this p-value for one-tailed test)
        ssu.glm <- glm(observed.richness.SSU ~ log(leech_qty), data = LabID.obsv.richness, family=poisson())
        ssu.glm %>% summary() # z = 9.988, Pr(>|z|) = < 2e-16 (--> halve this p-value for one-tailed test)

        lsu.predict <- data.frame(leech_qty = seq(1, 100, 1))
        lsu.predict.temp <- predict(lsu.glm, list(leech_qty = lsu.predict$leech_qty), type="response", se.fit = TRUE)
        lsu.predict$observed.richness <- lsu.predict.temp$fit
        lsu.predict$richness.lower <- lsu.predict.temp$fit - lsu.predict.temp$se.fit
        lsu.predict$richness.upper <- lsu.predict.temp$fit + lsu.predict.temp$se.fit
        rm(lsu.predict.temp)

        ssu.predict <- data.frame(leech_qty = seq(1, 100, 1))
        ssu.predict.temp <- predict(ssu.glm, list(leech_qty = ssu.predict$leech_qty), type="response", se.fit=TRUE)
        ssu.predict$observed.richness <- ssu.predict.temp$fit
        ssu.predict$richness.lower <- ssu.predict.temp$fit - ssu.predict.temp$se.fit
        ssu.predict$richness.upper <- ssu.predict.temp$fit + ssu.predict.temp$se.fit
        rm(ssu.predict.temp)

        richness.predict <- bind_rows(LSU  = lsu.predict, SSU = ssu.predict, .id="dataset")

    # what is predicted species richness in a replicate consisting of x leeches?

        # prediction (not confidence) intervals for 1:10 leeches
            ggpredict(lsu.glm, terms="leech_qty[1:10]", interval = "prediction")    # 0.33 (95% PI: [0.06, 1.78])
        # prediction (not confidence) intervals for 1:10 leeches
            ggpredict(ssu.glm, terms="leech_qty[1:10]", interval = "prediction")    # 0.29 (95% PI: [0.05, 1.73])

    # draw scatterplots
        LabID.obsv.richness %>%
            pivot_longer(c("observed.richness.LSU", "observed.richness.SSU"), names_to = "dataset", names_prefix = "observed.richness.", values_to = "observed.richness") %>%
            ggplot(aes(x = leech_qty, y = observed.richness)) +
                geom_jitter(width=0.25, height=0.2, alpha=0.33) + facet_wrap("dataset") +
                labs(x = "number of leeches per replicate ", y = "observed species richness") +
                geom_path(data= richness.predict, col= "blue", size=1) +
                geom_ribbon(data= richness.predict, aes(ymin = richness.lower, ymax = richness.upper),
                    fill= "blue", alpha = 0.4)
        ggsave(here("figures","FigS4b_richness_numleeches.pdf"), width = 7, height = 4)

    rm(lsu.glm, ssu.glm, lsu.predict, ssu.predict, richness.predict)

# calculate observed richness per Polygon_ID

    # determine which Lab_IDs are represented in both datasets
    Lab_IDs.to.keep <- leech %>% select(dataset, Lab_ID) %>% mutate(dummy = 0) %>% distinct() %>%
        pivot_wider(id_cols = Lab_ID, names_from = dataset, names_prefix = "dummy.", values_from = dummy) %>%
        filter(complete.cases(.)) %>% select(Lab_ID)

    # calculate per-Polygon_ID richness (only counts Lab_IDs shared by SSU and LSU datasets) -- all OTUs
    PolygonID.obsv.richness <- leech %>% select(dataset, Lab_ID, Polygon_ID, consensus.short, reads) %>%
        # filter to shared Lab_IDs only
            inner_join(Lab_IDs.to.keep, by = "Lab_ID") %>%
        # calculate observed richness per Polygon_ID
            mutate(occupied = ifelse(reads > 0, 1, 0)) %>%
            group_by(dataset, Polygon_ID, consensus.short) %>%
            summarize(occupied = ifelse(sum(occupied) > 0, 1, 0), .groups = "drop_last") %>%
            summarize(observed.richness = as.integer(sum(occupied, na.rm = TRUE)), .groups = "drop") %>%
        # convert to wide format
            pivot_wider(id_cols = Polygon_ID, names_from = dataset, values_from = observed.richness, names_prefix = "observed.richness.") %>%
            filter(complete.cases(.))

    # add counts of LabIDs (i.e. replicates, i.e. tubes) per PolygonID (only counts Lab_IDs shared by SSU and LSU datasets)
        LabIDs.per.PolygonID <- leech %>% inner_join(Lab_IDs.to.keep, by = "Lab_ID") %>% distinct(dataset, Polygon_ID, Lab_ID) %>%
            group_by(dataset, Polygon_ID) %>% summarize(replicates = n(), .groups = "drop") %>% distinct(Polygon_ID, replicates)
        PolygonID.obsv.richness <- PolygonID.obsv.richness %>% left_join(LabIDs.per.PolygonID, by = "Polygon_ID")

    rm(Lab_IDs.to.keep, LabIDs.per.PolygonID)

# FigS4c : observed richness per Polygon_ID

    PolygonID.obsv.richness %>% filter(complete.cases(.)) %>%
        # ggplot(aes(x = observed.richness.LSU, y = observed.richness.SSU, colour = replicates)) +
        ggplot(aes(x = observed.richness.LSU, y = observed.richness.SSU)) +
        geom_jitter(width = 0.2, height=0.2, alpha=0.4) +
        geom_abline(slope=1, intercept=0, linetype="dotted") +
        theme(aspect.ratio = 1) + labs(x="LSU observed species richness", y= "SSU observed species richness") +
        # scale_color_gradient(low = "blue", high = "red", breaks = seq(10, 30, 10)) + labs(col="# replicates") +
        expand_limits(x = 0, y = 0)
    ggsave(here("figures","FigS4c_obsvrich_PolygonID.pdf"), width=4, height=4)

    # test statistical significance of correlation
        temp <- PolygonID.obsv.richness %>% filter(complete.cases(.)) %>% select(observed.richness.LSU, observed.richness.SSU)
        cor.test(x = temp$observed.richness.LSU, y = temp$observed.richness.SSU, alternative = "greater")  # cor = 0.8851508; t = 20.839, df = 120, p-value < 2.2e-16
        rm(temp)

# FigS4d : observed richness per Polygon_ID scatterplots with GLMs

    # does observed richness increase with number of replicates per Polygon_ID?

    # estimate GLM regressions
        # z-value in summary is "... z ratio if the dispersion is known (or fixed by the family)."
        # p-value in summary is "two-tailed p-value corresponding to the t or z ratio based on a Student t or Normal reference distribution."
        lsu.glm <- glm(observed.richness.LSU ~ log(replicates), data = PolygonID.obsv.richness, family=poisson())
        lsu.glm %>% summary() # z = 10.197, Pr(>|z|) = < 2e-16 (--> halve this p-value for one-tailed test)
        ssu.glm <- glm(observed.richness.SSU ~ log(replicates), data = PolygonID.obsv.richness, family=poisson())
        ssu.glm %>% summary() # z = 14.937, Pr(>|z|) = < 2e-16 (--> halve this p-value for one-tailed test)

        lsu.predict <- data.frame(replicates = seq(0.5, 36, 0.5))
        lsu.predict.temp <- predict(lsu.glm, list(replicates = lsu.predict$replicates), type="response", se.fit=TRUE)
        lsu.predict$observed.richness <- lsu.predict.temp$fit
        lsu.predict$observed.richness.lower <- lsu.predict.temp$fit - lsu.predict.temp$se.fit
        lsu.predict$observed.richness.upper <- lsu.predict.temp$fit + lsu.predict.temp$se.fit
        rm(lsu.predict.temp)

        ssu.predict <- data.frame(replicates = seq(0.5, 36, 0.5))
        ssu.predict.temp <- predict(ssu.glm, list(replicates = ssu.predict$replicates), type="response", se.fit=TRUE)
        ssu.predict$observed.richness <- ssu.predict.temp$fit
        ssu.predict$observed.richness.lower <- ssu.predict.temp$fit - ssu.predict.temp$se.fit
        ssu.predict$observed.richness.upper <- ssu.predict.temp$fit + ssu.predict.temp$se.fit
        rm(ssu.predict.temp)

        richness.predict <- bind_rows(LSU  = lsu.predict, SSU = ssu.predict, .id="dataset")

    # draw scatterplots
        PolygonID.obsv.richness %>%
            pivot_longer(c("observed.richness.LSU", "observed.richness.SSU"), names_to = "dataset", names_prefix = "observed.richness.", values_to = "observed.richness") %>%
            ggplot(aes(x = replicates, y = observed.richness)) +
            geom_jitter() + facet_wrap("dataset") +
            labs(x = "number of replicates per patrol area", y = "observed species richness") +
            geom_path(data = richness.predict, col= "blue", size=1) +
            geom_ribbon(data = richness.predict, aes(ymin = observed.richness.lower, ymax = observed.richness.upper),
                fill= "blue", alpha = 0.4)
        ggsave(here("figures","FigS4d_observed_richness_sampling.pdf"), width = 7, height = 4)

    rm(lsu.glm, ssu.glm, lsu.predict, ssu.predict, richness.predict)

# calculate estimated richnesses in each dataset

    # get estimated richness values
        PolygonID.estimated.richness <- bind_rows(LSU = Nsite.output$LSU, SSU = Nsite.output$SSU, .id = "dataset") %>%
            select(dataset, Polygon_ID, everything())

    # add number of replicates (=tubes) per Polygon_ID in each dataset (actual replicates, not the dummy ones in leech.supplement)
        Lab_IDs.per.Polygon_ID <- leech %>%
            select(dataset, Polygon_ID, Lab_ID) %>% distinct() %>%
            group_by(dataset, Polygon_ID) %>% summarize(replicates = n(), .groups = "drop")
        PolygonID.estimated.richness <- PolygonID.estimated.richness %>%
            left_join(Lab_IDs.per.Polygon_ID, by = c("dataset", "Polygon_ID"))

    # reformat as wide data and calculate mean number of replicates
    # (use mean for Fig S2e because models use Lab_IDs regardless of whether they match across LSU/SSU, so these counts may differ)
        PolygonID.estimated.richness <- PolygonID.estimated.richness %>%
            pivot_wider(id_cols = Polygon_ID, names_from = dataset, names_sep = ".", values_from = c("mean", "replicates")) %>%
            rowwise() %>% mutate(replicates.mean = mean(c(replicates.LSU, replicates.SSU))) %>% ungroup()

    rm(Lab_IDs.per.Polygon_ID)

# FigS4e : estimated richness per Polygon_ID

    PolygonID.estimated.richness %>% filter(complete.cases(.)) %>%
        rename(`LSU estimated species richness` = mean.LSU, `SSU estimated species richness` = mean.SSU) %>%
        ggplot(aes(x = `LSU estimated species richness`, y = `SSU estimated species richness`)) +
        geom_point(alpha=0.4) +
        geom_abline(slope=1, intercept=0, linetype="dotted") +
        theme(aspect.ratio = 1) + expand_limits(x = 13, y = 13)
    ggsave(here("figures","FigS4e_estimated_richness_PolygonID.pdf"), width=4, height=4)

    # test statistical significance of correlation
        PolygonID.estimated.richness %>% filter(complete.cases(.)) %>%
            cor.test(~ mean.LSU + mean.SSU, data = ., alternative = "greater") # cor = 0.8717483, t = 19.491, df = 120, p-value < 2.2e-16

# FigS4f : estimated richness per Polygon_ID scatterplots

    # does estimated richness increase with number of replicates per Polygon_ID?

    # regressions are non-significant
    # (note use of linear regression rather than poisson GLM since estimates are not integer like observed values are)

        lsu.lm <- lm(mean.LSU ~ log(replicates.LSU), data = PolygonID.estimated.richness)
        lsu.lm %>% summary() # F-statistic: 0.005516 on 1 and 124 DF,  p-value: 0.9409
        ssu.lm <- lm(mean.SSU ~ log(replicates.SSU), data = PolygonID.estimated.richness)
        ssu.lm %>% summary() # F-statistic: 1.493 on 1 and 125 DF,  p-value: 0.224

    # so draw scatterplots without LOESS or regression lines

        PolygonID.estimated.richness %>%
            pivot_longer(c("mean.LSU", "mean.SSU"), names_to = "dataset", names_prefix = "mean.", values_to = "estimated richness") %>%
            mutate(replicates = ifelse(dataset == "LSU", replicates.LSU, replicates.SSU)) %>%
            select(-replicates.mean, -replicates.LSU, -replicates.SSU) %>%
            filter(complete.cases(.)) %>%
            ggplot(aes(x = replicates,  y= `estimated richness`)) +
            geom_jitter() + facet_wrap("dataset") +
            labs(x = "number of replicates per patrol area", y = "estimated species richness")
        ggsave(here("figures","FigS4f_estimated_richness_sampling.pdf"), width = 7, height = 4)
