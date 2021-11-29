# R code to investigate whether missing data are missing at random

library("here")
library("ggfortify") # for PCA plot, see https://cran.r-project.org/web/packages/ggfortify/vignettes/plot_pca.html
library("aod") # for Wald test
library("MASS")
library("klaR")
library("tidyverse")

# get data

    load(file=here("rdata","Ailaoshan_OTU_table.rdata"))

# prepare data for analysis

    replicates <- bind_rows(leech, leech.supplement) %>%
        # remove extraneous variables (mostly to do with OTUs)
            select(-OTU, -reads, -replicate_no, -starts_with("consensus"), -Chinese_common_name, -domestic, -starts_with("AdultBodyMass"), -ends_with("_std_dev"), -ends_with("_mean")) %>%
            distinct() %>%
        # put data from the two different datasets into different columns
            pivot_wider(names_from = dataset, values_from = c(leech_qty, fraction.reads.humans)) %>%
        # summarize at Polygon_ID level
            group_by_at(vars(-starts_with("leech_qty"), -starts_with("fraction.reads.humans"), -Lab_ID)) %>%
            summarize(no.replicates_LSU = sum(!is.na(leech_qty_LSU)),
                no.replicates_SSU = sum(!is.na(leech_qty_SSU)),
                mean.leech_qty_LSU = ifelse(sum(!is.na(leech_qty_LSU)) == 0, NA, mean(leech_qty_LSU, na.rm = TRUE)),
                mean.leech_qty_SSU = ifelse(sum(!is.na(leech_qty_SSU)) == 0, NA, mean(leech_qty_LSU, na.rm = TRUE)),
                mean.fraction.reads.humans_LSU = ifelse(sum(!is.na(leech_qty_LSU)) == 0, NA, mean(fraction.reads.humans_LSU, na.rm = TRUE)),
                mean.fraction.reads.humans_SSU = ifelse(sum(!is.na(leech_qty_SSU)) == 0, NA, mean(fraction.reads.humans_LSU, na.rm = TRUE)), .groups = "drop") %>%
        # add dummies to show whether Polygon_IDs have data or not
            mutate(sequence.data = !is.na(Ranger_ID), polygon.data = !is.na(elevation_median)) %>%
        # rearrange slightly
            select(Polygon_ID, sequence.data, polygon.data, no.replicates_LSU, no.replicates_SSU, mean.leech_qty_LSU, mean.leech_qty_SSU,
                mean.fraction.reads.humans_LSU, mean.fraction.reads.humans_SSU, everything())

        replicates %>% select(sequence.data, polygon.data) %>% table()

        names(replicates)

# is polygon data missing at random from the sequence data?
#  -- plots

    sequence_data.plot <- replicates %>% filter(sequence.data == TRUE) %>% ggplot(aes(x = polygon.data))

    sequence_data.plot + geom_boxplot(aes(y = no.replicates_LSU))
    sequence_data.plot + geom_boxplot(aes(y = no.replicates_SSU))

    sequence_data.plot + geom_boxplot(aes(y = mean.leech_qty_LSU))
    sequence_data.plot + geom_boxplot(aes(y = mean.leech_qty_SSU))

    sequence_data.plot + geom_boxplot(aes(y = mean.fraction.reads.humans_LSU))
    sequence_data.plot + geom_boxplot(aes(y = mean.fraction.reads.humans_SSU))

    rm(sequence_data.plot)

# is polygon data missing at random from the sequence data?
#  -- tests

    replicates %>% filter(sequence.data == TRUE) %>%
        kruskal.test(no.replicates_LSU ~ polygon.data, data = .)  # Kruskal-Wallis chi-squared = 0.49155, df = 1, p-value = 0.4832
    replicates %>% filter(sequence.data == TRUE) %>%
        kruskal.test(no.replicates_SSU ~ polygon.data, data = .)  # Kruskal-Wallis chi-squared = 0.047004, df = 1, p-value = 0.8284

    replicates %>% filter(sequence.data == TRUE) %>%
        kruskal.test(mean.leech_qty_LSU ~ polygon.data, data = .)  # Kruskal-Wallis chi-squared = 0.025819, df = 1, p-value = 0.8723
    replicates %>% filter(sequence.data == TRUE) %>%
        kruskal.test(mean.leech_qty_SSU ~ polygon.data, data = .)  # Kruskal-Wallis chi-squared = 0.12143, df = 1, p-value = 0.7275

    replicates %>% filter(sequence.data == TRUE) %>%
        kruskal.test(mean.fraction.reads.humans_LSU ~ polygon.data, data = .)  # Kruskal-Wallis chi-squared = 0.12076, df = 1, p-value = 0.7282
    replicates %>% filter(sequence.data == TRUE) %>%
        kruskal.test(mean.fraction.reads.humans_SSU ~ polygon.data, data = .)  # Kruskal-Wallis chi-squared = 0.00030374, df = 1, p-value = 0.9861

# is sequence data missing at random from the polygons
#  -- plots

    polygons.plot <- replicates %>% filter(polygon.data == TRUE) %>% ggplot(aes(x = sequence.data))

    polygons.plot + geom_boxplot(aes(y = shape_perimeter))
    polygons.plot + geom_boxplot(aes(y = shape_area_ha))

    polygons.plot + geom_boxplot(aes(y = longitude))
    polygons.plot + geom_boxplot(aes(y = latitude))

    polygons.plot + geom_boxplot(aes(y = elevation_median))
    polygons.plot + geom_boxplot(aes(y = tpi_median))
    polygons.plot + geom_boxplot(aes(y = distance_to_road_median))
    polygons.plot + geom_boxplot(aes(y = distance_to_stream_median))
    polygons.plot + geom_boxplot(aes(y = distance_to_nature_reserve_boundary))

    rm(polygons.plot)

# is sequence data missing at random from the polygons?
#  -- tests

    replicates %>% filter(polygon.data == TRUE) %>%
        kruskal.test(shape_perimeter ~ sequence.data, data = .)  # Kruskal-Wallis chi-squared = 3.4269, df = 1, p-value = 0.06414
    replicates %>% filter(polygon.data == TRUE) %>%
        kruskal.test(shape_area_ha ~ sequence.data, data = .)  # Kruskal-Wallis chi-squared = 3.4366, df = 1, p-value = 0.06377

    replicates %>% filter(polygon.data == TRUE) %>%
        kruskal.test(longitude ~ sequence.data, data = .)  # Kruskal-Wallis chi-squared = 0.39707, df = 1, p-value = 0.5286
    replicates %>% filter(polygon.data == TRUE) %>%
        kruskal.test(latitude ~ sequence.data, data = .)  # Kruskal-Wallis chi-squared = 0.69569, df = 1, p-value = 0.4042

    replicates %>% filter(polygon.data == TRUE) %>%
        kruskal.test(elevation_median ~ sequence.data, data = .)  # Kruskal-Wallis chi-squared = 5.5014, df = 1, p-value = 0.019
    replicates %>% filter(polygon.data == TRUE) %>%
        kruskal.test(tpi_median ~ sequence.data, data = .)  # Kruskal-Wallis chi-squared = 0.015313, df = 1, p-value = 0.9015
    replicates %>% filter(polygon.data == TRUE) %>%
        kruskal.test(distance_to_road_median ~ sequence.data, data = .)  # Kruskal-Wallis chi-squared = 9.6492, df = 1, p-value = 0.001894
    replicates %>% filter(polygon.data == TRUE) %>%
        kruskal.test(distance_to_stream_median ~ sequence.data, data = .)  # Kruskal-Wallis chi-squared = 0.52278, df = 1, p-value = 0.4697
    replicates %>% filter(polygon.data == TRUE) %>%
        kruskal.test(distance_to_nature_reserve_boundary ~ sequence.data, data = .)  # Kruskal-Wallis chi-squared = 3.2079, df = 1, p-value = 0.07328

# do environmental values vary between replicates with and without sequence data?

    # mean environmental values for replicates with and without sequence data
    replicates %>% filter(polygon.data == TRUE) %>% group_by(sequence.data) %>%
        summarize_at(vars(elevation_median, tpi_median, distance_to_road_median, distance_to_nature_reserve_boundary), .fun=mean, .groups = "drop")

    # PCA for environmental values, with and without sequence data
    my.data <- replicates %>% filter(polygon.data == TRUE) %>%
        select(sequence.data, elevation_median, tpi_median, distance_to_road_median, distance_to_stream_median, distance_to_nature_reserve_boundary) %>%
        mutate(elevation_median = scale(elevation_median),
            tpi_median = scale(tpi_median),
            distance_to_road_median = scale(distance_to_road_median),
            distance_to_stream_median = scale(distance_to_stream_median),
            distance_to_nature_reserve_boundary = scale(distance_to_nature_reserve_boundary))
    my.pca <- prcomp(my.data %>% select(-sequence.data))
    autoplot(my.pca, data=my.data, colour = "sequence.data", loadings = TRUE, loadings.colour = 'blue',
             loadings.label = TRUE, loadings.label.size = 3, loadings.label.colour = "blue")

# histograms for environmental variables, with and without sequence data

    replicates %>% filter(polygon.data == TRUE) %>% ggplot(aes(x = shape_perimeter)) + geom_histogram() + facet_wrap(~sequence.data, ncol = 1)
    replicates %>% filter(polygon.data == TRUE) %>% ggplot(aes(x = shape_area_ha)) + geom_histogram() + facet_wrap(~sequence.data, ncol = 1)
    replicates %>% filter(polygon.data == TRUE) %>% ggplot(aes(x = elevation_median)) + geom_histogram() + facet_wrap(~sequence.data, ncol = 1)
    replicates %>% filter(polygon.data == TRUE) %>% ggplot(aes(x = distance_to_road_median)) + geom_histogram() + facet_wrap(~sequence.data, ncol = 1)
    replicates %>% filter(polygon.data == TRUE) %>% ggplot(aes(x = distance_to_nature_reserve_boundary)) + geom_histogram() + facet_wrap(~sequence.data, ncol = 1)

# do values of environmental variables predict whether replicate has sequence data or not?
# logit regression
# see https://stats.idre.ucla.edu/r/dae/logit-regression/

    replicates %>% names()

    mylogit <- glm(sequence.data ~ elevation_median + tpi_median + distance_to_road_median + distance_to_stream_median + distance_to_nature_reserve_boundary, data = replicates %>% filter(polygon.data), family = "binomial")

    summary(mylogit)

    # CIs using profiled log-likelihood
    confint(mylogit)

    # CIs using standard errors
    confint.default(mylogit)

    wald.test(b = coef(mylogit), Sigma = vcov(mylogit), Terms = 2:6)

# do environmental variables predict whether replicate has sequence data or not?
#  -- LDA approach

    training_sample <- sample(c(TRUE, FALSE), nrow(replicates[replicates$polygon.data,]), replace = T, prob = c(0.6,0.4))
    train <- replicates[replicates$polygon.data,][training_sample, ]
    test <- replicates[replicates$polygon.data,][!training_sample, ]

    lda.sequence.data <- lda(sequence.data ~ shape_perimeter + shape_area_ha +longitude + latitude + elevation_median + tpi_median +
        distance_to_road_median + distance_to_stream_median + distance_to_nature_reserve_boundary, data=train)
    lda.sequence.data #show results
    plot(lda.sequence.data, col = as.integer(train$sequence.data))

    plot(lda.sequence.data, dimen = 1, type = "b")

    partimat(factor(sequence.data) ~ elevation_median + tpi_median + distance_to_road_median + distance_to_stream_median + distance_to_nature_reserve_boundary, data=train, method="lda")
