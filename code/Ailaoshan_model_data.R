# R code to prepare data for occupancy modelling

library("tidyverse")

# get command line arguments

    args <- commandArgs(trailingOnly = TRUE)
    M <- args[1]    # M is total number of species in augmented dataset i.e. including observed species
    input.filename <- args[2] # rdata/Ailaoshan_OTU_table.rdata
    output.filename <- args[3]

# get raw data

    load(file = input.filename)

# augment read data with extra Polygon_IDs and convert reads to true/false

    leech.full <- sapply(c("LSU", "SSU"), simplify = FALSE, FUN = function (X) {
        leech %>%
            # add extra rows for Polygon_IDs with missing data
                bind_rows(leech.supplement) %>%
            # split by dataset
                filter(dataset == X) %>%
            # convert read counts to true/false
                mutate(detected = (reads>0)) %>%
                select(-reads) %>%
            # scale variables
                mutate(leech_qty_scaled = scale(leech_qty)) %>%
            # convert class to lower case and replace reptiles with squamates
            # (should move this to OTU table code)
                mutate(consensus.class = tolower(consensus.class)) %>%
                mutate(consensus.class = ifelse(consensus.class == "reptiles", "squamates", consensus.class))
    })

######################################################
# augment data with extra species

    # observed species in LSU and SSU datasets (n)
    species.observed <- lapply(leech.full, function (X) {
        X %>% select(consensus.class, OTU) %>% distinct() %>%
            group_by(consensus.class) %>%
            summarize(observed = n(), .groups = "drop")
    })

    # species breakdowns to use depending on M supplied on command line
    M.species <- list(
        `200` = c("mammals" = 48, "birds" = 115, "squamates" = 14, "amphibians" = 23),
        `150` = c("mammals" = 41, "birds" = 77, "squamates" = 9, "amphibians" = 23),
        `100` = c("mammals" = 34, "birds" = 39, "squamates" = 5, "amphibians" = 22)
    )
    species.inventory <- M.species[[M]] %>% as_tibble_row() %>%
        pivot_longer(cols = everything(), names_to = "consensus.class", values_to = "M")

    # augment is number of extra species we need to augment each dataset by
    species.augment <- lapply(species.observed, function (X) {
        full_join(X, species.inventory, by = "consensus.class") %>%
            mutate(augment = M - observed) %>%
            select(consensus.class, observed, augment, M)
    })

    # data (i.e. zeros) to add to leech.full for modelling
    leech.zeros <- sapply(names(leech.full), simplify = FALSE, FUN = function (X) {
        leech.full[[X]] %>% select(dataset, Polygon_ID, replicate_no, leech_qty, elevation_median, tpi_median, distance_to_road_median, distance_to_stream_median, distance_to_nature_reserve_boundary) %>% distinct() %>%
            mutate(detected = FALSE) %>%
            group_by(dataset, Polygon_ID, replicate_no, leech_qty, detected, elevation_median, tpi_median, distance_to_road_median, distance_to_stream_median, distance_to_nature_reserve_boundary) %>%
            expand(consensus.class = species.augment[[X]]$consensus.class) %>%
            inner_join(species.augment[[X]] %>% select(consensus.class, augment) %>% filter(augment > 0), by = "consensus.class") %>%
            group_by(dataset, Polygon_ID, replicate_no, leech_qty, detected, consensus.class, augment, elevation_median, tpi_median, distance_to_road_median, distance_to_stream_median, distance_to_nature_reserve_boundary) %>%
            expand(OTU = paste0(X, "_", consensus.class,"_", sprintf("%03d", 1:augment))) %>%
            ungroup() %>% select(-augment)
    })

    leech.augmented <- sapply(names(leech.full), simplify = FALSE, FUN = function (X) bind_rows(detected = leech.full[[X]], augmented = leech.zeros[[X]], .id = "OTU.type") )

# make empty lists of arrays to be filled with data

    # get max dimensions for sizing data arrays

        # number of sites i.e. Polygon_IDs for each dataset
        num.sites <- lapply(leech.augmented,
            FUN = function (X) X %>% select(Polygon_ID) %>% distinct() %>% nrow() )

        # (maximum) number of replicates per Polygon_ID in each dataset
        num.reps <- lapply(leech.augmented,
            FUN = function (X) X %>% group_by(Polygon_ID, OTU) %>% tally() %>% group_by(Polygon_ID) %>% summarize(max.n = max(n), .groups = "drop") %>% pull(max.n) %>% max() )

        # total number of observed + augmented species for each dataset
        num.species <- lapply(leech.augmented,
            FUN = function (X) X %>% select(OTU) %>% distinct() %>% nrow() )

    # make empty data arrays to be filled later

        # species detections for each site, to be indexed as [num.sites, num.reps, num.species]
        y <- sapply(names(leech.augmented), simplify = FALSE, FUN = function (X) array(NA, dim=c(num.sites[[X]], num.reps[[X]], num.species[[X]])) )

        # number of leeches per replicate, to be indexed as [num.sites, num.reps]
        numleeches <- sapply(names(leech.augmented), simplify = FALSE, FUN = function (X) array(NA, dim=c(num.sites[[X]], num.reps[[X]])) )

        # group for each species (mammals/birds, amphibians/squamates)
        g <- sapply(names(leech.augmented), simplify = FALSE, FUN = function (X) rep(NA, num.species[[X]]) )

    # check data dimensions
    # (note that num.sites now includes Polygon_IDs that never had any reads but also extra Polygon_IDs
    # derived from Ranger_IDs, so exceeds number of actual polygons
        unlist(num.sites) #LSU: 209, SSU: 209
        unlist(num.reps) #LSU: 40, SSU: 38
        unlist(num.species) #LSU: 474, SSU: 474 (cf. observed species only = LSU: 59, SSU: 72)

# prepare data

    # names for ensuring that arrays are aligned

        # species
        # (note sorting by desc(OTU.type) puts all the detected species first and the augmented species last)
        OTU.labels <- lapply(leech.augmented,
            FUN = function (X) X %>% select(OTU.type, OTU) %>% distinct() %>% arrange(desc(OTU.type), OTU) %>% pull(OTU) )

        # sites
        polygon.labels <- lapply(leech.augmented,
            FUN = function (X) X %>% select(Polygon_ID) %>% distinct() %>% arrange(Polygon_ID) %>% pull(Polygon_ID) )

        # replicates
        replicate.labels <- lapply(leech.augmented,
            FUN = function (X) X %>% select(replicate_no) %>% distinct() %>% arrange(replicate_no) %>% pull(replicate_no) )

    num.species.groups <- list()    # 2 -- i.e. mammals/birds, amphibians/squamates
    species.groups <- list()        # i.e. mammals & birds = 1, amphibians & squamates = 2
    species.group.labels <- list()  # i.e. mammals/birds, amphibians/squamates

    for (i in names(leech.augmented)) {

        # fill y[[i]] with data
            temp <- leech.augmented[[i]] %>%
                select(Polygon_ID, OTU, replicate_no, detected) %>%
                pivot_wider(names_from = replicate_no, values_from = detected)
            for (k in seq_along(OTU.labels[[i]])) {
                # this is just the wide data from before
                    temp.k <- temp %>% filter(OTU == OTU.labels[[i]][k]) # filter to the kth OTU
                # arrange rows and then remove Polygon_ID column
                    temp.k <- temp.k[match(polygon.labels[[i]], temp.k$Polygon_ID),] %>% select(starts_with("replicate_"))
                # arrange columns
                    temp.k <- temp.k[,match(replicate.labels[[i]], names(temp.k))]
                # convert to numeric matrix
                    y[[i]][,,k] <- temp.k %>% data.matrix()
            }
            rm(temp, k, temp.k)
            dimnames(y[[i]]) <- list(polygon.labels[[i]], replicate.labels[[i]], OTU.labels[[i]])

        # fill numleeches[[i]] with data
            temp <- leech.augmented[[i]] %>%
                select(dataset, Polygon_ID, replicate_no, leech_qty) %>%
                distinct() %>% pivot_wider(names_from = replicate_no, values_from = leech_qty) ### N.B. using unscaled leech_qty ###
            temp <- temp[match(polygon.labels[[i]], temp$Polygon_ID),]
            temp <- temp %>% select(starts_with("replicate_"))
            temp <- temp[,match(replicate.labels[[i]], names(temp))]
            numleeches[[i]] <- temp %>% data.matrix()
            rm(temp)
            dimnames(numleeches[[i]]) <- list(polygon.labels[[i]], replicate.labels[[i]])

        # fill g[[i]] with species group data in numeric form (i.e. group that each species belongs to: mammals/birds, amphibians/squamates)
            num.species.groups[[i]] <- 2
            species.groups[[i]] <- list(mammals = 1, birds = 1, amphibians = 2, squamates = 2)
            species.group.labels[[i]] <- c("mammals/birds", "amphibians/squamates") # should be consistent with species.groups[[i]] above
            temp <- leech.augmented[[i]] %>% select(dataset, OTU, consensus.class) %>% distinct()
            temp <- temp[match(OTU.labels[[i]], temp$OTU),]
            g[[i]] <- temp %>%
                mutate(g = recode(consensus.class, !!!species.groups[[i]])) %>%
                mutate(g = as.integer(g)) %>% pull(g)
            rm(temp)
            names(g[[i]]) <- OTU.labels[[i]]

    }

    rm(i, OTU.labels, polygon.labels, replicate.labels)

# observed occurrence values as starting values for z

    # note that these vary by dataset

    # note also that we have some site/species combinations that don't even have one replicate
    #    any(is.na(y$LSU[,1,]))
    #    any(is.na(y$SSU[,1,]))

    # combinations of site,species that had at least one detection
    z.start <- lapply(y, FUN = function (X) {
        suppressWarnings( # for the 'all NA' sites, max returns -Inf with a warning
            zst <- apply(X, MARGIN = c(1,3), FUN = max, na.rm = TRUE) # observed occurrences as starting values for z
        )
        zst[!is.finite(zst)] <- 1 # performs replacement for values of -Inf (arising from sites with all NAs, i.e. extra Polygon_IDs without any replicates) as well as NA
        return(zst)
    })

    # species that had at least one detection (i.e. not the data augmentation species)
    w.start <- lapply(y, FUN = function (X) apply(X, MARGIN = 3, FUN = max, na.rm = TRUE) )

# extract and scale per-site occupancy covariates

    site.covariates.scaled <- sapply(names(leech.augmented), simplify = FALSE, function (X) {
        env.data = leech.augmented[[X]] %>%
            select(Polygon_ID, elevation_median, tpi_median, distance_to_road_median, distance_to_stream_median, distance_to_nature_reserve_boundary) %>%
            distinct()
        site.covs = tibble(Polygon_ID = y[[X]] %>% rownames()) %>% left_join(env.data, by = "Polygon_ID") # ensure row order is same as y
        site.covs.scaled <- scale(site.covs[,-1])
        rownames(site.covs.scaled) <- site.covs$Polygon_ID
        return(site.covs.scaled)
    })

# set b parameter to be supplied to JAGS for beta prior distribution on omega
# b is parameter for beta prior distribution on omega
# add b to model data depending on M from command line

    b.values <- list(`100` = 0.6, `150` = 3.3, `200` = 6.1)
    b <- b.values[[M]]

# save number of detected species to facilitate modelling only these species with occupancy covariates
#  - species 1:num.detected.species are the original detected species (e.g. first num.detected.species rows of y[[]])
#  - species (num.detected.species + 1):M are the augmented species
    num.detected.species <- sapply(names(leech.augmented), simplify = FALSE, function (X) sum(species.augment[[X]]$observed) )

# save prepared data to be packaged up for JAGS

    save(leech, leech.supplement, leech.augmented,
        num.sites, num.reps, num.species,
        numleeches, site.covariates.scaled, y,
        g, num.species.groups, species.groups, species.group.labels,
        z.start, w.start, b, num.detected.species,
        file = output.filename)
