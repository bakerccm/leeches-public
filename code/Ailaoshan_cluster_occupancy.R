# R code for getting species-level occupancies for each cluster of sites

# This file take inputs <modeldata.rdata>, <modeloutput.rds> and <siteclusters.rds>
# and saves estimated occupancy for each species in each cluster as <output.rds>

library("tidyverse")
library("jagsUI")
library("coda")

########################################################################################
# get command line arguments

    args <- commandArgs(trailingOnly = TRUE)

    modeldata.rdata.filename <- args[1]  # file containing model data
    modeloutput.rds.filename <- args[2]  # file containing modelling results
    siteclusters.rds.filename <- args[3]  # file containing site clusters
    output.rds.filename <- args[4]  # file name to save site cluster occupancy to

########################################################################################
# model input data
# (take polygon and OTU labels from jags.data$model.data$y or jags.data$model.data$z.start as stored in this rdata file)

    load(modeldata.rdata.filename)

# modelling output

    model.output <- readRDS(file = modeloutput.rds.filename)

# site clusters

    site.clusters <- readRDS(file = siteclusters.rds.filename)
    site.clusters <- site.clusters %>% mutate(jaccard.cut = as.numeric(jaccard.cut))

    num.clusters <- site.clusters$jaccard.cut %>% unique() %>% length()

########################################################################################
# parse posterior samples into an array

    # convert samples to mcmc.list so coda can access them

        output.mcmc <- as.mcmc.list(model.output$samples)
        rm(model.output)

        # notes:
        # output.mcmc %>% str() # list of mcmc chains; rows are iterations and columns are variables
        # output.mcmc[[1]] # first mcmc chain
        # output.mcmc[1:10, 5, drop=FALSE] # first 10 iterations from the fifth variable
        # varnames(output.mcmc) == colnames(output.mcmc[[1]]) # variable names

        # can summarize mcmc.like this (and can customize quantiles):
        # summary(output.mcmc[, 5:10, drop=FALSE]) -> x

        # to convert a subset of the samples to matrix or array:
        # output.mcmc[,selected.columns,drop=FALSE] %>% as.matrix() # as.matrix() concatenates all the chains, giving an (iterations x variables) matrix
        # output.mcmc[,selected.columns,drop=FALSE] %>% as.array() # as.array() keeps chains separate, giving an (iterations x variables x chains) array

    # get just the z matrix from the posterior samples
    # rows are iterations, columns are all combinations of site x species

        z.columns <- colnames(output.mcmc[[1]]) %>% grep('^z\\[', .)
        z.matrix <- output.mcmc[, z.columns, drop=FALSE] %>% as.matrix()
        rm(z.columns, output.mcmc)

    # parse z[.,.] parameter names from the posterior z matrix

        z.colnames <- colnames(z.matrix) %>%
            gsub("^z\\[", "", .) %>% gsub("\\]$", "", .) %>%
            str_split(",", simplify = TRUE) %>%
            as.numeric() %>% matrix(ncol=2)
        colnames(z.colnames) <- c("site", "species")

    # determine dimensions of array to store z matrix posterior samples in

        num.iterations <- nrow(z.matrix) # note, this is the number of iterations saved, which is less than actual iterations if thinning is used
        num.sites <- max(z.colnames[,"site"])
        num.species <- max(z.colnames[,"species"])

    # convert z matrix to array indexed as iteractions, sites, species

        z.array <- array(z.matrix, dim=c(num.iterations, num.sites, num.species))
        rm(z.matrix)

    # add dimnames to z array

    dimnames(z.array) <- list(
        1:num.iterations,
        dimnames(jags.data$model.data$y)[[1]], # sites, i.e. PolygonIDs
        dimnames(jags.data$model.data$y)[[3]]  # species, i.e. OTU names
    )

########################################################################################
# which sites are in each cluster?

    cluster.PolygonIDs <- lapply(1:num.clusters, function (X) site.clusters$Polygon_ID[site.clusters$jaccard.cut == X])

# which columns in z.array correspond to these sites?

    cluster.indexes <- lapply(1:num.clusters, function (X) which(dimnames(z.array)[[2]] %in% cluster.PolygonIDs[[X]]))

# how many sites are in each cluster?

    cluster.sitecounts <- sapply(1:num.clusters, function (X) sum(site.clusters$jaccard.cut == X))

    cluster.sitecounts

########################################################################################
# for each MCMC iteration, get fraction of sites in each cluster occupied by each species

    cluster.occ <- array(NA, dim = c(num.iterations, num.clusters, num.species))

    dimnames(cluster.occ) <- list(1:num.iterations, 1:num.clusters, dimnames(z.array)[[3]])

    for (cluster in 1:num.clusters) {
        for (species in 1:num.species) {
            cluster.occ[,cluster,species] <- (z.array[,cluster.indexes[[cluster]],species] / cluster.sitecounts[cluster]) %>% apply(MARGIN = 1, sum)
        }
    }

    # inspect results
    # cluster.occ[1:10,,1:2]

########################################################################################
# summarize results and save to file

    # get means and quantiles
    cluster.occ.summary = list(
        mean = apply(cluster.occ, MARGIN = c(2,3), mean),
        BCI = apply(cluster.occ, MARGIN = c(2,3), function (X) quantile(X, probs = c(0.025, 0.25, 0.50, 0.75, 0.975)))
    )

    # dim(cluster.occ.summary$mean)
    # [1]  3 59

    # dim(cluster.occ.summary$BCI)
    # [1]  2  3 59

    # convert summary results to list of data.frames
    cluster.occ.summary.dflist <- list(
        mean = cluster.occ.summary$mean[,] %>%
            as.data.frame() %>%
            rownames_to_column(var = "cluster") %>%
            pivot_longer(!cluster, names_to = "OTU", values_to = "mean"),
        "2.5%" = cluster.occ.summary$BCI["2.5%",,] %>%
            as.data.frame() %>%
            rownames_to_column(var = "cluster") %>%
            pivot_longer(!cluster, names_to = "OTU", values_to = "2.5%"),
        "25%" = cluster.occ.summary$BCI["25%",,] %>%
            as.data.frame() %>%
            rownames_to_column(var = "cluster") %>%
            pivot_longer(!cluster, names_to = "OTU", values_to = "25%"),
        "50%" = cluster.occ.summary$BCI["50%",,] %>%
            as.data.frame() %>%
            rownames_to_column(var = "cluster") %>%
            pivot_longer(!cluster, names_to = "OTU", values_to = "50%"),
        "75%" = cluster.occ.summary$BCI["75%",,] %>%
            as.data.frame() %>%
            rownames_to_column(var = "cluster") %>%
            pivot_longer(!cluster, names_to = "OTU", values_to = "75%"),
        "97.5%" = cluster.occ.summary$BCI["97.5%",,] %>%
            as.data.frame() %>%
            rownames_to_column(var = "cluster") %>%
            pivot_longer(!cluster, names_to = "OTU", values_to = "97.5%")
    )

    # join data.frames ready for export
    cluster.occ.summary.df <- cluster.occ.summary.dflist[["mean"]] %>%
        full_join(cluster.occ.summary.dflist[["2.5%"]], by = c("cluster", "OTU")) %>%
        full_join(cluster.occ.summary.dflist[["25%"]], by = c("cluster", "OTU")) %>%
        full_join(cluster.occ.summary.dflist[["50%"]], by = c("cluster", "OTU")) %>%
        full_join(cluster.occ.summary.dflist[["75%"]], by = c("cluster", "OTU")) %>%
        full_join(cluster.occ.summary.dflist[["97.5%"]], by = c("cluster", "OTU")) %>%
        select(OTU, cluster, everything())

    # save output
    saveRDS(cluster.occ.summary.df, file = output.rds.filename)
