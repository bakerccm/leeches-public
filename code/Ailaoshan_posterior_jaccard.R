# R code for calculating posterior mean jaccard matrices from occupancy model output

library("tidyverse")
library("jagsUI")
library("coda")

########################################################################################
# get command line arguments

    args <- commandArgs(trailingOnly = TRUE)

    filenames <- list(
        model.data = args[1],    # .rdata file containing model input data
        model.output = args[2],    # .rds file containing model output
        jaccard.similarities = args[3]    # .rds file to store jaccard similarities between sites
    )

########################################################################################
# get data

    # model input data
    # (take polygon and OTU labels from jags.data$model.data$y or jags.data$model.data$z.start as stored in this rdata file)

        load(filenames$model.data)

    # modelling output

        model.output <- readRDS(file = filenames$model.output)

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

        # sites, i.e. PolygonIDs
        dimnames(z.array)[[2]] <- dimnames(jags.data$model.data$y)[[1]]

        # species, i.e. OTU names
        dimnames(z.array)[[3]] <- dimnames(jags.data$model.data$y)[[3]]


########################################################################################
# jaccard calculations

    # calculate total species per site to save having to calculate these totals more than once below

        z.sum <- apply(z.array, MARGIN = c(1,2), FUN = sum) # returns matrix with dimensions interations x sites

    # function to calculate posterior mean jaccard similarities between two sites i and j
        pairwise.jaccard <- function (site.i, site.j) {
            shared.species = apply(z.array[,site.i,] * z.array[,site.j,], MAR = 1, sum) # vector of length num.iterations
            total.species = z.sum[,site.i] + z.sum[,site.j] - shared.species # vector of length num.iterations
            jaccard = shared.species / total.species # vector of length num.iterations
            mean(jaccard) # scalar
            # quantiles or other posterior statistics could also be calculated here
        }

# calculate posterior mean jaccard similarities pairwise for sites

    # get PolygonIDs to calculate Jaccard distances for
    # i.e. only polygons with information on location, i.e. excludes those that just represent a ranger without location information

        all.PolygonIDs <- dimnames(z.array)[[2]]
        real.PolygonID.indexes <- which(!is.na(as.numeric(all.PolygonIDs))) # this gives a warning about NAs being introduced by coercion, but this is the desired behavior
        names(real.PolygonID.indexes) <- all.PolygonIDs[real.PolygonID.indexes]

    # array of NAs to store jaccard similarities in

        jaccard.posterior.means <- array(NA, dim = c(length(real.PolygonID.indexes), length(real.PolygonID.indexes)), dimnames = list(names(real.PolygonID.indexes), names(real.PolygonID.indexes)))

    # loop over each pair of sites

        cat ("Calculating Jaccard distances\n")
        cat ("Progress:\n")
        for (i in seq_along(real.PolygonID.indexes)) {
            if ((i %% 5) == 0) cat(paste0("Index i = ", i, "\n"))
            for (j in seq_len(i-1)) {
                jaccard.posterior.means[i,j] <- pairwise.jaccard(real.PolygonID.indexes[i], real.PolygonID.indexes[j])
            }
            jaccard.posterior.means[i,i] <- 1
        }

    # duplicate calculated values since the distance matrix is symmetric

        jaccard.posterior.means[upper.tri(jaccard.posterior.means)] <- t(jaccard.posterior.means)[upper.tri(jaccard.posterior.means)]

########################################################################################
# output results

    saveRDS(jaccard.posterior.means, file = filenames$jaccard.similarities)
