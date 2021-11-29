# R code to package up data for initial occupancy models to gauge values for parameters

library("tidyverse")

# get command line arguments

    args <- commandArgs(trailingOnly = TRUE)
    input.filename <- args[1]
    output.filename <- list(
        LSU = args[2], # LSU jags data
        SSU = args[3]  # SSU jags data
    )

# get pre-prepared data

    load(file = input.filename)

# package data for JAGS

    model.data <- sapply(names(leech.augmented), simplify = FALSE, function (X) {
        list(y = y[[X]],
            num.sites = num.sites[[X]],
            num.reps = num.reps[[X]],
            num.species = num.species[[X]],
            # occupancy covariates
                occupancy.slopes = 5,  # number of occupancy slope coefficients
                occ = rbind(elev = site.covariates.scaled[[X]][,"elevation_median"],
                    tpi = site.covariates.scaled[[X]][,"tpi_median"],
                    road = site.covariates.scaled[[X]][,"distance_to_road_median"],
                    stream = site.covariates.scaled[[X]][,"distance_to_stream_median"],
                    reserve = site.covariates.scaled[[X]][,"distance_to_nature_reserve_boundary"]),
            # sampling covariate
                numleeches = numleeches[[X]],
            # species groups
                num.species.groups = num.species.groups[[X]], # number of different species groups, i.e. 2
                g = g[[X]], # group for each species (1 = mammals/birds, 2 = amphibians/squamates)
            # b parameter for beta prior distribution on omega
                b = b,
            # number of detected species
                num.detected.species = num.detected.species[[X]]
        )
    })

# initial values

    # function supplies z.start and w.start as present in jags.data at runtime
    # (using a function allows same values of z.start and w.start to be used for all chains; JAGS generates
    # random values for everything else, but this function can also be used to generate and provide random
    # values at runtime if desired)
    inits <- function () list(z = jags.data$z.start, w = jags.data$w.start)

# parameters to monitor

    params <- c("mu.beta", "sigma.beta", "beta.keep", "omega", "Ntotal")

# MCMC settings

    mcmc.settings <- list(
        ni = 60000, # total iterations (i.e. inclusive of burnin in both jagsUI and R2jags)
        nt = 1,     # no thinning
        nb = 30000, # burnin
        nc = 3      # number of chains
    )

# save data to files for modelling

    datasets <- c("LSU", "SSU")

    for (d in datasets) {
        jags.data <- list(
            model.data = model.data[[d]],
            z.start = z.start[[d]], w.start = w.start[[d]],
            inits = inits,
            params = params,
            mcmc.settings = mcmc.settings,
            species.groups = species.groups[[d]], species.group.labels = species.group.labels[[d]] # not used by JAGS
        )
        save(jags.data, file = output.filename[[d]])
    }
