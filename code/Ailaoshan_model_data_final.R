# R code to package up data for final occupancy models

library("tidyverse")

# get command line arguments

    args <- commandArgs(trailingOnly = TRUE)
    input.filename <- args[1]
    output.filename <- list(
        LSU = args[2], # LSU jags data without thinning (z not saved)
        LSU.z = args[3], # LSU jags data with thinning
        SSU = args[4]  # SSU jags data without thinning (z not saved)
        SSU.z = args[5]  # SSU jags data with thinning
    )

# get pre-prepared data

    load(file = input.filename)

# package data for JAGS

    model.data <- sapply(names(leech.augmented), simplify = FALSE, function (X) {
        list(y = y[[X]],
            num.sites = num.sites[[X]],
            num.reps = num.reps[[X]],
            num.species = num.species[[X]],
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

    # add environmental covariates
    # note these differ between datasets

        model.data$LSU$occupancy.slopes <- 2 # number of occupancy slope coefficients
        model.data$LSU$occ <- rbind(elev = site.covariates.scaled[["LSU"]][,"elevation_median"],
                                    reserve = site.covariates.scaled[["LSU"]][,"distance_to_nature_reserve_boundary"])

        model.data$SSU$occupancy.slopes <- 1 # number of occupancy slope coefficients
        model.data$SSU$occ <- rbind(elev = site.covariates.scaled[["SSU"]][,"elevation_median"])

# initial values

    # function supplies z.start and w.start as present in jags.data at runtime
    # (using a function allows same values of z.start and w.start to be used for all chains; JAGS generates
    # random values for everything else, but this function can also be used to generate and provide random
    # values at runtime if desired)
    inits <- function () list(z = jags.data$z.start, w = jags.data$w.start)

# parameters to monitor

    # for version without thinning (z not saved)
        params <- c("estimated.occupancy", "Nsite", "beta0", "beta", "gamma0", "mu.eta", "mu.beta", "rho.beta0.gamma0", "Ntotal")
    # for version with thinning (z saved)
        params.z <- c("estimated.occupancy", "Nsite", "beta0", "beta", "gamma0", "mu.eta", "mu.beta", "rho.beta0.gamma0", "Ntotal", "z")

# MCMC settings

    # for version without thinning (z not saved)
        mcmc.settings <- list(
            ni = 100000, # total iterations (i.e. inclusive of burnin in both jagsUI and R2jags)
            nt = 1,     # no thinning
            nb = 50000, # burnin
            nc = 5      # number of chains
        )
    # for version with thinning (z saved)
        mcmc.settings.z <- list(
            ni = 100000, # total iterations (i.e. inclusive of burnin in both jagsUI and R2jags)
            nt = 10,     # 10-fold thinning
            nb = 50000, # burnin
            nc = 5      # number of chains
        )

# save data to files for modelling

    datasets <- c("LSU", "SSU")

    # for versions without thinning (z not saved)
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

    # for versions with thinning (z saved)
        for (d in datasets) {
            jags.data <- list(
                model.data = model.data[[d]],
                z.start = z.start[[d]], w.start = w.start[[d]],
                inits = inits,
                params = params.z,
                mcmc.settings = mcmc.settings.z,
                species.groups = species.groups[[d]], species.group.labels = species.group.labels[[d]] # not used by JAGS
            )
            save(jags.data, file = output.filename[[paste0(d,".z")]])
        }
