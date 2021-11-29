# R code for generating community and detection occupancy predictions
# LSU dataset

library("here")
library("tidyverse")
library("jagsUI")

########################################################################################
# input and output file names

    filenames <- list(

        modeldata.all.rdata = here("rdata", "Ailaoshan_model_data_200.rdata"), # rdata file containing data for all models

        modeldata.LSU.rdata = here("rdata", "Ailaoshan_model_data_200_final_LSU.rdata"), # rdata file containing data for final LSU model

        modeloutput.rds = here("rds", "Ailaoshan_model_output_200_final_LSU.rds"), # file containing modelling results

        output.rdata = here("rdata", "Ailaoshan_occupancy_covariates_LSU_community_predictions.rdata") # file name to save predictions to

    )

########################################################################################
# get data

    # file containing LSU model data
    # (take polygon and OTU labels from jags.data$model.data$y or jags.data$model.data$z.start as stored in this rdata file)
    load(filenames$modeldata.LSU.rdata)

    # file containing model data prior to pulling out only the data required for the models
    # (use leech.augmented from this .rdata file for the unscaled covariate values)
    load(filenames$modeldata.all.rdata)

    # modelling output
    model.output <- readRDS(file = filenames$modeloutput.rds)

########################################################################################
# get number of chains and samples

    nchains <- length(model.output$samples) # each element of model.output$samples represents one mcmc chain

    nsamp <- nrow(model.output$samples[[1]]) # nrow should be the same for any of these chains

########################################################################################
# occupancy estimates

# elevation

    # predictor values
    # jags.data$model.data$occ[1,] is elev
    elev.vals <- data.frame(
        scaled = seq(min(jags.data$model.data$occ["elev",], na.rm = TRUE), max(jags.data$model.data$occ["elev",], na.rm = TRUE), length.out = 500),
        unscaled = seq(min(leech.augmented$LSU$elevation_median, na.rm = TRUE), max(leech.augmented$LSU$elevation_median, na.rm = TRUE), length.out = 500)
    )

    # generate community occupancy predictions based on the mcmc iterations
    elev.pred <- array(NA, dim = c(nrow(elev.vals), nchains * nsamp, 2)) # empty array to be filled
    for(c in 1:nchains){
        for(s in 1:nsamp){
            elev.pred[, s + (c - 1) * nsamp, 1] <- plogis(model.output$samples[[c]][s, "mu.eta[1,1]"] + model.output$samples[[c]][s, "mu.beta[1]"] * elev.vals$scaled)  # psi ~ elev for group 1 (i.e. mammals/birds)
            elev.pred[, s + (c - 1) * nsamp, 2] <- plogis(model.output$samples[[c]][s, "mu.eta[1,2]"] + model.output$samples[[c]][s, "mu.beta[1]"] * elev.vals$scaled)  # psi ~ elev for group 2 (i.e. amphibians/reptiles)
        }
    }

    #calculate means and 95% credible intervals
    elev.pred.mean <- apply(elev.pred, c(1,3), mean)
    elev.pred.BCI <- apply(elev.pred, c(1,3), function(x) quantile(x, prob = c(0.025, 0.975), na.rm=TRUE))
    rm(elev.pred)

    # turn occupancy predictions into dataframes
    occupancy.elevation <- tibble(
        `elevation (m)` = c(elev.vals$unscaled, elev.vals$unscaled), # note switch to unscaled values
        group = factor(c(rep("mammals/birds",nrow(elev.vals)), rep("amphibians/reptiles",nrow(elev.vals))),
            levels = c("mammals/birds","amphibians/reptiles")),
        mean = c(elev.pred.mean[,1], elev.pred.mean[,2]),
        `2.5%` = c(elev.pred.BCI[1,,1], elev.pred.BCI[1,,2]),
        `97.5%` = c(elev.pred.BCI[2,,1], elev.pred.BCI[2,,2])
    )

# reserve

    # predictor values
    # jags.data$model.data$occ[2,] is reserve
    reserve.vals <- data.frame(
        scaled = seq(min(jags.data$model.data$occ["reserve",], na.rm = TRUE), max(jags.data$model.data$occ["reserve",], na.rm = TRUE), length.out = 500),
        unscaled = seq(min(leech.augmented$LSU$distance_to_nature_reserve_boundary, na.rm = TRUE), max(leech.augmented$LSU$distance_to_nature_reserve_boundary, na.rm = TRUE), length.out = 500)
    )

    # generate community occupancy predictions based on the mcmc iterations
    reserve.pred <- array(NA, dim = c(nrow(reserve.vals), nchains * nsamp, 2)) # empty array to be filled
    for(c in 1:nchains){
        for(s in 1:nsamp){
            reserve.pred[, s + (c - 1) * nsamp, 1] <- plogis(model.output$samples[[c]][s, "mu.eta[1,1]"] + model.output$samples[[c]][s, "mu.beta[2]"] * reserve.vals$scaled)  # psi ~ reserve for group 1 (i.e. mammals/birds)
            reserve.pred[, s + (c - 1) * nsamp, 2] <- plogis(model.output$samples[[c]][s, "mu.eta[1,2]"] + model.output$samples[[c]][s, "mu.beta[2]"] * reserve.vals$scaled)  # psi ~ reserve for group 2 (i.e. amphibians/reptiles)
        }
    }

    #calculate means and 95% credible intervals
    reserve.pred.mean <- apply(reserve.pred, c(1,3), mean)
    reserve.pred.BCI <- apply(reserve.pred, c(1,3), function(x) quantile(x, prob = c(0.025, 0.975), na.rm=TRUE))
    rm(reserve.pred)

    # turn occupancy predictions into dataframes
    occupancy.reserve <- tibble(
        `distance to reserve edge (m)` = c(reserve.vals$unscaled, reserve.vals$unscaled),
        group = factor(c(rep("mammals/birds",nrow(reserve.vals)), rep("amphibians/reptiles",nrow(reserve.vals))),
            levels = c("mammals/birds","amphibians/reptiles")),
        mean = c(reserve.pred.mean[,1], reserve.pred.mean[,2]),
        `2.5%` = c(reserve.pred.BCI[1,,1], reserve.pred.BCI[1,,2]),
        `97.5%` = c(reserve.pred.BCI[2,,1], reserve.pred.BCI[2,,2])
    )


########################################################################################
# detection estimates

    # predictor values
    numleeches.vals <- seq(1,100, length.out = 500)

    # generate community detection predictions based on the mcmc iterations
    detection.pred <- array(NA, dim = c(length(numleeches.vals), nchains * nsamp, 2)) # empty array to be filled
    for(c in 1:nchains){
        for(s in 1:nsamp){
            detection.pred[, s + (c - 1) * nsamp, 1] <- 1 - ((1 - plogis(model.output$samples[[c]][s, "mu.eta[2,1]"])) ^ (numleeches.vals/100))     # p ~ numleeches for group 1 (i.e. mammals/birds)
            detection.pred[, s + (c - 1) * nsamp, 2] <- 1 - ((1 - plogis(model.output$samples[[c]][s, "mu.eta[2,2]"])) ^ (numleeches.vals/100))     # p ~ numleeches for group 2 (i.e. amphibians/reptiles)
        }
    }

    #calculate means and 95% credible intervals
    detection.pred.mean <- apply(detection.pred, c(1,3), mean)
    detection.pred.BCI <- apply(detection.pred, c(1,3), function(x) quantile(x, prob = c(0.025, 0.975), na.rm=TRUE))
    rm(detection.pred)

    # turn detection predictions into dataframe
    detection <- tibble(
        `number of leeches` = c(numleeches.vals, numleeches.vals),
        group = factor(c(rep("mammals/birds", length(numleeches.vals)), rep("amphibians/reptiles", length(numleeches.vals))),
                       levels = c("mammals/birds", "amphibians/reptiles")),
        mean = c(detection.pred.mean[,1], detection.pred.mean[,2]),
        `2.5%` = c(detection.pred.BCI[1,,1], detection.pred.BCI[1,,2]),
        `97.5%` = c(detection.pred.BCI[2,,1], detection.pred.BCI[2,,2])
    )

########################################################################################
# save predictions to file

LSU.community.predictions <- list(
    occupancy.elevation = occupancy.elevation,
    occupancy.reserve = occupancy.reserve,
    detection = detection
)

save(LSU.community.predictions, file = filenames$output.rdata)

########################################################################################
