# R code for generating occupancy and detection predictions for individual species
# SSU dataset

library("here")
library("tidyverse")
library("jagsUI")
#library("coda")

########################################################################################
# input and output file names

    filenames <- list(

        modeldata.all.rdata = here("rdata", "Ailaoshan_model_data_200.rdata"), # rdata file containing data for all models

        modeldata.SSU.rdata = here("rdata", "Ailaoshan_model_data_200_final_SSU.rdata"), # rdata file containing data for final SSU model

        modeloutput.rds = here("rds", "Ailaoshan_model_output_200_final_SSU.rds"), # file containing SSU modelling results

        output.rdata = here("rdata", "Ailaoshan_occupancy_covariates_SSU_individual_predictions.rdata") # file name to save predictions to

    )

########################################################################################
# get data

    # file containing SSU model data
    # (take polygon and OTU labels from jags.data$model.data$y or jags.data$model.data$z.start as stored in this rdata file)
    load(filenames$modeldata.SSU.rdata)

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
# extract info from model data

    # num.species <- jags.data$model.data$num.species # don't use this: it includes the augmented species which do not need to be plotted
    num.detected.species <- jags.data$model.data$num.detected.species

    species.groups <- tibble(OTU = names(jags.data$model.data$g), group.no = jags.data$model.data$g) %>%
            mutate(group.name = recode(group.no, "1" = "mammals/birds", "2" = "amphibians/reptiles"))

    species.names <- leech.augmented$SSU %>% select(OTU,consensus.short, consensus.class) %>% distinct() %>% mutate(consensus.class = tolower(consensus.class))

########################################################################################
# range of unscaled and equivalent scaled values to generate predictions for

    pred.vals <- list(
        elev = data.frame(
            scaled = seq(min(jags.data$model.data$occ["elev",], na.rm = TRUE), max(jags.data$model.data$occ["elev",], na.rm = TRUE), length.out = 500),
            unscaled = seq(min(leech.augmented$SSU$elevation_median, na.rm = TRUE), max(leech.augmented$SSU$elevation_median, na.rm = TRUE), length.out = 500)),
        numleeches = seq(1, 100, length.out = 100)
    )

# lists of empty arrays to fill with predictions

    # note that z.start is arranged so that the detected species are the first 1:num.detected.species columns
    predictions <- list(
        elev = array(NA, dim = c(500, num.detected.species), dimnames = list(1:500, colnames(jags.data$z.start)[1:num.detected.species])),
        numleeches = array(NA, dim = c(100, num.detected.species), dimnames = list(1:100, colnames(jags.data$z.start)[1:num.detected.species]))
    )

# generate predictions
    for(k in 1:num.detected.species){
        # empty arrays to fill with predictions for species k
            elev.pred <- array(NA, dim = c(nrow(pred.vals$elev), nchains * nsamp)) # empty array to be filled
            numleeches.pred <- array(NA, dim = c(length(pred.vals$numleeches), nchains * nsamp)) # empty array to be filled
        # generate predictions based on the mcmc iterations
            for(c in 1:nchains){
                for(s in 1:nsamp){
                    elev.pred[, s + (c - 1) * nsamp] <- plogis(model.output$samples[[c]][s, paste0("beta0[", k, "]")] + model.output$samples[[c]][s, paste0("beta[1,", k, "]")] * pred.vals$elev$scaled) # psi ~ elev for species k
                    numleeches.pred[, s + (c - 1) * nsamp] <- 1 - ((1 - plogis(model.output$samples[[c]][s, paste0("gamma0[", k, "]")])) ^ (pred.vals$numleeches/100))     # p ~ numleeches for species k
                }
            }
        #calculate means and store in list of prediction arrays
            predictions$elev[,k] <- apply(elev.pred, 1, mean)
            predictions$numleeches[,k] <- apply(numleeches.pred, 1, mean)
        rm(elev.pred, numleeches.pred)
    }

# convert arrays to dataframes

    SSU.individual.predictions <- list(
        elev = predictions$elev %>% as.data.frame() %>%
            mutate(elev.scaled = pred.vals$elev$scaled, elev.unscaled = pred.vals$elev$unscaled) %>%
            pivot_longer(-c(elev.scaled, elev.unscaled), names_to = "OTU", values_to="estimated occupancy"),
        numleeches = predictions$numleeches %>% as.data.frame() %>%
            mutate(leeches = pred.vals$numleeches) %>%
            pivot_longer(-leeches, names_to = "OTU", values_to="estimated detection")
    )

    SSU.individual.predictions <- lapply(SSU.individual.predictions, function (X) {
        X %>% left_join(species.groups %>% select(OTU, group.name), by = "OTU") %>%
            left_join(species.names, by = "OTU")
    })

########################################################################################
# save predictions to file

    save(SSU.individual.predictions, file = filenames$output.rdata)

########################################################################################
