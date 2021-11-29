# R code to run occupancy model with all variables and with Bayesian variable selection

library("jagsUI")

# get command line arguments

    args <- commandArgs(trailingOnly = TRUE)

    data.filename <- args[1]    # file containing jags.data
    jags.filename <- args[2]    # file containing jags model
    output.filename <- args[3]  # file to save model output to as data stream

# read in data

    load(file = data.filename)

########################################################################################
# get start time
starttime <- Sys.time()
print(paste("Starting:", starttime))

########################################################################################
# run model using jagsUI

set.seed(123)

model.output <- jags(
    model.file = jags.filename,
    data = jags.data$model.data,
    inits = jags.data$inits,
    parameters.to.save = jags.data$params,
    n.chains = jags.data$mcmc.settings$nc,
    n.adapt = NULL, # default NULL runs adaptive phase until JAGS determines model is adapted
    n.iter = jags.data$mcmc.settings$ni, # total iterations including burnin (but excluding adaptive phase)
    n.burnin = jags.data$mcmc.settings$nb,
    n.thin = jags.data$mcmc.settings$nt,
    parallel = TRUE,
    n.cores = NULL, # should we supply this? e.g. = jags.data$mcmc.settings$nc
    DIC = FALSE,
    verbose = FALSE)
    # other options:
    # store.data=FALSE, codaOnly=FALSE, seed=NULL, bugs.format=FALSE, verbose=TRUE

# save model results as data stream
saveRDS(model.output, file = output.filename)

########################################################################################
# produce summaries and trace plots

    # summary info on model run
    # (view code with jagsUI:::summary.jagsUI to see how the summary is derived)
    summary(model.output)

    # summaries of saved parameters
    model.output$summary

    # trace and density plots
    # plot(model.output)
    # densityplot(model.output)

########################################################################################
# get end time
endtime <- Sys.time()
print(paste("Finished:", endtime))
print(paste("Elapsed:", endtime - starttime))

########################################################################################
