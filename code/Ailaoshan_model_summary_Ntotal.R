# R code for postprocessing occupancy modelling results
#  - this code gets Ntotal results from all 6 model runs (3 values of M for each dataset)

library("here")
library("tidyverse")
library("jagsUI")
library("coda")

########################################################################################
# get command line arguments

    # files containing modelling results
    input.filenames <- list(
        LSU = list(
            `100` = here("rds", "Ailaoshan_model_output_100_final_LSU.rds"),
            `150` = here("rds", "Ailaoshan_model_output_150_final_LSU.rds"),
            `200` = here("rds", "Ailaoshan_model_output_200_final_LSU.rds")
        ),
        SSU = list(
            `100` = here("rds", "Ailaoshan_model_output_100_final_SSU.rds"),
            `150` = here("rds", "Ailaoshan_model_output_150_final_SSU.rds"),
            `200` = here("rds", "Ailaoshan_model_output_200_final_SSU.rds")
        )
    )

    output.rds.filename <- here("rds", "Ailaoshan_model_summary_Ntotal.rds") # file name to save Ntotal results summary to

########################################################################################
# set up empty tibble for summary stats

    Ntotal.summary <- tibble(
        dataset = rep(c("LSU","SSU"), each=3),
        M = rep(c(100,150,200), times=2),
        mean = NA,
        sd = NA,
        min = NA,
        `2.5%` = NA,
        `25%` = NA,
        `50%` = NA,
        `75%` = NA,
        `97.5%` = NA,
        max = NA
    )

########################################################################################
# loop through results files one at a time to limit memory use

    for (d in c("LSU","SSU")) {
        for (m in c("100","150","200")) {

            # get modelling output

                model.output <- readRDS(file = input.filenames[[d]][[m]])

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

            # parse model estimates from JAGS output

                # there is just one Ntotal for each model run
                    Ntotal.columns <- grep("^Ntotal", varnames(output.mcmc))
                    Ntotal.mtx <- output.mcmc[,Ntotal.columns,drop=FALSE] %>% as.matrix()

                # which row of Ntotal.summary should results go in?
                    Ntotal.summary.row <- which(Ntotal.summary$dataset == d & Ntotal.summary$M == as.numeric(m))

                # extract summary stats and save to Ntotal.summary
                    Ntotal.summary[Ntotal.summary.row, "mean"] <- mean(Ntotal.mtx)
                    Ntotal.summary[Ntotal.summary.row, "sd"] <- sd(Ntotal.mtx)
                    Ntotal.summary[Ntotal.summary.row, "min"] <- min(Ntotal.mtx)
                    Ntotal.summary[Ntotal.summary.row, "2.5%"] <- quantile(Ntotal.mtx, 0.025)
                    Ntotal.summary[Ntotal.summary.row, "25%"] <- quantile(Ntotal.mtx, 0.25)
                    Ntotal.summary[Ntotal.summary.row, "50%"] <- quantile(Ntotal.mtx, 0.50)
                    Ntotal.summary[Ntotal.summary.row, "75%"] <- quantile(Ntotal.mtx, 0.75)
                    Ntotal.summary[Ntotal.summary.row, "97.5%"] <- quantile(Ntotal.mtx, 0.975)
                    Ntotal.summary[Ntotal.summary.row, "max"] <- max(Ntotal.mtx)

                # clean up
                    rm(Ntotal.mtx)

        }
    }

########################################################################################
# saves processed model output to RDS file
saveRDS(Ntotal.summary, file = output.rds.filename)

########################################################################################
