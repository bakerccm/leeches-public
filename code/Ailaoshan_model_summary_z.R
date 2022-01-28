# R code for postprocessing occupancy modelling results (model with z saved)

library("tidyverse")
library("jagsUI")
library("coda")

########################################################################################
# get command line arguments

    args <- commandArgs(trailingOnly = TRUE)

    modeldata.rdata.filename <- args[1]    # file containing model data

    modeloutput.rds.filename <- args[2]    # file containing modelling results

    modelsummary.rds.filename <- args[3]    # file name to save results summary to

########################################################################################
# get model data
# (take polygon and OTU labels from jags.data$model.data$y or jags.data$model.data$z.start as stored in this rdata file)

    load(modeldata.rdata.filename)

# get data labels from jags.data

    jags.data.Polygon_ID <- dimnames(jags.data$model.data$y)[[1]] # Polygon_ID labels
    # identical(dimnames(jags.data$model.data$y)[[1]], rownames(jags.data$z.start)) # TRUE

    jags.data.replicate <- dimnames(jags.data$model.data$y)[[2]] # replicate labels

    jags.data.OTU <- dimnames(jags.data$model.data$y)[[3]] # OTU labels
    # identical(dimnames(jags.data$model.data$y)[[3]], colnames(jags.data$z.start)) # TRUE
    # identical(dimnames(jags.data$model.data$y)[[3]], names(jags.data$w.start)) # TRUE

    jags.data.occupancy.covariate <- rownames(jags.data$model.data$occ) # occupancy covariates

########################################################################################
# get modelling output

    model.output <- readRDS(file = modeloutput.rds.filename)

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

########################################################################################
# parse model estimates from JAGS output

# saved params:
# c("estimated.occupancy", "Nsite", "z", "beta0", "beta", "gamma0", "mu.eta", "mu.beta", "rho.beta0.gamma0")

# estimated.occupancy varies by species
# estimated.occupancy = proportion of sites occupied, for each species
    estocc.columns <- grep("^estimated.occupancy\\[", varnames(output.mcmc))
    estocc.output <- tibble(columns = estocc.columns) %>%
        mutate(OTU_ = varnames(output.mcmc)[columns])
    for (i in seq_along(estocc.output$columns)) {
        temp <- output.mcmc[,estocc.output$columns[i],drop=FALSE] %>% as.matrix()
        estocc.output[i,"mean"] <- mean(temp)
        estocc.output[i,"sd"] <- sd(temp)
        estocc.output[i,"2.5%"] <- quantile(temp, 0.025)
        estocc.output[i,"50%"] <- quantile(temp, 0.5)
        estocc.output[i,"97.5%"] <- quantile(temp, 0.975)
    }
    estocc.output <- estocc.output %>% select(-columns) %>%
        # reduce OTU_ to index number
            mutate(OTU_ = sub("estimated.occupancy\\[","",OTU_)) %>% mutate(OTU_ = sub("\\]","",OTU_)) %>%
            mutate(OTU_ = as.integer(OTU_)) %>%
        # convert OTU_ index numbers back to our original labels
            mutate(OTU = jags.data.OTU[OTU_]) %>%
            select(OTU, everything(), -OTU_)

# Nsite varies by site
# Nsite = number of species occurring at each site
    Nsite.columns <- grep("^Nsite\\[", varnames(output.mcmc))
    Nsite.output <- tibble(columns = Nsite.columns) %>%
        mutate(Polygon_ID_ = varnames(output.mcmc)[columns])
    for (i in seq_along(Nsite.output$columns)) {
        temp <- output.mcmc[,Nsite.output$columns[i],drop=FALSE] %>% as.matrix()
        Nsite.output[i,"mean"] <- mean(temp)
        Nsite.output[i,"sd"] <- sd(temp)
    }
    Nsite.output <- Nsite.output %>% select(-columns) %>%
        # reduce Polygon_ID_ to index number
            mutate(Polygon_ID_ = sub("Nsite\\[","",Polygon_ID_)) %>% mutate(Polygon_ID_ = sub("\\]","",Polygon_ID_)) %>%
            mutate(Polygon_ID_ = as.integer(Polygon_ID_)) %>%
        # convert Polygon_ID_ index numbers back to our original labels
            mutate(Polygon_ID = jags.data.Polygon_ID[Polygon_ID_]) %>%
            select(Polygon_ID, everything(), -Polygon_ID_)

# z varies by species and site
    z.columns <- grep("^z\\[", varnames(output.mcmc))
    z.output <- tibble(columns = z.columns) %>%
        mutate(var = varnames(output.mcmc)[columns])
    for (i in seq_along(z.output$columns)) {
        temp <- output.mcmc[,z.output$columns[i],drop=FALSE] %>% as.matrix()
        z.output[i,"mean"] <- mean(temp)
        z.output[i,"sd"] <- sd(temp)
    }
    z.output <- z.output %>% select(-columns) %>%
        # parse rownames into indexes for Polygon_ID and OTU
            mutate(var = sub("z\\[","",var)) %>% mutate(var = sub("\\]","",var)) %>%
            separate(var, sep=',', into=c("Polygon_ID_","OTU_")) %>%
            mutate(Polygon_ID_ = as.integer(Polygon_ID_), OTU_ = as.integer(OTU_)) %>%
        # convert Polygon_ID_ and OTU_ indexes back to our original labels
            mutate(Polygon_ID = jags.data.Polygon_ID[Polygon_ID_]) %>%
            mutate(OTU = jags.data.OTU[OTU_]) %>%
            select(Polygon_ID, OTU, everything(), -Polygon_ID_, -OTU_)
#    # other stats?
#    #     mutate(`2.5%`=as.integer(`2.5%`), `25%`=as.integer(`25%`), `50%`=as.integer(`50%`),
#    #     `75%`=as.integer(`75%`), `97.5%`=as.integer(`97.5%`), n.eff = as.integer(n.eff)) %>%

# beta0 varies by species
# beta0 is the species-specific constant for occupancy
# converting to probability scale should give occupancy estimate (=psi) for constant (i.e. mean) values
# for environmental covariates included in model, since predictors are scaled and centered
    beta0.columns <- grep("^beta0\\[", varnames(output.mcmc))
    beta0.output <- tibble(columns = beta0.columns) %>%
        mutate(OTU_ = varnames(output.mcmc)[columns])
    for (i in seq_along(beta0.output$columns)) {
        temp <- output.mcmc[,beta0.output$columns[i],drop=FALSE] %>% as.matrix()
        beta0.output[i,"mean"] <- mean(temp)
        beta0.output[i,"sd"] <- sd(temp)
        beta0.output[i,"2.5%"] <- quantile(temp, 0.025)
        beta0.output[i,"50%"] <- quantile(temp, 0.5)
        beta0.output[i,"97.5%"] <- quantile(temp, 0.975)
    }
    beta0.output <- beta0.output %>% select(-columns) %>%
        # reduce OTU_ to index number
            mutate(OTU_ = sub("beta0\\[","",OTU_)) %>% mutate(OTU_ = sub("\\]","",OTU_)) %>%
            mutate(OTU_ = as.integer(OTU_)) %>%
        # convert OTU_ index numbers back to our original labels
            mutate(OTU = jags.data.OTU[OTU_]) %>%
            select(OTU, everything(), -OTU_) %>%
    # add columns for beta0 values converted to probability scale (plogis = inverse logit)
            mutate(prob_mean = plogis(mean), `prob_2.5%` = plogis(`2.5%`), `prob_50%` = plogis(`50%`), `prob_97.5%` = plogis(`97.5%`)) %>%
            select(OTU, prob_mean, `prob_2.5%`, `prob_50%`, `prob_97.5%`, everything())

# beta varies by species
# beta is indexed as beta[x,k] where x={1,2,3,...} denotes which occupancy covariate and k denotes the species
    beta.columns <- grep("^beta\\[", varnames(output.mcmc))
    beta.output <- tibble(columns = beta.columns) %>%
        mutate(var = varnames(output.mcmc)[columns])
    for (i in seq_along(beta.output$columns)) {
        temp <- output.mcmc[,beta.output$columns[i],drop=FALSE] %>% as.matrix()
        beta.output[i,"mean"] <- mean(temp)
        beta.output[i,"sd"] <- sd(temp)
        beta.output[i,"2.5%"] <- quantile(temp, 0.025)
        beta.output[i,"25%"] <- quantile(temp, 0.25)
        beta.output[i,"50%"] <- quantile(temp, 0.50)
        beta.output[i,"75%"] <- quantile(temp, 0.75)
        beta.output[i,"97.5%"] <- quantile(temp, 0.975)
    }
    beta.output <- beta.output %>% select(-columns) %>%
        # parse rownames to indexes for occupancy.covariate and OTU
            mutate(var = sub("beta\\[","",var)) %>% mutate(var = sub("\\]","",var)) %>%
            separate(var, sep=',', into=c("occupancy.covariate_","OTU_")) %>%
            mutate(occupancy.covariate_ = as.integer(occupancy.covariate_), OTU_ = as.integer(OTU_)) %>%
        # convert occupancy.covariate_ and OTU_ index numbers back to our original labels
            mutate(occupancy.covariate = jags.data.occupancy.covariate[occupancy.covariate_]) %>%
            mutate(OTU = jags.data.OTU[OTU_]) %>%
            select(OTU, occupancy.covariate, everything(), -OTU_, -occupancy.covariate_)
    beta.output <- split(beta.output, beta.output$occupancy.covariate) # split into list of dataframes, one for each occupancy covariate

# gamma0 varies by species
# gamma0 is the species-specific constant for detection probability per 100 leeches on the logit scale
# converting to probability scale should give estimate of detection probability per 100 leeches (=r)
    gamma0.columns <- grep("^gamma0\\[", varnames(output.mcmc))
    gamma0.output <- tibble(columns = gamma0.columns) %>%
        mutate(OTU_ = varnames(output.mcmc)[columns])
    for (i in seq_along(gamma0.output$columns)) {
        temp <- output.mcmc[,gamma0.output$columns[i],drop=FALSE] %>% as.matrix()
        gamma0.output[i,"mean"] <- mean(temp)
        gamma0.output[i,"sd"] <- sd(temp)
        gamma0.output[i,"2.5%"] <- quantile(temp, 0.025)
        gamma0.output[i,"50%"] <- quantile(temp, 0.5)
        gamma0.output[i,"97.5%"] <- quantile(temp, 0.975)
    }
    gamma0.output <- gamma0.output %>% select(-columns) %>%
        # reduce OTU_ to index number
            mutate(OTU_ = sub("gamma0\\[","",OTU_)) %>% mutate(OTU_ = sub("\\]","",OTU_)) %>%
            mutate(OTU_ = as.integer(OTU_)) %>%
        # convert OTU_ index numbers back to our original labels
            mutate(OTU = jags.data.OTU[OTU_]) %>%
            select(OTU, everything(), -OTU_) %>%
    # add columns for gamma0 values converted to probability scale (plogis = inverse logit)
            mutate(prob_mean = plogis(mean), `prob_2.5%` = plogis(`2.5%`), `prob_50%` = plogis(`50%`), `prob_97.5%` = plogis(`97.5%`)) %>%
            select(OTU, prob_mean, `prob_2.5%`, `prob_50%`, `prob_97.5%`, everything())

# "mu.eta"
    mu.eta.columns <- grep("^mu.eta\\[", varnames(output.mcmc))
    mu.eta.output <- tibble(columns = mu.eta.columns) %>%
        mutate(var = varnames(output.mcmc)[columns])
    for (i in seq_along(mu.eta.output$columns)) {
        temp <- output.mcmc[,mu.eta.output$columns[i],drop=FALSE] %>% as.matrix()
        mu.eta.output[i,"mean"] <- mean(temp)
        mu.eta.output[i,"sd"] <- sd(temp)
    }
    mu.eta.output <- mu.eta.output %>% select(-columns) %>%
        # reduce rownames to index numbers for intercept (beta0 or gamma0) and taxonomic.group
            mutate(var = sub("mu.eta\\[","",var)) %>% mutate(var = sub("\\]","",var)) %>%
            separate(var, sep=',', into=c("intercept_","taxonomic.group_")) %>%
            mutate(intercept_ = as.integer(intercept_), taxonomic.group_ = as.integer(taxonomic.group_)) %>%
        # convert intercept_ and taxonomic.group_ index numbers back to our original labels
            mutate(intercept = ifelse(intercept_ == 1, "beta0", "gamma0"), taxonomic.group = jags.data$species.group.labels[taxonomic.group_]) %>%
            select(intercept, taxonomic.group, everything(), -intercept_, -taxonomic.group_)

# mu.beta varies by occupancy covariate
    mu.beta.columns <- grep("^mu.beta", varnames(output.mcmc)) # note no bracket at end
    mu.beta.output <- tibble(columns = mu.beta.columns) %>%
        mutate(occupancy.covariate_ = varnames(output.mcmc)[columns]) %>%
        # add index number if var name is just mu.beta, which occurs when there is only one occupancy covariate
            mutate(occupancy.covariate_ = ifelse(occupancy.covariate_ == "mu.beta", "mu.beta[1]", occupancy.covariate_))
    for (i in seq_along(mu.beta.output$columns)) {
        temp <- output.mcmc[,mu.beta.output$columns[i],drop=FALSE] %>% as.matrix()
        mu.beta.output[i,"mean"] <- mean(temp)
        mu.beta.output[i,"sd"] <- sd(temp)
    }
    mu.beta.output <- mu.beta.output %>% select(-columns) %>%
        # reduce occupancy.covariate_ to index number
            mutate(occupancy.covariate_ = sub("mu.beta\\[","",occupancy.covariate_)) %>% mutate(occupancy.covariate_ = sub("\\]","",occupancy.covariate_)) %>%
            mutate(occupancy.covariate_ = as.integer(occupancy.covariate_)) %>%
        # convert occupancy.covariate_ index numbers back to our original labels
            mutate(occupancy.covariate = jags.data.occupancy.covariate[occupancy.covariate_]) %>%
            select(occupancy.covariate, everything(), -occupancy.covariate_)

# rho.beta0.gamma0 varies by taxonomic group
    rho.columns <- grep("^rho.beta0.gamma0\\[", varnames(output.mcmc))
    rho.output <- tibble(columns = rho.columns) %>%
        mutate(taxonomic.group_ = varnames(output.mcmc)[columns])
    for (i in seq_along(rho.output$columns)) {
        temp <- output.mcmc[,rho.output$columns[i],drop=FALSE] %>% as.matrix()
        rho.output[i,"mean"] <- mean(temp)
        rho.output[i,"sd"] <- sd(temp)
    }
    rho.output <- rho.output %>% select(-columns) %>%
        # reduce taxonomic.group_ to index number
            mutate(taxonomic.group_ = sub("rho.beta0.gamma0\\[","",taxonomic.group_)) %>% mutate(taxonomic.group_ = sub("\\]","",taxonomic.group_)) %>%
            mutate(taxonomic.group_ = as.integer(taxonomic.group_)) %>%
        # convert taxonomic.group_ index numbers back to our original labels
            mutate(taxonomic.group = jags.data$species.group.labels[taxonomic.group_]) %>%
            select(taxonomic.group, everything(), -taxonomic.group_)

########################################################################################
# saves processed model output as a list
model.summary <- list(
    estocc.output = estocc.output,
    Nsite.output = Nsite.output,
    z.output = z.output,
    beta0.output = beta0.output,
    beta.output = beta.output,
    gamma0.output = gamma0.output,
    mu.eta.output = mu.eta.output,
    mu.beta.output = mu.beta.output,
    rho.output = rho.output
)
saveRDS(model.summary, file = modelsummary.rds.filename)

########################################################################################
