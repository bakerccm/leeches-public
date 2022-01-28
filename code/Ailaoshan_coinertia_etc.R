# R code to look at coinertia between LSU and SSU datasets

library("here")
library("ade4") # for co-inertia
library("tidyverse")

# get data

    load(file=here("rdata","Ailaoshan_OTU_table.rdata"))

    # model summaries generated with Ailaoshan_model_summary.R
    model.summary <- list(
        LSU = readRDS(file = here("rds", "Ailaoshan_model_summary_200_final_LSU_z.rds")),
        SSU = readRDS(file = here("rds", "Ailaoshan_model_summary_200_final_SSU_z.rds"))
    )
    z.output <- lapply(model.summary, FUN = function (X) as.data.frame(X$z.output))
    rm(model.summary)

# prepare data

    # mean z estimates from model (i.e. per-polygon z estimates)
    # species are columns, sites are rows
    z <- lapply(z.output, function (X) {
        X %>% select(mean, Polygon_ID, OTU) %>%
            pivot_wider(id_cols = Polygon_ID, names_from = OTU, values_from = mean) %>%
            select(-Polygon_ID) %>% as.data.frame()
    })

# co-inertia analysis between sites

    # prepare PCAs as input for co-inertia analysis
    # (N.B. coinertia function below retrieves the standardized data from the dudi.pca objects. The PCA results themselves are not used.)
    #   sites.dudi.pca <- lapply(z, function (X) dudi.pca(X, scale = FALSE, scan = FALSE))

    # co-inertia analysis
    #   coinertia.plot <- coinertia(sites.dudi.pca$LSU, sites.dudi.pca$SSU, scan=FALSE, nf=3)

    # plot the results
    #    pdf(here("figures", "coinertia.pdf"), width=12, height = 12)
    #    plot(coinertia.plot)
    #    dev.off()

# test coinertia

    RV.rtest(df1 = z$LSU, df2 = z$SSU, nrepet = 999)

    # Monte-Carlo test
    # Call: RV.rtest(df1 = z$LSU, df2 = z$SSU, nrepet = 999)
    #
    # Observation: 0.769695
    #
    # Based on 999 replicates
    # Simulated p-value: 0.001
    # Alternative hypothesis: greater
    #
    # Std.Obs  Expectation     Variance
    # 1.128773e+02 1.531611e-02 4.466486e-05
