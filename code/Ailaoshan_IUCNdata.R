# R code to get IUCN data for species in dataset

# usage: Rscript Ailaoshan_IUCNdata.R <APIkey>
# use rl_use_iucn() to get API key first if you don't have one

library("here")
library("tidyverse")
library("rredlist")

# get data

    load(file=here("rdata","Ailaoshan_OTU_table.rdata"))

# get taxon lists

    taxa <- leech %>%
        # simplify naming
            rename(class = consensus.class, order = consensus.order, family = consensus.family, genus = consensus.genus, species = consensus.species, `Chinese common name` = Chinese_common_name) %>%
        # calculate detection (=true/false) by Polygon_ID
            group_by(dataset, OTU, consensus.short, class, order, family, genus, species, `Chinese common name`, Polygon_ID) %>%
            summarize(detected = (sum(reads) > 0), .groups = "drop_last") %>%
        # calculate fraction of sites occupied
            summarize(`observed occupancy` = sum(detected) / n(), .groups = "drop") %>%
            mutate(`observed occupancy` = `observed occupancy` %>% round(3) %>% format(nsmall=3)) %>%
        # reorder rows so that sp1, sp2 etc go after any named taxa
            mutate(species = ifelse(grepl("(sp\\d)", species), paste0("zzz ", species), species)) %>%
            arrange(dataset, class, order, family, genus, species) %>%
            mutate(species = sub("zzz ", "", species)) %>%
        # add some empty columns to be filled with occupancy model results later
            mutate(`estimated occupancy` = NA, `estimated occupancy 95%` = NA, `estimated detection` = NA, `estimated detection 95%` = NA) %>%
        # prepare names to search for in IUCN database
            mutate(`IUCN scientific name` = paste(genus, species)) %>%
            mutate(`IUCN scientific name` = ifelse(grepl("[0-9]$", species), NA, `IUCN scientific name`)) %>%
            ### search for Capricornis milneedwardsii under the name Capricornis sumatraensis ###
            mutate(`IUCN scientific name` = ifelse(`IUCN scientific name` == "Capricornis milneedwardsii", "Capricornis sumatraensis", `IUCN scientific name`)) %>%
            select(dataset, OTU, `IUCN scientific name`, everything())

# prepare to get IUCN data

    iucn.data <- taxa %>%
        # consolidate names to search by removing NAs and also any duplicates across datasets
            select(`IUCN scientific name`) %>% filter(!is.na(`IUCN scientific name`)) %>% distinct() %>%
        # set up empty fields to populate
            mutate(
                `IUCN common name` = NA,
                `IUCN category` = NA,
                `IUCN criteria` = NA,
                `IUCN population trend` = NA,
                `IUCN assessment date` = NA,
            ) %>%
        # save as data.frame, not tibble
            as.data.frame()

# get IUCN data
# this is slow, mostly due to the requirement that queries take place one at a time with a couple of seconds between.

    for (i in 1:nrow(iucn.data)) {
        curr_sp_data <- rl_search(iucn.data$`IUCN scientific name`[i], key = iucn.key)$result
        if (!(is_empty(curr_sp_data))) {
            iucn.data[i, "IUCN category"] <- curr_sp_data$category
            iucn.data[i, "IUCN population trend"] <- curr_sp_data$population_trend
            iucn.data[i, "IUCN assessment date"] <- curr_sp_data$assessment_date
            iucn.data[i, "IUCN criteria"] <- curr_sp_data$criteria
            iucn.data[i, "IUCN common name"] <- curr_sp_data$main_common_name
        }
        rm(curr_sp_data)
        Sys.sleep(3) # try increasing this value to slow down requests if they are getting rejected by the server
    }

# merge IUCN data back in with taxon data

    taxa.iucn <- taxa %>%
        left_join(iucn.data, by = "IUCN scientific name") %>%
        select(-starts_with("IUCN "), starts_with("IUCN "))

# save to file

    save(taxa, iucn.data, taxa.iucn, file=here("rdata","Ailaoshan_IUCNdata.rdata"))
