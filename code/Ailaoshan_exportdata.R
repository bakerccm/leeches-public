# R code to export data for publication

library("here")
library("tidyverse")
library("writexl")

# load data

    load(file=here("rdata","Ailaoshan_OTU_table.rdata"))

    load(file=here("rdata","Ailaoshan_IUCNdata.rdata"))

    # model summaries generated with Ailaoshan_model_summary.R
    model.summary <- list(
        LSU = readRDS(file = here("rds", "Ailaoshan_model_summary_200_final_LSU.rds")),
        SSU = readRDS(file = here("rds", "Ailaoshan_model_summary_200_final_SSU.rds"))
    )
    estocc.output <- lapply(model.summary, FUN = function (X) as.data.frame(X$estocc.output))
    gamma0.output <- lapply(model.summary, FUN = function (X) as.data.frame(X$gamma0.output))
    rm(model.summary)

# relabel (non-avian) reptiles as squamates
    
    leech <- leech %>% mutate(consensus.class = ifelse(consensus.class == "Reptiles", "Squamates", consensus.class))

    leech.supplement <- leech.supplement %>% mutate(consensus.class = ifelse(consensus.class == "Reptiles", "Squamates", consensus.class))

    taxa <- taxa %>% mutate(class = ifelse(class == "Reptiles", "Squamates", class))

    taxa.iucn <- taxa.iucn %>% mutate(class = ifelse(class == "Reptiles", "Squamates", class))

# prepare occupancy and detection estimates

    occupancy.estimates <- bind_rows(LSU = estocc.output$LSU, SSU = estocc.output$SSU, .id = "dataset") %>%
        select(dataset, OTU, mean, `2.5%`, `97.5%`) %>%
        rename(occupancy_mean = mean, `occupancy_2.5%` = `2.5%`, `occupancy_97.5%` = `97.5%`)

    detection.estimates <- bind_rows(LSU = gamma0.output$LSU, SSU = gamma0.output$SSU, .id = "dataset") %>%
        select(dataset, OTU, prob_mean, `prob_2.5%`, `prob_97.5%`) %>%
        rename(detection_prob_mean = prob_mean, `detection_prob_2.5%` = `prob_2.5%`, `detection_prob_97.5%` = `prob_97.5%`)

    model.estimates <- full_join(occupancy.estimates, detection.estimates, by = c("dataset", "OTU"))

# add occupancy model estimates to original dataframes

    taxa.iucn.estimates <- taxa.iucn %>%
        # add occupancy model estimates
            left_join(model.estimates, by = c("dataset", "OTU")) %>%
        # add rank by occupancy estimate within each dataset
            group_by(dataset) %>%
            mutate(`estimated occupancy rank` = as.integer(rank(desc(occupancy_mean)))) %>%
            ungroup() %>%
            select (dataset, `estimated occupancy rank`, everything()) %>%
        # round estimates to 3 decimal places
            mutate(occupancy_mean = format(round(occupancy_mean,3), nsmall = 3),
                `occupancy_2.5%` = format(round(`occupancy_2.5%`,3), nsmall = 3),
                `occupancy_97.5%` = format(round(`occupancy_97.5%`,3), nsmall = 3)) %>%
            mutate(detection_prob_mean = format(round(detection_prob_mean,3), nsmall = 3),
                `detection_prob_2.5%` = format(round(`detection_prob_2.5%`,3), nsmall = 3),
                `detection_prob_97.5%` = format(round(`detection_prob_97.5%`,3), nsmall = 3)) %>%
        # concatenate 95% credible intervals to a single string
            mutate(`occupancy_95%CRI` = paste0(`occupancy_2.5%`, " - ", `occupancy_97.5%`),
                `detection_prob_95%CRI` = paste0(`detection_prob_2.5%`, " - ", `detection_prob_97.5%`)) %>%
            select(-`occupancy_2.5%`, -`occupancy_97.5%`, -`detection_prob_2.5%`, -`detection_prob_97.5%`) %>%
        # move occupancy model estimates to existing (empty) columns
            mutate(`estimated occupancy` = occupancy_mean, `estimated detection` = detection_prob_mean) %>%
            mutate(`estimated occupancy 95%` = `occupancy_95%CRI`, `estimated detection 95%` = `detection_prob_95%CRI`) %>%
            select(-occupancy_mean, -detection_prob_mean, -`occupancy_95%CRI`, -`detection_prob_95%CRI`)

# export top 20 most abundant OTUs in each dataset for insertion in LaTeX
# note that these need some manual editing after export
# e.g. add common names for species without IUCN common names ('NA'), and splitting into two sub-tables,
#     add '=' to tied ranks, check capitalisation on common names

    Table_1_2_data <- taxa.iucn.estimates %>%
        # only keep some fields
            select(dataset, `estimated occupancy rank`, consensus.short, genus, species, `IUCN common name`, `Chinese common name`,
                `IUCN category`, `estimated occupancy`, `estimated occupancy 95%`) %>%
            mutate(`Occupancy (95\\% BCI)` = paste0(`estimated occupancy`, " (", `estimated occupancy 95%`, ")")) %>%
        # arrange by estimated occupancy rank and filter to top 20 within each dataset
            arrange(dataset, `estimated occupancy rank`) %>%
            filter(`estimated occupancy rank` <= 20) %>%
            rename(Rank = `estimated occupancy rank`) %>%
        # make species names italics unless it's only identified to family or higher
        # this gives warnings but it's OK: substrings that are not numeric return NA, which is the desired behavior
            # find names that don't end with a number and make them italics
                mutate(consensus.short.italics = is.na(as.numeric(substr(consensus.short, nchar(consensus.short), nchar(consensus.short))))) %>%
                mutate(consensus.short = ifelse(consensus.short.italics, paste0("\\textit{", consensus.short, "}"), consensus.short)) %>%
            # find names that do end with a number but have genus name, and make only genus italics
                mutate(consensus.genus.italics = !consensus.short.italics & !is.na(genus)) %>%
                mutate(consensus.short = ifelse(consensus.genus.italics, paste0("\\textit{", genus, "} ", species), consensus.short)) %>%
            rename(`Scientific name` = consensus.short) %>%
        # join common names
            rename(`Common name` = `IUCN common name`) %>%
            mutate(`Common name` = paste0(`Common name`, " (", `Chinese common name`, ")")) %>%
        # make common name '--' if not identified to species
            mutate(`Common name` = ifelse(!consensus.short.italics, "--", `Common name`)) %>%
        # make IUCN category '--' if not available
            mutate(`IUCN category` = ifelse(is.na(`IUCN category`), "--", `IUCN category`)) %>%
        select(-consensus.short.italics, -consensus.genus.italics, -genus, -species, -`Chinese common name`, -`estimated occupancy`, -`estimated occupancy 95%`)
        
        Table_1_2_data %>% filter(dataset == "LSU") %>% select(-dataset) %>%
            write.table(here("tables","Table1_top20_OTUs_LSU_latex.txt"), quote = FALSE, sep = " & ", row.names = FALSE, eol = " \\\\\n")
        Table_1_2_data %>% filter(dataset == "SSU") %>% select(-dataset) %>%
            write.table(here("tables","Table2_top20_OTUs_SSU_latex.txt"), quote = FALSE, sep = " & ", row.names = FALSE, eol = " \\\\\n")
        
 # examine domesticated species and IUCN categories

    # domesticated animals
        taxa.iucn.estimates %>% filter(consensus.short %in% c("Bos taurus", "Capra hircus", "Ovis aries"))

    # IUCN categories
        taxa.iucn.estimates %>% pull(`IUCN category`) %>% unique() # [1] NA   "NT" "VU" "EN" "LC"

    # IUCN threatened mammals
        taxa.iucn.estimates %>% filter(`IUCN category` %in% c("EN", "VU")) %>% filter(class == "Mammals") %>% arrange(`IUCN scientific name`)

    # IUCN threatened amphibians
        taxa.iucn.estimates %>% filter(`IUCN category` %in% c("EN", "VU")) %>% filter(class != "Mammals") %>% arrange(`IUCN scientific name`)

    # IUCN threatened or near-threatened OTUs (some duplicates because of two different datasets)
        taxa.iucn %>% filter(`IUCN category` %in% c("NT","VU","EN")) %>% arrange(consensus.short)

    # IUCN threatened or near-threatened OTUS (datasets combined)
        taxa.iucn %>% filter(`IUCN category` %in% c("NT","VU","EN")) %>%
            select(class, consensus.short, `IUCN common name`, `IUCN category`) %>% distinct()

# export IUCN threatened or near threatened species for insertion in LaTeX
# note that these need some manual editing after export
# e.g. add common names for species without IUCN common names ('NA'), and check capitalisation on common names

    taxa.iucn.estimates %>%
        # only keep threatened and near-threatened species
            filter(`IUCN category` %in% c("NT","VU","EN")) %>%
        # arrange rows
            arrange(class, genus, species) %>%
        # only keep some fields
            select(dataset, class, consensus.short, `IUCN common name`, `Chinese common name`,
                `IUCN category`, `estimated occupancy`, `estimated occupancy 95%`) %>%
        # concatenate point and interval estimates for occupancy
            mutate(occupancy = paste0(`estimated occupancy`, " (", `estimated occupancy 95%`, ")")) %>%
            mutate(occupancy = ifelse(is.na(occupancy), '--', occupancy)) %>%
            select(-`estimated occupancy`, -`estimated occupancy 95%`) %>%
        # make scientific name italics
            mutate(consensus.short = paste0("\\textit{", consensus.short, "} ")) %>%
            rename(`Scientific name` = consensus.short) %>%
        # make IUCN common name '--' if not identified to species
            mutate(`IUCN common name` = ifelse(is.na(`IUCN common name`), '--', `IUCN common name`)) %>%
        # join common names
            rename(`Common name` = `IUCN common name`) %>%
            mutate(`Common name` = paste0(`Common name`, " (", `Chinese common name`, ")")) %>%
            select(-`Chinese common name`) %>%
        # pivot wider
            pivot_wider(id_cols=c(class, `Scientific name`, `Common name`, `IUCN category`),
                names_from = dataset, names_glue = "{dataset} {.value}", values_from = `occupancy`, values_fill = "--") %>%
        # capitalize heading
            rename(`Class` = class) %>%
        # save to file
            write.table(here("tables","Table3_IUCN_OTUs_latex.txt"), quote = FALSE, sep = " & ", row.names = FALSE, eol = " \\\\\n")

# write taxon data to xlsx file

    xlsx.export.taxa = list(
        "LSU taxa" = taxa.iucn.estimates %>% filter(dataset == "LSU") %>% select(-dataset) %>% rename(`OTU name` = consensus.short),
        "SSU taxa" = taxa.iucn.estimates %>% filter(dataset == "SSU") %>% select(-dataset) %>% rename(`OTU name` = consensus.short)
    )

    write_xlsx(xlsx.export.taxa,
        path = here("tables","Supplementary_Data_1.xlsx"),
        col_names = TRUE, format_headers = TRUE)

    rm(xlsx.export.taxa)

# prepare OTU tables and metadata for export

    LSU.OTU.table <- leech %>% filter(dataset == "LSU") %>%
        select(OTU, Lab_ID, reads) %>%
        arrange(OTU, Lab_ID) %>%
        pivot_wider(id_cols = "OTU", names_from = "Lab_ID", values_from = "reads")

    SSU.OTU.table <- leech %>% filter(dataset == "SSU") %>%
        select(OTU, Lab_ID, reads) %>%
        arrange(OTU, Lab_ID) %>%
        pivot_wider(id_cols = "OTU", names_from = "Lab_ID", values_from = "reads")

    LabID.metadata <- leech %>%
        select(Lab_ID, Polygon_ID, leech_qty) %>%
        distinct() %>%
        arrange(Lab_ID) %>%
        rename(`replicate ID` = Lab_ID, `patrol area ID` = Polygon_ID, `number of leeches` = leech_qty)

    PolygonID.metadata <- leech %>%
        select(Polygon_ID, region_English, Ranger_ID, longitude, latitude, elevation_median, tpi_median,
            distance_to_road_median, distance_to_stream_median, distance_to_nature_reserve_boundary) %>%
        distinct() %>%
        arrange(Polygon_ID) %>%
        rename(`patrol area ID` = Polygon_ID, `region name` = region_English, `ranger ID` = Ranger_ID,
            `median elevation (m)` = elevation_median, `median topograpic position index (TPI)` = tpi_median,
            `median distance to nearest road (m)` = distance_to_road_median, `median distance to nearest stream (m)` = distance_to_stream_median,
            `centroid distance to reserve boundary (m)` = distance_to_nature_reserve_boundary)

# export OTU tables and metadata as xlsx files

    xlsx.export.OTUtables = list(
        "LSU OTU table" = LSU.OTU.table,
        "SSU OTU table" = SSU.OTU.table,
        "replicate metadata" = LabID.metadata,
        "patrol area metadata" = PolygonID.metadata
    )

    write_xlsx(xlsx.export.OTUtables,
        path = here("tables","Supplementary_Data_6.xlsx"),
        col_names = TRUE, format_headers = TRUE)

    rm(xlsx.export.OTUtables)

# write data to datafile

    save(model.estimates, taxa.iucn.estimates, LSU.OTU.table, SSU.OTU.table, LabID.metadata, PolygonID.metadata, file=here("rdata","Ailaoshan_dataexport.rdata"))
