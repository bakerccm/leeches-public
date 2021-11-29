# R code to process the preOTU data in preparation for making network plots

library("here")
library("tidyverse")
library("readxl")

# get read tables and protax assignments

    preOTUtables.filepath <- here("data","OTUs_protax_outputs_12S_16S_weighted_unweighted_20190624.xlsx")

# get read tables and convert to tidy

    read.tbl <- list(
        LSU = read_excel(path = preOTUtables.filepath, sheet = "16S_otu_table_swarm_lulu_201906", na = "NA", col_names = TRUE),
        SSU = read_excel(path = preOTUtables.filepath, sheet = "12S_otu_table_swarm_lulu_201906", na = "NA", col_names = TRUE)
    )
    read.tidy <- lapply(read.tbl, FUN = function (X) gather(X, -OTU, key = "Lab_ID", value = "reads"))
    rm(read.tbl)

# get protax output tables

    taxon.tbl <- list(
        LSU = read_excel(path = preOTUtables.filepath, sheet = "protaxout_swarm_16S_weighted_un", na = "NA", col_names = TRUE),
        SSU = read_excel(path = preOTUtables.filepath, sheet = "protaxout_swarm_12S_weighted_un", na = "NA", col_names = TRUE)
    )

# get collection metadata

    collections.filepath <- here("data", "Ailaoshan_leech_2016_collections_20180406.xlsx")

    # note that collections$Polygon_ID still contains non-numeric fields like '25 or 24' that get converted to NA upon import and give rise to warnings here
    suppressWarnings(
        collections <- read_excel(path = collections.filepath, sheet = "All", col_names = TRUE, na=c("NA","？？？")) %>%
            # select(region_English, Ranger_ID, Polygon_ID, Lab_ID.aggregated, Big_tubes, Small_tubes, leech_qty, starts_with("Tube_"), lysis_buffer, leech_volume) %>%
            select(region_English, Ranger_ID, Polygon_ID, Lab_ID.aggregated) %>%
            arrange(Lab_ID.aggregated)
    )

# get lab metadata

    leech_qty.filepath <- here("data", "2016_AilaoshanLeeches_info.xlsx")

    leech_qty <- read_excel(path = leech_qty.filepath, sheet = "leech_qty", col_names = TRUE, na=c("NA","？？？")) %>%
        arrange(Lab_ID) %>%
        mutate(Lab_ID.aggregated = str_replace(Lab_ID, "\\.\\d*", "")) %>%
        mutate(Lab_ID.extraction = ifelse(Lab_ID.aggregated == Lab_ID, 1, str_replace(Lab_ID, ".*\\.", "")))

# join collections and lab metadata

    replicate.metadata <- leech_qty %>% left_join(collections, by = "Lab_ID.aggregated") %>% select(Lab_ID.aggregated, Lab_ID, Lab_ID.extraction, everything())

# compare Lab_IDs in read.tidy and replicate.metadata

    # LSU data

        # Lab_IDs that are in read.tidy but not in replicate.metadata
        read.tidy$LSU %>% filter(!(Lab_ID %in% replicate.metadata$Lab_ID)) %>% pull(Lab_ID) %>% unique()
        #  [1] "HUMAN1"  "HUMAN2"  "HUMAN3"  "HUMAN4"  "HUMAN5"  "HUMAN6"  "HUMAN7"  "HUMAN8"  "HUMAN10" "HUMAN12" "CXNCE"   "HUMAN9"  "HUMAN11" "HUMAN13" "JDNCE"   "NC11"    "SBNCD1"  "XPNCD1"  "ZYNCD"

        # Lab_IDs that are in replicate.metadata but not in read.tidy
        replicate.metadata %>% filter(!(Lab_ID %in% read.tidy$LSU$Lab_ID)) %>% pull(Lab_ID) %>% unique()
        #   [1] "CX09"    "CX27"    "CX30"    "CX31"    "CX32"    "CX33"    "CX34"    "CX35"    "CX36"    "CX37"    "CX38"    "CX39"    "JD02"    "JD04"    "JD06"    "JD07"    "JD08"    "JD09.1"  "JD09.2"  "JD09.3"  "JD09.4"  "JD09.5"
        #  [23] "JD12.1"  "JD12.4"  "JD12.5"  "JD15"    "JD16"    "JD20"    "JD23"    "JD24"    "JD25"    "JD28.1"  "JD28.3"  "JD28.5"  "JD29.3"  "JD29.5"  "JD30"    "JD35"    "JD36"    "JD42"    "JD44"    "JD45"    "JD53"    "JD56"
        #  [45] "JD57.1"  "JD57.4"  "JD58.1"  "JD58.2"  "JD58.4"  "JD58.5"  "JD59.1"  "JD59.2"  "JD59.4"  "JD59.5"  "JD60.1"  "JD60.2"  "JD60.3"  "JD60.5"  "JD61.1"  "JD61.2"  "JD61.3"  "JD61.4"  "JD61.5"  "JD64"    "JD66"    "JD71"
        #  [67] "JD73.1"  "JD73.3"  "JD73.4"  "JD75"    "JD79.1"  "JD79.2"  "JD79.3"  "JD79.4"  "JD79.5"  "JD82"    "JD86.1"  "JD86.2"  "JD86.3"  "JD86.4"  "JD89.1"  "JD89.2"  "JD89.3"  "JD89.4"  "JD89.5"  "JD91.2"  "JD91.3"  "JD91.4"
        #  [89] "JD92.3"  "JD92.4"  "JD92.5"  "JD96.3"  "NH03"    "NH04"    "NH07"    "NH08"    "NH09"    "NH10"    "NH14"    "NH17"    "NH21"    "NH22"    "NH23"    "NH24"    "NH25"    "NH26"    "NH28"    "NH29"    "SB007"   "SB008"
        # [111] "SB014.1" "SB014.3" "SB014.5" "SB017.2" "SB017.3" "SB020"   "SB024.1" "SB024.2" "SB024.4" "SB025.4" "SB029"   "SB032"   "SB034"   "SB039.3" "SB039.4" "SB039.5" "SB041"   "SB044"   "SB045"   "SB049"   "SB053"   "SB058.1"
        # [133] "SB060.3" "SB060.5" "SB072"   "SB080.1" "SB081.5" "SB082.4" "SB083"   "SB086"   "SB087"   "SB091"   "SB092"   "SB098.2" "SB098.4" "SB098.5" "SB100"   "SB106"   "SB108"   "SB110"   "SB114"   "SB118"   "SB119"   "SB120"
        # [155] "SB122"   "SB137"   "SB142"   "SB144"   "SB146"   "SB148"   "SB154"   "XP008"   "XP013.2" "XP015.2" "XP032"   "XP042"   "XP044.1" "XP047.3" "XP062"   "XP065.4" "XP092.5" "XP093.3" "XP094.3" "XP094.4" "XP095.1" "XP095.3"
        # [177] "XP095.5" "XP096.2" "XP096.3" "XP097.2" "XP097.3" "XP097.4" "XP097.5" "XP100"   "XP109.1" "XP109.2" "XP109.3" "XP109.4" "XP109.5" "XP111.1" "XP111.3" "XP111.5" "ZY02"    "ZY12.1"  "ZY12.2"  "ZY12.5"  "ZY14.3"  "ZY14.4"
        # [199] "ZY15.2"  "ZY16.2"  "ZY16.4"  "ZY17.5"  "ZY18.1"  "ZY19.2"  "ZY19.3"  "ZY19.4"  "ZY20.2"  "ZY20.4"  "ZY20.5"  "ZY21.2"  "ZY21.4"  "ZY21.5"  "ZY22.1"  "ZY22.2"  "ZY22.3"  "ZY24.2"  "ZY25"    "ZY29"    "ZY30"    "ZY31"
        # [221] "ZY32"    "ZY33"    "ZY38"    "ZY41"    "ZY42"    "ZY43"    "ZY45"    "ZY46"    "ZY47"    "ZY48"    "ZY49"    "ZY50"    "ZY53"    "ZY54.4"  "ZY54.5"  "ZY56"    "ZY57"    "ZY61"    "ZY62.3"

    # SSU data

        # Lab_IDs that are in read.tidy but not in replicate.metadata
        read.tidy$SSU %>% filter(!(Lab_ID %in% replicate.metadata$Lab_ID)) %>% pull(Lab_ID) %>% unique()
        #     [1] "HUMAN4"  "HUMAN14" "HUMAN10" "HUMAN13" "HUMAN5"  "HUMAN6"  "HUMAN7"  "HUMAN8"  "HUMAN9"  "HUMAN11" "HUMAN12" "HUMAN1"  "HUMAN2"  "HUMAN3"

        replicate.metadata %>% filter(!(Lab_ID %in% read.tidy$SSU$Lab_ID)) %>% pull(Lab_ID) %>% unique()
        #   [1] "CX07"    "CX09"    "CX10"    "CX18"    "CX27"    "CX30"    "CX31"    "CX32"    "CX33"    "CX34"    "CX35"    "CX36"    "CX37"    "CX38"    "CX39"    "JD06"    "JD07"    "JD08"    "JD09.1"  "JD09.2"  "JD09.3"  "JD09.4"
        #  [23] "JD09.5"  "JD12.1"  "JD12.5"  "JD15"    "JD16"    "JD28.1"  "JD28.2"  "JD28.3"  "JD28.4"  "JD29.5"  "JD36"    "JD38"    "JD42"    "JD44"    "JD45"    "JD57.1"  "JD57.5"  "JD58.2"  "JD58.4"  "JD58.5"  "JD59.2"  "JD59.4"
        #  [45] "JD59.5"  "JD60.5"  "JD61.1"  "JD61.2"  "JD61.3"  "JD61.4"  "JD61.5"  "JD64"    "JD66"    "JD73.1"  "JD73.3"  "JD75"    "JD79.5"  "JD89.2"  "JD89.4"  "JD89.5"  "JD91.3"  "JD91.4"  "JD92.5"  "NH03"    "NH04"    "NH06"
        #  [67] "NH07"    "NH08"    "NH09"    "NH10"    "NH14"    "NH21"    "NH22"    "NH23"    "NH24"    "NH25"    "NH26"    "NH28"    "NH29"    "SB007"   "SB008"   "SB014.3" "SB017.3" "SB024.1" "SB025.4" "SB039.3" "SB039.4" "SB060.3"
        #  [89] "SB060.5" "SB061.2" "SB072"   "SB078.2" "SB078.3" "SB078.4" "SB078.5" "SB081.5" "SB082.4" "SB083"   "SB086"   "SB091"   "SB092"   "SB098.2" "SB100"   "SB109"   "SB110"   "SB112"   "SB114"   "SB118"   "SB119"   "SB120"
        # [111] "SB122"   "SB133.2" "SB136.4" "SB137"   "SB139"   "SB146"   "SB154"   "XP026"   "XP028.5" "XP034.2" "XP034.3" "XP042"   "XP049"   "XP052"   "XP082"   "XP083"   "XP094.3" "XP100"   "XP109.5" "XP111.1" "XP111.3" "XP111.5"
        # [133] "XP120.2" "XP121.2" "XP123.2" "ZY02"    "ZY06"    "ZY07"    "ZY12.1"  "ZY20.2"  "ZY20.4"  "ZY20.5"  "ZY21.5"  "ZY22.3"  "ZY24.3"  "ZY27"    "ZY32"    "ZY40"    "ZY50"    "ZY56"    "ZY62.3"  "ZY62.5"

# remove control samples

    # Remove Lab_IDs that are not in collections metadata. The samples excluded are control samples (see code block above).

    for (i in names(read.tidy)) {
        read.tidy[[i]] <- read.tidy[[i]] %>%
            filter(Lab_ID %in% replicate.metadata$Lab_ID)
    }
    rm(i)

# add class to taxon.tbl

    mammal.orders <- c("Primates", "Artiodactyla", "Rodentia", "Carnivora", "Soricomorpha", "Erinaceomorpha")

    bird.orders <- c("Passeriformes", "Galliformes", "Charadriiformes")

    amphibian.orders <- c("Anura", "Caudata")

    reptile.orders <- "Squamata"

    for (i in names(taxon.tbl)) {
        taxon.tbl[[i]] <- taxon.tbl[[i]] %>% mutate(class.inferred = NA) %>%
            mutate(class.inferred = ifelse(order %in% mammal.orders, "Mammals", class.inferred)) %>%
            mutate(class.inferred = ifelse(order %in% bird.orders, "Birds", class.inferred)) %>%
            mutate(class.inferred = ifelse(order %in% amphibian.orders, "Amphibians", class.inferred)) %>%
            mutate(class.inferred = ifelse(order %in% reptile.orders, "Reptiles", class.inferred))
    }
    rm(i)

# export data

    objects.to.save <- c("replicate.metadata", "read.tidy", "taxon.tbl")

    save(list = objects.to.save, file = here("rdata","Ailaoshan_preOTUs.rdata"))
