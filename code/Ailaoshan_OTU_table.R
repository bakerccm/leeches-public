# R code to generate final OTU tables from preOTU data
#  -- preOTUs with the same final taxonomy assigned manually by DWY are collapsed to give final OTUs

library("here")
library("readxl")
library("tidyverse")

# input filenames

    # this is the preOTU file augmented with the final taxon assignment from Doug
    preOTUtables.filepath <- here("data","final taxonomy assignments","OTUs_protax_outputs_12S_16S_weighted_unweighted_20190624_20190830_20190924.xlsx")

    # this is the collections metadata
    collections.filepath <- here("data","Ailaoshan_leech_2016_collections_20180406.xlsx")

    # this is the lab metadata (number of leeches per extraction)
    leech_qty.filepath <- here("data","2016_AilaoshanLeeches_info.xlsx")

    # environmental metadata from Ailaoshan_environmental.Rmd
    environmental.filepath <- here("rdata","Ailaoshan_environmental.rdata")

    # pantheria mammals data filepath
    pantheria.filepath <- here("data", "pantheria", "PanTHERIA_1-0_WR05_Aug2008.txt")

# get read tables and protax assignments

    # read tables

        # Note: we use SSU and LSU to denote the (eukaryotic mitochondrial) 12S and 16S datasets respectively, since names like data$12S are not allowed in R, but names like data$ssu are fine.
        read.tbl <- list(
            LSU = read_excel(path = preOTUtables.filepath, sheet = "16S_otu_table_swarm_lulu_201906", na = "NA", col_names = TRUE),
            SSU = read_excel(path = preOTUtables.filepath, sheet = "12S_otu_table_swarm_lulu_201906", na = "NA", col_names = TRUE)
        )
        read.tbl <- lapply(read.tbl, FUN = function (X) X %>% rename(queryID = OTU)) # rename for consistency with other input data (these are preOTUs at this stage, not OTUs)
        # convert to tidy
        read.tidy <- lapply(read.tbl, FUN = function (X) gather(X, -queryID, key = "Lab_ID", value = "reads") %>% mutate(reads = as.integer(reads)))

    # total reads

        lapply(read.tidy, FUN = function (X) X %>% select(reads) %>% sum()) # LSU: 82653202, SSU: 91132342

    # protax output tables

        taxon.tbl <- list(
            LSU = read_excel(path = preOTUtables.filepath, sheet = "protaxout_swarm_16S_weighted_un", na = "NA", col_names = TRUE),
            SSU = read_excel(path = preOTUtables.filepath, sheet = "protaxout_swarm_12S_weighted_un", na = "NA", col_names = TRUE)
        )
        taxon.tbl$LSU <- taxon.tbl$LSU %>% mutate(dataset = "LSU") %>% select(dataset, everything())
        taxon.tbl$SSU <- taxon.tbl$SSU %>% mutate(dataset = "SSU") %>% select(dataset, everything())

# get collections and lab metadata

    ## note that collections$Polygon_ID still contains non-numeric fields like '25 or 24' that get converted to NA upon import and give rise to warnings here
    suppressWarnings(
        collections <- read_excel(path = collections.filepath, sheet = "All", col_names = TRUE, na=c("NA","？？？")) %>%
            select(region_English, Ranger_ID, Polygon_ID, Lab_ID.aggregated) %>%
            arrange(Lab_ID.aggregated)
    )

    leech_qty <- read_excel(path = leech_qty.filepath, sheet = "leech_qty", col_names = TRUE, na=c("NA","？？？")) %>%
        arrange(Lab_ID) %>%
        # extracted aggregated Lab_IDs from the more detailed ones in order to match to collections
        mutate(Lab_ID.aggregated = str_replace(Lab_ID, "\\.\\d*", "")) %>%
        mutate(leech_qty = as.integer(leech_qty))

    leech_qty %>% summarize(total.leeches = sum(leech_qty), .groups = "drop") # 30468 leeches total

    # join collections and lab metadata
    replicates <- leech_qty %>% left_join(collections, by = "Lab_ID.aggregated") %>% select(Lab_ID, everything())

    replicates %>% summarize(total.leeches = sum(leech_qty), .groups = "drop") # still have 30468 leeches after left_join

    # how many Lab_ID.aggregated's had (leech_qty > 100)? (these bags were supposed to be split up in lab processing)
    replicates %>%
        group_by(Lab_ID.aggregated) %>% summarize(leech_qty = sum(leech_qty), .groups = "drop") %>%
        mutate(`leech_qty > 100` = ifelse(leech_qty > 100, TRUE, FALSE)) %>%
        group_by(`leech_qty > 100`) %>% summarize(count = n(), .groups = "drop") %>% mutate(fraction = count / sum(count)) # fraction of bags, not collecting areas

    # how many Polygon_ID's had at least one Lab_ID.aggregated (i.e. bag) with (leech_qty > 100)?
    replicates %>%
        group_by(Lab_ID.aggregated, Polygon_ID) %>% summarize(leech_qty = sum(leech_qty), .groups = "drop") %>%
        group_by(Polygon_ID) %>% summarize(`at least one bag leech_qty > 100` = ifelse(any(leech_qty > 100), TRUE, FALSE), .groups = "drop") %>%
        group_by(`at least one bag leech_qty > 100`) %>% summarize(count = n(), .groups = "drop") %>% mutate(fraction = count / sum(count)) # fraction of bags , bot collect ing earas

    # discard Lab_ID.aggregated as we don't need it any further
    replicates <- replicates %>% select(-Lab_ID.aggregated)

# compare Lab_IDs in read.tidy and replicates metadata

    # LSU data

        # Lab_IDs that are in read.tidy but not in replicates metadata
        read.tidy$LSU %>% filter(!(Lab_ID %in% replicates$Lab_ID)) %>% pull(Lab_ID) %>% unique()
        #  [1] "HUMAN1"  "HUMAN2"  "HUMAN3"  "HUMAN4"  "HUMAN5"  "HUMAN6"  "HUMAN7"  "HUMAN8"  "HUMAN10" "HUMAN12" "CXNCE"   "HUMAN9"  "HUMAN11" "HUMAN13" "JDNCE"   "NC11"    "SBNCD1"  "XPNCD1"  "ZYNCD"

        # Lab_IDs that are in replicates metadata but not in read.tidy
        replicates %>% filter(!(Lab_ID %in% read.tidy$LSU$Lab_ID)) %>% pull(Lab_ID) %>% unique()
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

        # Lab_IDs that are in read.tidy but not in replicates metadata
        read.tidy$SSU %>% filter(!(Lab_ID %in% replicates$Lab_ID)) %>% pull(Lab_ID) %>% unique()
        #     [1] "HUMAN4"  "HUMAN14" "HUMAN10" "HUMAN13" "HUMAN5"  "HUMAN6"  "HUMAN7"  "HUMAN8"  "HUMAN9"  "HUMAN11" "HUMAN12" "HUMAN1"  "HUMAN2"  "HUMAN3"

        replicates %>% filter(!(Lab_ID %in% read.tidy$SSU$Lab_ID)) %>% pull(Lab_ID) %>% unique()
        #   [1] "CX07"    "CX09"    "CX10"    "CX18"    "CX27"    "CX30"    "CX31"    "CX32"    "CX33"    "CX34"    "CX35"    "CX36"    "CX37"    "CX38"    "CX39"    "JD06"    "JD07"    "JD08"    "JD09.1"  "JD09.2"  "JD09.3"  "JD09.4"
        #  [23] "JD09.5"  "JD12.1"  "JD12.5"  "JD15"    "JD16"    "JD28.1"  "JD28.2"  "JD28.3"  "JD28.4"  "JD29.5"  "JD36"    "JD38"    "JD42"    "JD44"    "JD45"    "JD57.1"  "JD57.5"  "JD58.2"  "JD58.4"  "JD58.5"  "JD59.2"  "JD59.4"
        #  [45] "JD59.5"  "JD60.5"  "JD61.1"  "JD61.2"  "JD61.3"  "JD61.4"  "JD61.5"  "JD64"    "JD66"    "JD73.1"  "JD73.3"  "JD75"    "JD79.5"  "JD89.2"  "JD89.4"  "JD89.5"  "JD91.3"  "JD91.4"  "JD92.5"  "NH03"    "NH04"    "NH06"
        #  [67] "NH07"    "NH08"    "NH09"    "NH10"    "NH14"    "NH21"    "NH22"    "NH23"    "NH24"    "NH25"    "NH26"    "NH28"    "NH29"    "SB007"   "SB008"   "SB014.3" "SB017.3" "SB024.1" "SB025.4" "SB039.3" "SB039.4" "SB060.3"
        #  [89] "SB060.5" "SB061.2" "SB072"   "SB078.2" "SB078.3" "SB078.4" "SB078.5" "SB081.5" "SB082.4" "SB083"   "SB086"   "SB091"   "SB092"   "SB098.2" "SB100"   "SB109"   "SB110"   "SB112"   "SB114"   "SB118"   "SB119"   "SB120"
        # [111] "SB122"   "SB133.2" "SB136.4" "SB137"   "SB139"   "SB146"   "SB154"   "XP026"   "XP028.5" "XP034.2" "XP034.3" "XP042"   "XP049"   "XP052"   "XP082"   "XP083"   "XP094.3" "XP100"   "XP109.5" "XP111.1" "XP111.3" "XP111.5"
        # [133] "XP120.2" "XP121.2" "XP123.2" "ZY02"    "ZY06"    "ZY07"    "ZY12.1"  "ZY20.2"  "ZY20.4"  "ZY20.5"  "ZY21.5"  "ZY22.3"  "ZY24.3"  "ZY27"    "ZY32"    "ZY40"    "ZY50"    "ZY56"    "ZY62.3"  "ZY62.5"

# remove control samples
#  -- removes Lab_IDs that are not in replicates metadata; the samples excluded are control samples (see code block above)

    for (i in names(read.tidy)) {
        read.tidy[[i]] <- read.tidy[[i]] %>%
            filter(Lab_ID %in% (replicates %>% pull(Lab_ID)))
    }
    rm(i)

# reads after removing control samples

    lapply(read.tidy, FUN = function (X) X %>% select(reads) %>% sum()) # LSU: 77772174, SSU: 90182071

# clean up protax output tables

    # LSU: 71 preOTUs without consensus taxonomy, 98 preOTUs with consensus taxonomy
    #    taxon.tbl$LSU %>% mutate(consensus.taxonomy.available = is.na(consensus_taxonomy_v3_KIZ)) %>% group_by(consensus.taxonomy.available) %>% tally()
    # SSU: 78 preOTUs without consensus taxonomy, 92 preOTUs with consensus taxonomy
    #    taxon.tbl$SSU %>% mutate(consensus.taxonomy.available = is.na(consensus_taxonomy_v3_KIZ)) %>% group_by(consensus.taxonomy.available) %>% tally()

    taxon.parsed <- list()

    taxon.parsed <- lapply(taxon.tbl, function (X) {
        X %>%
        # use only weighted rows
            filter(protaxmod == "weighted") %>% select(-protaxmod) %>%
        # remove rows without taxonomic info
            filter(!is.na(consensus_taxonomy_v3_KIZ)) %>%
        # split consensus_taxonomy_v3 string
            rename(consensus = consensus_taxonomy_v3_KIZ) %>%
            mutate(consensus.class = NA) %>% # fill in later
            separate(consensus, into = c("consensus.order", "consensus.family", "consensus.genus", "consensus.species"), remove = FALSE, fill = "right") %>%
        # make short consensus name
            mutate(consensus.short = paste0(consensus.genus, " ", consensus.species)) %>%
            mutate(consensus.short = ifelse(is.na(consensus.species), paste0(consensus.family, " ", consensus.genus), consensus.short)) %>%
            mutate(consensus.short = ifelse(is.na(consensus.genus), paste0(consensus.order, " ", consensus.family), consensus.short)) %>%
        # parse split consensus names
            mutate(consensus.genus = ifelse(is.na(consensus.genus), consensus.family, consensus.genus)) %>%
            mutate(consensus.species = ifelse(is.na(consensus.species), consensus.genus, consensus.species)) %>%
            mutate(consensus.genus = ifelse(consensus.genus == consensus.species & grepl("^sp\\d", consensus.genus), NA, consensus.genus)) %>%
            mutate(consensus.family = ifelse(consensus.family == consensus.species & grepl("^sp\\d", consensus.family), NA, consensus.family))
    })

    # fill in consensus.class

        (orders.in.dataset <- c(taxon.parsed$LSU %>% pull(consensus.order), taxon.parsed$SSU %>% pull(consensus.order)) %>% unique())

        # [1] "Primates"        "Anura"           "Artiodactyla"    "Caudata"         "Rodentia"        "Carnivora"       "Soricomorpha"    "Passeriformes"   "Galliformes"
        #[10] "Squamata"        "Erinaceomorpha"  "Charadriiformes"

        mammal.orders <- c("Primates", "Artiodactyla", "Rodentia", "Carnivora", "Soricomorpha", "Erinaceomorpha")
        bird.orders <- c("Passeriformes", "Galliformes", "Charadriiformes")
        amphibian.orders <- c("Anura", "Caudata")
        reptile.orders <- "Squamata"

        # check we have all the orders in the dataset
        c(mammal.orders, bird.orders, amphibian.orders, reptile.orders) %in% orders.in.dataset
        orders.in.dataset %in% c(mammal.orders, bird.orders, amphibian.orders, reptile.orders)

        for (i in names(taxon.tbl)) {
            taxon.parsed[[i]] <- taxon.parsed[[i]] %>%
                mutate(consensus.class = ifelse(consensus.order %in% mammal.orders, "Mammals", consensus.class)) %>%
                mutate(consensus.class = ifelse(consensus.order %in% bird.orders, "Birds", consensus.class)) %>%
                mutate(consensus.class = ifelse(consensus.order %in% amphibian.orders, "Amphibians", consensus.class)) %>%
                mutate(consensus.class = ifelse(consensus.order %in% reptile.orders, "Reptiles", consensus.class))
        }

        rm(i, mammal.orders, bird.orders, amphibian.orders, reptile.orders, orders.in.dataset)

    # create unique taxon labels

        # first calculate total reads per preOTU
        total_reads <- list()
        for (i in names(read.tidy)) {
            total_reads[[i]] <- read.tidy[[i]] %>%
                group_by(queryID) %>%
                summarize(total.reads = sum(reads), .groups="drop")
        }

        # now use total read counts to grab queryID number from the most abundant preOTU for each taxonomy
        for (i in names(taxon.parsed)) {
            taxon.parsed[[i]] <- taxon.parsed[[i]] %>%
                # get numeric bit of queryID
                    mutate(queryID.num = str_replace(queryID, "OTU", "") %>% as.numeric()) %>%
                # attach total reads
                    left_join(total_reads[[i]], by = "queryID") %>%
                # get queryID for preOTU with most reads within each unique taxonomic assignment
                    group_by(dataset, consensus) %>%
                    mutate(OTU = ifelse(total.reads == max(total.reads), queryID.num, -999)) %>%
                    mutate(OTU = max(OTU)) %>%
                    mutate(OTU = paste0(dataset, sprintf("%03d", OTU))) %>%
                    ungroup() %>%
                # rearrange tables
                    select(dataset, OTU, queryID, queryID.num, total.reads, consensus, consensus.short,
                           consensus.class, consensus.order, consensus.family, consensus.genus, consensus.species, everything()) %>%
                    arrange(OTU) %>%
                # remove intermediate fields
                    select(-queryID.num, -total.reads) %>%
                # remove old protax fields
                    select(-class, -prob_class,
                           -order, -prob_order,
                           -family, -prob_family,
                           -genus, -prob_genus,
                           -species, -prob_species,
                           -bestHit_similarity, -bestHit,
                           -evidence, -consensus_taxonomy_v3, -compare_v3_v2, -consensus_taxonomy_v2, -consensus_taxonomy_v1,
                           -BLAST, -pct_ID)
        }

        rm(i, total_reads)

# attach taxonomic info to tidy read data and merge preOTUs to make OTUs

    # attach taxonomic info and split preOTUs according to whether taxonomic info is available
    preOTUs <- list()
    preOTUs.matched <- list()
    preOTUs.unmatched <- list()
    for (i in names(read.tidy)) {
        # attach taxonomic info to preOTUs
            preOTUs[[i]] <- read.tidy[[i]] %>% left_join(taxon.parsed[[i]], by = "queryID")
        # split preOTU data (unmatched data will be discarded from analysis)
            preOTUs.matched[[i]] <- preOTUs[[i]] %>% filter(!is.na(OTU))
            preOTUs.unmatched[[i]] <- preOTUs[[i]] %>% filter(is.na(OTU))
        # Doug says Muntiacus muntjak is actually Muntiacus vaginalis due to some upstream error, so fixing it here:
            preOTUs.matched[[i]]$consensus[preOTUs.matched[[i]]$consensus == "Artiodactyla,Cervidae,Muntiacus,muntjak"] <- "Artiodactyla,Cervidae,Muntiacus,vaginalis"
            preOTUs.matched[[i]]$consensus.short[preOTUs.matched[[i]]$consensus.short == "Muntiacus muntjak"] <- "Muntiacus vaginalis"
            preOTUs.matched[[i]]$consensus.species[preOTUs.matched[[i]]$consensus.species == "muntjak"] <- "vaginalis"
    }
    rm(i, preOTUs)

    # collapse matched preOTUs to generate our final OTUs
    OTUs <- list()
    for (i in names(read.tidy)) {
        OTUs[[i]] <- preOTUs.matched[[i]] %>%
            group_by_at(vars(-queryID, -reads)) %>%
            summarize(reads = sum(reads), .groups="drop") %>%
            select(Lab_ID, dataset, OTU, reads, everything())
    }

    # How many reads were assigned taxonomic information (matched) or not (unmatched)?
    # N.B. unmatched reads, i.e. those without taxonomic information, are discarded from analysis.
        unmatched <- bind_rows(LSU= preOTUs.unmatched$LSU %>% select(reads), SSU = preOTUs.unmatched$SSU %>% select(reads), .id = "dataset")
        matched <- bind_rows(LSU= OTUs$LSU %>% select(reads), SSU = OTUs$SSU %>% select(reads), .id = "dataset")
        bind_rows(unmatched = unmatched, matched = matched, .id = "matched") %>%
            group_by(dataset, matched) %>%
            summarize(reads = sum(reads), .groups="drop") %>%
            mutate(percent = round(reads/sum(reads) * 100, 1))
        rm(matched, unmatched)

# generate data object for export

    # concatenate reads data from LSU and SSU dataset
        leech <- bind_rows(OTUs$LSU, OTUs$SSU)

    # attach selected replicates metadata

        leech <- leech %>% left_join(replicates, by = "Lab_ID")

    # remove Lab_IDs that have no reads (these should be Lab_IDs that only had unmatched OTUs)
    # -- otherwise the calculated fraction of human reads for those samples is NaN
    # Note that samples with only human reads are still in here, but will get reduced to zero reads when we remove Homo sapiens later.

        # leech %>% group_by(dataset, Lab_ID) %>%
        #    summarize(total.reads = sum(reads), .groups="drop") %>%
        #    filter(total.reads == 0)  # LSU Lab_ID ZY22.4 has no reads and will be excluded in the line below

        leech <- leech %>% group_by(dataset, Lab_ID) %>% filter(sum(reads) > 0) %>% ungroup()

    # remove samples without Ranger_ID or Polygon_ID information
    # which samples are going to be removed?
        # leech %>% filter(is.na(Polygon_ID) & is.na(Ranger_ID))  # one sample is missing both: Lab_ID = ZY47
    # actually remove reads
        leech <- leech %>% filter(!is.na(Polygon_ID) | !is.na(Ranger_ID))
    # how many reads now?
        # leech %>% group_by(dataset) %>% summarize(reads = sum(reads), .groups="drop") # LSU: 77092321, SSU: 89943206

    # use Ranger_ID in place of missing Polygon_IDs
        leech <- leech %>%
            add_column(Polygon_ID_original = NA, .after = "Polygon_ID") %>%
            mutate(Polygon_ID_original = Polygon_ID) %>%
            mutate(Polygon_ID = ifelse(is.na(Polygon_ID), Ranger_ID, Polygon_ID))

    # add id numbers -- arbitrary but unique within each Polygon_ID (which we are going to treat
    # as a unique 'location') -- to facilitate converting data for occupancy modelling
    # (Note: replicate IDs will not necessarily match across the datasets, which is probably
    # appropriate since replicates may have worked for one dataset but not the other)
        replicate.numbers <- leech %>%
            select(dataset, Polygon_ID, Lab_ID) %>% distinct() %>%
            group_by(dataset, Polygon_ID) %>% mutate(replicate_no = seq_along(Lab_ID)) %>% ungroup() %>%
            mutate(replicate_no = paste0("replicate_",sprintf("%02d", replicate_no)))

        leech <- leech %>% left_join(replicate.numbers, by=c("dataset", "Polygon_ID", "Lab_ID")) %>%
            select(dataset, Polygon_ID, Lab_ID, replicate_no, OTU, everything())

    # add fraction of reads from humans in each Lab_ID

        leech <- leech %>%
            group_by(dataset, Lab_ID) %>%
            mutate(fraction.reads.humans = sum(ifelse(consensus.short == "Homo sapiens", reads, 0) / sum(reads))) %>%
            ungroup()

# get environmental metadata

    load(file = environmental.filepath)

# augment leech with extra rows for Polygon_IDs with missing data

    leech.data.to.duplicate <- list()
    Polygon_IDs.to.add <- list()
    leech.extrarows <- list()

    for (i in c("LSU", "SSU")) {

        leech.data.to.duplicate[[i]] <- leech %>% filter(dataset == i) %>%
            # get all the taxon-specific columns
                select(dataset, OTU, consensus, consensus.short, consensus.class, consensus.order,
                   consensus.family, consensus.genus, consensus.species, Chinese_common_name) %>% distinct() %>%
            # add back in columns that will need to be filled
                mutate(replicate_no = "replicate_01", reads = NA, leech_qty = NA, Ranger_ID = NA, fraction.reads.humans = NA, region_English = NA) %>%
            # save as data frame
                as.data.frame()

        Polygon_IDs.to.add[[i]] <- env.data %>%
            filter(!(Polygon_ID %in% (leech %>% filter(dataset == i) %>% pull(Polygon_ID)))) %>%
            select(Polygon_ID) %>% as.data.frame()

        leech.extrarows[[i]] <- bind_cols(
                as_tibble(Polygon_IDs.to.add[[i]][rep(seq_along(1:nrow(Polygon_IDs.to.add[[i]])), each = nrow(leech.data.to.duplicate[[i]])), , drop = FALSE]),
                as_tibble(leech.data.to.duplicate[[i]][rep(seq_along(1:nrow(leech.data.to.duplicate[[i]])), times = nrow(Polygon_IDs.to.add[[i]])),])
            ) %>% mutate(Lab_ID = paste0(".", Polygon_ID), Polygon_ID_original = Polygon_ID) %>%
            select(names(leech))

    }

    leech <- leech %>%
        bind_rows(leech.extrarows$LSU) %>%
        bind_rows(leech.extrarows$SSU) %>%
        arrange(dataset, Polygon_ID, Lab_ID, OTU)

    rm(leech.data.to.duplicate, Polygon_IDs.to.add, leech.extrarows)

# attach environmental metadata

    leech <- leech %>% left_join(env.data, by="Polygon_ID")

# add data on domestic vs non-domestic and size to leech

    # load pantheria database
    pantheria <- read_tsv(file = pantheria.filepath)

    # make table of mammal taxa to add domestic/wild and body size info to
    mammals <- leech %>% filter(consensus.class == "Mammals") %>%
            select(consensus.order, consensus.family, consensus.genus, consensus.short) %>% distinct() %>%
            arrange(consensus.order, consensus.family, consensus.genus) %>%
            mutate(domestic = NA)

    # add domestic/non-domestic info
    mammals <- mammals %>%
        mutate(domestic = ifelse(consensus.short %in% c("Bos taurus", "Capra hircus", "Ovis aries"), "domestic", "non-domestic"))

    # add body weight for species that match
    mammals.pantheria <- mammals %>%
        left_join(pantheria %>% select(MSW05_Binomial, `5-1_AdultBodyMass_g`), by=c("consensus.short" = "MSW05_Binomial")) %>%
        rename(AdultBodyMass_g = `5-1_AdultBodyMass_g`) %>%
        mutate(AdultBodyMass_g = ifelse(AdultBodyMass_g == -999, NA, AdultBodyMass_g)) %>%
        mutate(AdultBodyMass_g_source = ifelse(is.na(AdultBodyMass_g), NA, "species"))

    # add genus body weight means where species are missing
    pantheria.genus.summary <- pantheria %>% select(MSW05_Genus, `5-1_AdultBodyMass_g`) %>%
        mutate(`5-1_AdultBodyMass_g` = ifelse(`5-1_AdultBodyMass_g` == -999, NA, `5-1_AdultBodyMass_g`)) %>%
        group_by(MSW05_Genus) %>%
        summarize(`5-1_AdultBodyMass_g.mean` = mean(`5-1_AdultBodyMass_g`, na.rm = TRUE),
            `5-1_AdultBodyMass_g.median` = median(`5-1_AdultBodyMass_g`, na.rm = TRUE), .groups="drop")
    mammals.pantheria <- mammals.pantheria %>%
        left_join(pantheria.genus.summary %>% select(MSW05_Genus, `5-1_AdultBodyMass_g.mean`), by=c("consensus.genus" = "MSW05_Genus")) %>%
        mutate(AdultBodyMass_g_source = ifelse(!is.na(AdultBodyMass_g_source), AdultBodyMass_g_source, ifelse(!is.na(`5-1_AdultBodyMass_g.mean`), "genus mean", NA))) %>%
        mutate(AdultBodyMass_g = ifelse(!is.na(AdultBodyMass_g), AdultBodyMass_g, ifelse(!is.na(`5-1_AdultBodyMass_g.mean`), `5-1_AdultBodyMass_g.mean`, NA))) %>%
        select(-`5-1_AdultBodyMass_g.mean`)

    # add family body weight means where species and genus are missing
    pantheria.family.summary <- pantheria %>% select(MSW05_Family, `5-1_AdultBodyMass_g`) %>%
        mutate(`5-1_AdultBodyMass_g` = ifelse(`5-1_AdultBodyMass_g` == -999, NA, `5-1_AdultBodyMass_g`)) %>%
        group_by(MSW05_Family) %>%
        summarize(`5-1_AdultBodyMass_g.mean` = mean(`5-1_AdultBodyMass_g`, na.rm = TRUE),
            `5-1_AdultBodyMass_g.median` = median(`5-1_AdultBodyMass_g`, na.rm = TRUE), .groups="drop")
    mammals.pantheria <- mammals.pantheria %>%
        left_join(pantheria.family.summary %>% select(MSW05_Family, `5-1_AdultBodyMass_g.mean`), by=c("consensus.family" = "MSW05_Family")) %>%
        mutate(AdultBodyMass_g_source = ifelse(!is.na(AdultBodyMass_g_source), AdultBodyMass_g_source, ifelse(!is.na(`5-1_AdultBodyMass_g.mean`), "family mean", NA))) %>%
        mutate(AdultBodyMass_g = ifelse(!is.na(AdultBodyMass_g), AdultBodyMass_g, ifelse(!is.na(`5-1_AdultBodyMass_g.mean`), `5-1_AdultBodyMass_g.mean`, NA))) %>%
        select(-`5-1_AdultBodyMass_g.mean`)

    # add family body weight means where species and genus are missing
    pantheria.order.summary <- pantheria %>% select(MSW05_Order, `5-1_AdultBodyMass_g`) %>%
        mutate(`5-1_AdultBodyMass_g` = ifelse(`5-1_AdultBodyMass_g` == -999, NA, `5-1_AdultBodyMass_g`)) %>%
        group_by(MSW05_Order) %>%
        summarize(`5-1_AdultBodyMass_g.mean` = mean(`5-1_AdultBodyMass_g`, na.rm = TRUE),
            `5-1_AdultBodyMass_g.median` = median(`5-1_AdultBodyMass_g`, na.rm = TRUE), .groups="drop")
    mammals.pantheria <- mammals.pantheria %>%
        left_join(pantheria.order.summary %>% select(MSW05_Order, `5-1_AdultBodyMass_g.mean`), by=c("consensus.order" = "MSW05_Order")) %>%
        mutate(AdultBodyMass_g_source = ifelse(!is.na(AdultBodyMass_g_source), AdultBodyMass_g_source, ifelse(!is.na(`5-1_AdultBodyMass_g.mean`), "order mean", NA))) %>%
        mutate(AdultBodyMass_g = ifelse(!is.na(AdultBodyMass_g), AdultBodyMass_g, ifelse(!is.na(`5-1_AdultBodyMass_g.mean`), `5-1_AdultBodyMass_g.mean`, NA))) %>%
        select(-`5-1_AdultBodyMass_g.mean`)

    rm(pantheria.genus.summary, pantheria.family.summary, pantheria.order.summary)

    # add everything to leech
    leech <- leech %>%
        left_join(mammals.pantheria %>% select(consensus.short, domestic, AdultBodyMass_g, AdultBodyMass_g_source), by = "consensus.short")

# OTU level overlap between datasets

    # 46 OTUs overlap ... N.B. INCLUDES humans

    LSU.OTU.names <- leech %>% filter(dataset == "LSU") %>% pull(consensus.short) %>% unique()
    SSU.OTU.names <- leech %>% filter(dataset == "SSU") %>% pull(consensus.short) %>% unique()

    # shared OTU names
        LSU.OTU.names[LSU.OTU.names %in% SSU.OTU.names]
    # LSU only names
        LSU.OTU.names[!(LSU.OTU.names %in% SSU.OTU.names)]
    # SSU only names
        SSU.OTU.names[!(SSU.OTU.names %in% LSU.OTU.names)]

    rm(LSU.OTU.names, SSU.OTU.names)

# genus level overlap between datasets

    # 37 genera overlap ... N.B. INCLUDES Homo

    LSU.genus.names <- leech %>% filter(dataset == "LSU") %>% filter(consensus.genus != "NA") %>%
        pull(consensus.genus) %>% unique()
    SSU.genus.names <- leech %>% filter(dataset == "SSU") %>% filter(consensus.genus != "NA") %>%
        pull(consensus.genus) %>% unique()

    # shared OTU names
        LSU.genus.names[LSU.genus.names %in% SSU.genus.names]
    # LSU only names
        LSU.genus.names[!(LSU.genus.names %in% SSU.genus.names)]
    # SSU only names
        SSU.genus.names[!(SSU.genus.names %in% LSU.genus.names)]

    rm(LSU.genus.names, SSU.genus.names)

# ranger IDs

    # about three-quarters of the Ranger_IDs that were used in place of missing Polygon_IDs
    # were not present among the Ranger_IDs that were associated with Polygon_IDs

    LSU.rangers.with.PolygonID.missing <- leech %>% filter(is.na(Polygon_ID_original)) %>% pull(Ranger_ID) %>% unique()
    LSU.rangers.with.PolygonID.present <- leech %>% filter(!is.na(Polygon_ID_original)) %>% pull(Ranger_ID) %>% unique()
    LSU.rangers.with.PolygonID.missing %in% LSU.rangers.with.PolygonID.present %>% table()

    SSU.rangers.with.PolygonID.missing <- leech %>% filter(is.na(Polygon_ID_original)) %>% pull(Ranger_ID) %>% unique()
    SSU.rangers.with.PolygonID.present <- leech %>% filter(!is.na(Polygon_ID_original)) %>% pull(Ranger_ID) %>% unique()
    SSU.rangers.with.PolygonID.missing %in% SSU.rangers.with.PolygonID.present %>% table()

    rm(LSU.rangers.with.PolygonID.missing, LSU.rangers.with.PolygonID.present, SSU.rangers.with.PolygonID.missing, SSU.rangers.with.PolygonID.present)

# how many samples are missing environmental data?

    # all samples have Polygon_ID (though some are just the Ranger_ID
    # and thus do not connect up with the environmental data, and some are
    # ambiguous e.g. "28 or 27", "25 or 24" or "89 or 90"

        leech %>% filter(is.na(Polygon_ID))

    # how many Lab_IDs have latitude? (if they have latitude they have everything else)

        leech %>% select(-OTU, -reads, -starts_with("consensus"), -Chinese_common_name, -domestic, -AdultBodyMass_g, -AdultBodyMass_g_source) %>% distinct() %>%
            mutate(has.latitude = !is.na(latitude)) %>%
            select(dataset, has.latitude) %>% table()

    #        has.latitude
    # dataset FALSE TRUE
    #     LSU   179  557
    #     SSU   198  624

    # how many Polygon_IDs have latitude and/or reads?
    # - polygons with latitude but without reads are the extra rows added to represent polygons with missing data
    # - polygons without latitude but with reads are the dummy polygons (=rangerID) used to represent data without location information
    # - polygons with latitude and with reads are the properly labelled data

            leech %>%
            mutate(has.reads = !is.na(reads)) %>%
            select(-OTU, -reads, -starts_with("consensus"),  -Chinese_common_name, -domestic, -AdultBodyMass_g, -AdultBodyMass_g_source) %>% distinct() %>%
            select(dataset, latitude, Polygon_ID, has.reads) %>%
            group_by(dataset, Polygon_ID, has.reads) %>% summarize(has.latitude = !any(is.na(latitude)), .groups="drop") %>%
            select(dataset, has.reads, has.latitude) %>% table()

    # , , has.latitude = FALSE
    #
    #        has.reads
    # dataset FALSE TRUE
    #     LSU     0   37
    #     SSU     0   37
    #
    # , , has.latitude = TRUE
    #
    #        has.reads
    # dataset FALSE TRUE
    #     LSU    83   89
    #     SSU    82   90

# ubiquity of human reads

    # fraction of LabIDs with human reads
        leech %>%
            filter (!is.na(reads)) %>%
            filter(consensus.short == "Homo sapiens") %>%
            group_by(dataset) %>% summarize(`Lab_IDs with human reads` = sum(reads > 0) / n(), .groups="drop")

    # fraction of PolygonIDs with human reads
        leech %>%
            filter (!is.na(reads)) %>%
            filter(consensus.short == "Homo sapiens") %>%
            group_by(dataset, Polygon_ID) %>%
            summarize(humans.detected = ifelse(sum(reads) > 0, 1, 0), .groups="drop_last") %>%
            summarize(`Polygon_IDs with human reads` = sum(humans.detected) / n(), .groups="drop")

# fraction of total reads from humans

    # human reads are a much larger problem in the LSU dataset
    leech %>% mutate(human = ifelse(consensus.short == "Homo sapiens", "human reads", "non-human reads")) %>%
        group_by(dataset, human) %>%
        summarize(total.reads = sum(reads, na.rm = TRUE), .groups="drop_last") %>%
        spread(key = human, value = total.reads) %>%
        mutate(`total reads` = `human reads` + `non-human reads`) %>%
        mutate(`fraction human reads` = `human reads` / `total reads`)

# humans vs non-humans

    # non-human reads vs human reads
    # there is some tradeoff between number of human reads and number of non-human reads
    leech %>%
        filter (!is.na(reads)) %>%
        mutate(human = ifelse(consensus.short == "Homo sapiens", "human reads", "non-human reads")) %>%
        group_by(dataset, Lab_ID, human) %>%
        summarize(reads = sum(reads), .groups="drop_last") %>%
        spread(key = human, value = reads) %>%
        ggplot(aes(x = `human reads` + 1, y = `non-human reads` + 1)) +
            geom_jitter() + scale_x_log10() + scale_y_log10() + facet_wrap("dataset")

    # non-human species richness vs non-human reads
    leech %>%
        filter (!is.na(reads)) %>%
        filter(consensus.short != "Homo sapiens") %>%
        group_by(dataset, Lab_ID) %>%
        summarize(`non-human reads` = sum(reads), `non-human species richness` = sum(reads > 0), .groups="drop") %>%
        ggplot(aes(x = `non-human reads` + 1, y = `non-human species richness`)) +
            geom_jitter() + scale_x_log10() + facet_wrap("dataset")

    # non-human species richness vs human reads
    # the relationship between nonhuman species richness and number of human reads is less obvious,
    # but presumably there is a negative relationship given the tradeoff between human and nonhuman
    # reads, and the relationship between nonhuman species richness vs nonhuman reads
    nonhuman.richness <- leech %>% filter (!is.na(reads)) %>% filter(consensus.short != "Homo sapiens") %>%
        group_by(dataset, Lab_ID) %>% summarize(`non-human species richness` = sum(reads > 0), .groups="drop")
    human.reads <- leech %>% filter (!is.na(reads)) %>% filter(consensus.short == "Homo sapiens") %>%
        group_by(dataset, Lab_ID) %>% summarize(`human reads` = sum(reads), .groups="drop")
    full_join(nonhuman.richness, human.reads, by = c("dataset","Lab_ID")) %>%
        ggplot(aes(x = `human reads` + 1, y = `non-human species richness`)) +
            geom_jitter() + scale_x_log10() + facet_wrap("dataset")
    rm(nonhuman.richness, human.reads)

# split out human reads and extra Polygon_IDs without data

    # remove human reads
    leech <- leech %>% filter(consensus.short!= "Homo sapiens")

    # separate the extra rows we added earlier (minus Homo sapiens)
    leech.supplement <- leech %>% filter(is.na(reads)) #10,801 rows
    leech <- leech %>% filter(!is.na(reads)) #91807 rows

# basic info about raw data

    # number of OTUs in each dataset
        leech %>% filter(dataset == "LSU") %>% select(OTU) %>% unique() %>% nrow() # 59
        leech %>% filter(dataset == "SSU") %>% select(OTU) %>% unique() %>% nrow() # 72

    # number of Lab_IDs in each dataset
        leech %>% filter(dataset == "LSU") %>% select(Lab_ID) %>% unique() %>% nrow() # 653
        leech %>% filter(dataset == "SSU") %>% select(Lab_ID) %>% unique() %>% nrow() # 740

    # number of Polygon_IDs in each dataset
        leech %>% filter(dataset == "LSU") %>% select(Polygon_ID) %>% unique() %>% nrow() # 126
        leech %>% filter(dataset == "SSU") %>% select(Polygon_ID) %>% unique() %>% nrow() # 127

# how many reads now?

    leech %>% group_by(dataset) %>% summarize(reads = sum(reads), .groups="drop") # LSU: 18502593 (from 77092321!), SSU: 84951011 (from 89943206)

# accounting for all the reads and the leeches in our dataset ...

    # read data as initially imported
        read.tbl$LSU %>% select(-queryID) %>% sum() # LSU: 82653202 reads
        read.tbl$SSU %>% select(-queryID) %>% sum() # SSU: 91132342 reads
    # leech quantity data as initially imported
        leech_qty %>% select(leech_qty) %>% sum() # 30468 leeches

    # read data after removing control samples
        read.tidy$LSU %>% select(reads) %>% sum() # LSU: 77772174 reads
        read.tidy$SSU %>% select(reads) %>% sum() # SSU: 90182071 reads
    # leech quantities after matching to reads as imported
        leech_qty[leech_qty$Lab_ID %in% colnames(read.tbl$LSU),] %>% select(leech_qty) %>% sum() # LSU: 23818 leeches
        leech_qty[leech_qty$Lab_ID %in% colnames(read.tbl$SSU),] %>% select(leech_qty) %>% sum() # SSU: 26804 leeches
    # leech quantities after matching to non-control-sample reads (excluding controls shouldn't change anything here)
        leech_qty[leech_qty$Lab_ID %in% read.tidy$LSU$Lab_ID,] %>% select(leech_qty) %>% sum() # LSU: 23818 leeches
        leech_qty[leech_qty$Lab_ID %in% read.tidy$SSU$Lab_ID,] %>% select(leech_qty) %>% sum() # SSU: 26804 leeches
    # leech quantities that don't match to read data (i.e. these samples failed)
        leech_qty[!(leech_qty$Lab_ID %in% colnames(read.tbl$LSU)),] %>% select(leech_qty) %>% sum() # LSU: 6650 leeches
        leech_qty[!(leech_qty$Lab_ID %in% colnames(read.tbl$SSU)),] %>% select(leech_qty) %>% sum() # SSU: 3664 leeches

    # read data after removing unmatched OTUs (i.e. no taxonomy information)
        OTUs$LSU %>% select(reads) %>% sum() # LSU: 77092321 reads
        OTUs$SSU %>% select(reads) %>% sum() # SSU: 89984501 reads
    # and as a check, read data belonging to unmatched OTUs (i.e. no taxonomy information) - these should add to totals above
        preOTUs.unmatched$LSU %>% select(reads) %>% sum() # LSU: 679853 reads
        preOTUs.unmatched$SSU %>% select(reads) %>% sum() # SSU: 197570 reads
    # matching OTUs with taxonomy information has no bearing on leech quantities included
        leech_qty[leech_qty$Lab_ID %in% OTUs$LSU$Lab_ID,] %>% select(leech_qty) %>% sum() # LSU: 23818  leeches
        leech_qty[leech_qty$Lab_ID %in% OTUs$SSU$Lab_ID,] %>% select(leech_qty) %>% sum() # SSU: 26804 leeches
        leech_qty[leech_qty$Lab_ID %in% preOTUs.unmatched$LSU$Lab_ID,] %>% select(leech_qty) %>% sum() # LSU: 23818 leeches
        leech_qty[leech_qty$Lab_ID %in% preOTUs.unmatched$SSU$Lab_ID,] %>% select(leech_qty) %>% sum() # SSU: 26804 leeches

    # excluding ZY47 (missing polygon and ranger information) from both datasets and ZY22.4 from LSU (no reads in LSU dataset
    # after removing OTUs without taxonomic match, though there are reads in SSU) brings down the leech count included in the two datasets set slightly:
        leech_qty[leech_qty$Lab_ID %in% (OTUs$LSU %>% pull(Lab_ID)),] %>%
            filter(Lab_ID != "ZY47" & Lab_ID != "ZY22.4") %>% select(leech_qty) %>% sum() # LSU: 23774 leeches
        leech_qty[leech_qty$Lab_ID %in% (OTUs$SSU %>% pull(Lab_ID)),] %>%
            filter(Lab_ID != "ZY47") %>% select(leech_qty) %>% sum() # LSU: 26775 leeches
    # and this now matches the leech counts in the final dataset:
        leech %>% group_by(dataset) %>% distinct(dataset, Lab_ID, leech_qty) %>%
            summarize(total.leeches = sum(leech_qty), .groups="drop") # LSU: 23774 leeches, SSU: 26775 leeches
    # excluding ZY47 also brings down the SSU read count slightly (ZY22.4 makes no difference since it had zero reads at this point):
        OTUs$LSU %>% filter(Lab_ID != "ZY47" & Lab_ID != "ZY22.4") %>% select(reads) %>% sum() # LSU: 77092321 reads still
        OTUs$SSU %>% filter(Lab_ID != "ZY47") %>% select(reads) %>% sum() # LSU: 89943206 reads

    # additionally excluding humans also brings down the LSU and SSU read counts further:
        OTUs$LSU %>% filter(Lab_ID != "ZY47" & Lab_ID != "ZY22.4") %>%
            filter(consensus.short != "Homo sapiens") %>% select(reads) %>% sum() # LSU: 18502593 reads
        OTUs$SSU %>% filter(Lab_ID != "ZY47") %>%
            filter(consensus.short != "Homo sapiens") %>% select(reads) %>% sum() # LSU: 84951011 reads
    # and this now matches the read counts in the final dataset:
        leech %>% group_by(dataset) %>% summarize(total.reads = sum(reads), .groups="drop") # LSU: 18502593 reads, SSU: 84951011 reads

# export data

    # save as R objects

    # leech is the main data tibble; leech.supplement is extra rows to represent Polygon_IDs with missing data
    # the two tibbles can be joined with bind_rows to create the full augmented dataset
        save(leech, leech.supplement, file = here("rdata","Ailaoshan_OTU_table.rdata"))

    # preOTUs.matched is the list of preOTUs with taxonomic data
    # -- use this to relabel the FASTA file of representative sequences
        preOTUs.matched.summary <- lapply(preOTUs.matched, FUN = function (X) X %>% select(-Lab_ID, -reads) %>% distinct())
        save(preOTUs.matched.summary, file = here("rdata","Ailaoshan_OTU_table_preOTUs.rdata"))
