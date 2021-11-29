# R code to export FASTA files for publication

library("here")
library("tidyverse")
library("seqinr")

# get preOTU data

    load(file = here("rdata", "Ailaoshan_OTU_table_preOTUs.rdata"))  # generated with Ailaoshan_OTU_table.R

# read in sequence data

    rep.seqs <- list(
        LSU = read.fasta(file = here("data", "fasta", "LSU_16S_otu_table_swarm_lulu_20190624.fas"),
            as.string = TRUE, forceDNAtolower = FALSE, strip.desc = TRUE),
        SSU = read.fasta(file = here("data", "fasta", "SSU_12S_otu_table_swarm_lulu_20190624.fas"),
            as.string = TRUE, forceDNAtolower = FALSE, strip.desc = TRUE)
    )

# filter sequences from FASTA files
# remove sequences from FASTA files without matches in preOTUs.matched.summary
# (these were discarded from the analysis, so don't include them in FASTA)

    filtered.rep.seqs <- sapply(names(rep.seqs), simplify = FALSE, function (X) {
        rep.seqs[[X]][names(rep.seqs[[X]]) %in% preOTUs.matched.summary[[X]]$queryID]
    })

# get new metadata to be used on output FASTA file

    parsed.rep.seqs <- sapply(names(filtered.rep.seqs), simplify = FALSE, function (X) {
        # start with queryID from the FASTA file (queryID is the preOTU identifier)
            tibble(queryID = getAnnot(filtered.rep.seqs[[X]]) %>% unlist(),
                sequence = getSequence(filtered.rep.seqs[[X]], as.string= TRUE) %>% unlist()) %>%
        # add correponding OTU number (OTU is the OTU identifier) and taxonomy
            left_join(preOTUs.matched.summary[[X]] %>% select(queryID, OTU, consensus.short), by = "queryID") %>%
        # add suffixes to OTU number for OTUs present multiple times (OTU.suffix is augmented OTU identifier)
            # add dummy and counts within each OTU
                group_by(OTU) %>%
                mutate(dummy = 1, count = n()) %>%
            # determine suffixes for OTUs present multiple times
                mutate(suffix = cumsum(dummy)) %>%
                mutate(suffix = paste0(".", suffix)) %>%
                mutate(suffix = ifelse(count == 1, "", suffix)) %>%
                ungroup() %>%
            # add suffix to existing OTU name
                mutate(OTU.suffix = paste0(OTU, suffix)) %>%
        # make deflines for new FASTA file
            mutate(defline = paste0(OTU.suffix, " ", consensus.short)) %>%
        # clean up
            select(queryID, OTU, OTU.suffix, consensus.short, defline, sequence) %>%
            arrange(OTU.suffix)
    })

#  write sequences to FASTA format files

    ## LSU ##

        LSU.filename <- here("fasta", "Supplementary_Data_4_LSU_representative_sequences.fasta")

        LSU.comment1 <- "# Ailaoshan representative sequences (LSU dataset)"
        LSU.comment2 <- "# The first six characters of each sequence identifier (e.g. LSU001) indicate the OTU that the sequence comes from. Identifiers with a suffix (e.g. LSU005.1) refer to OTUs that were formed by merging two or more pre-OTUs. In such cases, the representative sequence for each pre-OTU is provided, with the suffix distinguishing between them (e.g. LSU005.1 and LSU005.2 represent the two pre-OTUs that were merged to form the OTU LSU005)."

        cat(LSU.comment1, LSU.comment2, file = LSU.filename, sep = "\n", append = FALSE)
        write.fasta(sequences = as.list(parsed.rep.seqs$LSU$sequence),
            names = parsed.rep.seqs$LSU$defline, file.out = LSU.filename, open = "a")

    ## SSU ##

        SSU.filename <- here("fasta", "Supplementary_Data_5_SSU_representative_sequences.fasta")

        SSU.comment1 <- "# Ailaoshan representative sequences (SSU dataset)"
        SSU.comment2 <- "# The first six characters of each sequence identifier (e.g. SSU001) indicate the OTU that the sequence comes from. Identifiers with a suffix (e.g. SSU002.1) refer to OTUs that were formed by merging two or more pre-OTUs. In such cases, the representative sequence for each pre-OTU is provided, with the suffix distinguishing between them (e.g. SSU002.1 and SSU002.2 represent the two pre-OTUs that were merged to form the OTU SSU002)."

        cat(SSU.comment1, SSU.comment2, file = SSU.filename, sep="\n", append=FALSE)
        write.fasta(sequences = as.list(parsed.rep.seqs$SSU$sequence),
            names = parsed.rep.seqs$SSU$defline, file.out = SSU.filename, open = "a")
