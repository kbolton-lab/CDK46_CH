# Code used to generate the minimum required dataframes for the CH project

library(data.table)
library(tidyverse)

raw_data_file_path <- "/Users/irenaeuschan/Documents/Irenaeus/CDK46_CH/data/RawData/"

sclc_df <- fread(paste0(raw_data_file_path, "g1_sclc_df.csv"))                       # Small Cell Lung Cancer
tnbc_df <- fread(paste0(raw_data_file_path, "g1_tnbc_df.csv"))                       # Triple Negative Breast Cancer
crc_df <- fread(paste0(raw_data_file_path, "g1_crc_df.csv"))                         # Colorectal Cancer
untreated_df <- fread(paste0(raw_data_file_path, "Untreated.csv"))                   # Untreated

# What are the minimum information for the CH dataframe?
sclc_df %>%
    separate(SampleDeID, c("patientID", "sampleID"), sep = "_") %>%
    dplyr::rename(SampleDeID = sampleID) %>%
    mutate(
        oncoKB_reviewed = NA,
        n.loci.truncating.vep = NA,
        source.totals.loci = NA,
        source.totals.loci.truncating = NA,
        source.totals.p = NA,
        source.totals.c = NA,
        whereincycle = case_when(
            whereincycle == "pre" ~ "C1D1",
            whereincycle == "post" ~ "C5D1",
            whereincycle == "post-post" ~ "90D Post-Maintenance"
        ),
        whichdraw = case_when(
            whereincycle == "C1D1" ~ "PreTx",
            whereincycle == "C5D1" ~ "PostTx",
            whereincycle == "90D Post-Maintenance" ~ "PostMaintenance"
        )
    ) %>%
    select(
        key, SampleDeID, DeID, sample_key, Gene, average_AF, average_AD, pd_reason, putative_driver, 
        whereincycle, Trilaciclib, whichdraw,
        gene_aachange, gene_cDNAchange, n.loci.vep, n.loci.truncating.vep, n.HGVSp, n.HGVSc, COSMIC_ID, CosmicCount, heme_cosmic_count, myeloid_cosmic_count, oncoKB, oncoKB_reviewed, isOncogenic, isTSG, isTruncatingHotSpot, VariantClass, nsamples,
        source.totals.loci, source.totals.loci.truncating, source.totals.p, source.totals.c
    ) %>%
    mutate(
        key = gsub(" |>", ":", key),
        sample_key = paste0(SampleDeID, "_", key)
    ) %>% fwrite("~/Documents/irenaeus/CDK46_CH/data/g1_sclc_df.minimum.csv")

untreated_df %>%
    mutate(
        Trilaciclib = "Untreated",
        pd_reason = "Old Data Pass",
        putative_driver = ch_pd_ic,
        n.HGVSp = n.HGVSp.x,
        n.HGVSc = n.HGVSc.x,
        oncoKB_reviewed = NA,
        n.loci.truncating.vep = NA,
        source.totals.loci = NA,
        source.totals.loci.truncating = NA,
        source.totals.p = NA,
        source.totals.c = NA,
        whereincycle = case_when(
            whereincycle == "pre" ~ "PreSample",
            whereincycle == "post" ~ "PostSample"
        ),
        whichdraw = if_else(whichdraw == 1, "PreTx", "PostTx")
    ) %>%
    select(
        key, SampleDeID, DeID, sample_key, Gene, average_AF, average_AD, pd_reason, putative_driver, 
        whereincycle, Trilaciclib, whichdraw,
        gene_aachange, gene_cDNAchange, n.loci.vep, n.loci.truncating.vep, n.HGVSp, n.HGVSc, COSMIC_ID, CosmicCount, heme_cosmic_count, myeloid_cosmic_count, oncoKB, oncoKB_reviewed, isOncogenic, isTSG, isTruncatingHotSpot, VariantClass, nsamples,
        source.totals.loci, source.totals.loci.truncating, source.totals.p, source.totals.c
    ) %>%
    mutate(
        key = gsub(" |>", ":", key),
        sample_key = paste0(SampleDeID, "_", key)
    ) %>% fwrite("~/Documents/irenaeus/CDK46_CH/data/untreated_df.minimum.csv")

tnbc_df %>%
    mutate(
        oncoKB_reviewed = NA,
        whereincycle = case_when(
            whereincycle == "pre" ~ "C1D1",
            whereincycle == "during" ~ "C7D1",
            whereincycle == "post" ~ "Post Tx",
            whereincycle == "postpost" ~ "60D Post-Maintenance"
        ),
        whichdraw = case_when(
            whereincycle == "C1D1" ~ "PreTx",
            whereincycle == "C7D1" ~ "MidTx",
            whereincycle == "Post Tx" ~ "PostTx",
            whereincycle == "60D Post-Maintenance" ~ "PostMaintenance"
        ),
        Trilaciclib = case_when(
            Trilaciclib == "trila" ~ "Trilaciclib",
            Trilaciclib == "placebo" ~ "Placebo"
        )
    ) %>%
    select(
        key, SampleDeID, DeID, sample_key, Gene, average_AF, average_AD, pd_reason, putative_driver, 
        whereincycle, Trilaciclib, whichdraw,
        gene_aachange, gene_cDNAchange, n.loci.vep, n.loci.truncating.vep, n.HGVSp, n.HGVSc, COSMIC_ID, CosmicCount, heme_cosmic_count, myeloid_cosmic_count, oncoKB, oncoKB_reviewed, isOncogenic, isTSG, isTruncatingHotSpot, VariantClass, nsamples,
        source.totals.loci, source.totals.loci.truncating, source.totals.p, source.totals.c
    ) %>% fwrite("~/Documents/irenaeus/CDK46_CH/data/g1_tnbc_df.minimum.csv")

crc_df %>% 
    filter(!grepl("NTN", patientID)) %>%
    dplyr::rename(SampleDeID = sampleID, DeID = patientID) %>%
    mutate(
        whereincycle = if_else(whereincycle == "MAINTENANCE C1D1", "C1D1 Maintenance", whereincycle),
        whichdraw = case_when(
            whereincycle == "C1D1" ~ "PreTx",
            whereincycle == "C5D1" ~ "MidTx",
            whereincycle == "C1D1 Maintenance" ~ "PostTx"
        )
    ) %>%
    select(
        key, SampleDeID, DeID, sample_key, Gene, average_AF, average_AD, pd_reason, putative_driver, 
        whereincycle, Trilaciclib, whichdraw,
        gene_aachange, gene_cDNAchange, n.loci.vep, n.loci.truncating.vep, n.HGVSp, n.HGVSc, COSMIC_ID, CosmicCount, heme_cosmic_count, myeloid_cosmic_count, oncoKB, oncoKB_reviewed, isOncogenic, isTSG, isTruncatingHotSpot, VariantClass, nsamples,
        source.totals.loci, source.totals.loci.truncating, source.totals.p, source.totals.c
    ) %>% fwrite("~/Documents/irenaeus/CDK46_CH/data/g1_crc_df.minimum.csv")

sclc_tp <- fread(paste0(raw_data_file_path, "sclc_timepoints.csv"))                       # Small Cell Lung Cancer
tnbc_tp <- fread(paste0(raw_data_file_path, "tnbc_timepoints.csv"))                       # Triple Negative Breast Cancer
crc_tp <- fread(paste0(raw_data_file_path, "crc_timepoints.csv"))                         # Colorectal Cancer

sclc_tp %>%
    mutate(
        duringAge = NA,
    ) %>%
    select(
        deid_key, DeID, key, Gene, Trilaciclib,
        all_time_points, prepreVAF, preVAF, duringVAF, postVAF, postpostVAF,
        rescued, missing, called, at_limit,
        preAge, duringAge, postAge, postpostAge, in_cohort,
        DDR, Splice, DTA, Gene_Class,
        change_in_days, growth_rate, change_in_VAF, compare_group
    ) %>% 
    mutate(
        key = gsub(" |>", ":", key),
        compare_group = case_when(
            compare_group == "PrePost" ~ "C1D1vsC5D1",
            compare_group == "PrePostPost" ~ "C1D1vs90DPost-Maintenance",
            compare_group == "" ~ "PrevsPost"
        )    
    ) %>% fwrite("~/Documents/irenaeus/CDK46_CH/data/sclc_timepoints.minimum.csv")
tnbc_tp %>%
    filter(fromwhere == "G1 Therapeutics") %>%
    mutate(
        Trilaciclib = case_when(
            Trilaciclib == "trila" ~ "Trilaciclib",
            Trilaciclib == "placebo" ~ "Placebo"
        )
    ) %>%
    select(
        deid_key, DeID, key, Gene, Trilaciclib,
        all_time_points, prepreVAF, preVAF, duringVAF, postVAF, postpostVAF,
        rescued, missing, called, at_limit,
        preAge, duringAge, postAge, postpostAge, in_cohort,
        DDR, Splice, DTA, Gene_Class,
        change_in_days_during, growth_rate_during, change_in_VAF, compare_group
    ) %>% mutate(
        key = gsub(" |>", ":", key),
        compare_group = case_when(
            compare_group == "PreDuring" ~ "C1D1vsC7D1",
            compare_group == "PrePost" ~ "C1D1vsMaintenance"
        )    
    ) %>%
    dplyr::rename(
        growth_rate = growth_rate_during,
        change_in_days = change_in_days_during
    ) %>% fwrite("~/Documents/irenaeus/CDK46_CH/data/tnbc_timepoints.minimum.csv")
crc_tp %>%
    filter(fromwhere == "G1 Therapeutics") %>%
    mutate(
        duringAge = NA,
        postpostAge = NA,
    ) %>%
    select(
        deid_key, DeID, key, Gene, Trilaciclib,
        all_time_points, prepreVAF, preVAF, duringVAF, postVAF, postpostVAF,
        rescued, missing, called, at_limit,
        preAge, duringAge, postAge, postpostAge, in_cohort,
        DDR, Splice, DTA, Gene_Class,
        change_in_days, growth_rate, change_in_VAF, compare_group
    ) %>% fwrite("~/Documents/irenaeus/CDK46_CH/data/crc_timepoints.minimum.csv")
#fwrite(g1_tp, "~/Documents/irenaeus/CDK46_CH/data/g1_timepoints.csv")

g1_demo <- fread("/Users/irenaeuschan/Documents/Irenaeus/CDK4:6/OldData/g1.csv")
g1_demo %>%
select(
    DEID, OurID,
    Age, SEX, RACE, ETHNIC
) %>%
dplyr::rename(SampleDeID = OurID, DeID = DEID)

g1_info_extra <- read_excel("~/Documents/Irenaeus/ArcherDX/G1_manifests/WU CH sample pull list 13Feb23.xlsx", sheet = 2)
g1_info_extra <- g1_info_extra %>% 
    select(`Originating ID`, `Site #`, `Subject #`, `Time Point`, AGE, Trilaciclib) %>%
    dplyr::rename(
        sample_id = `Originating ID`,
        site = `Site #`,
        subject = `Subject #`,
        timepoint = `Time Point`
    )
g1_info_extra <- g1_info_extra %>%
    mutate(
        DeID = gsub("-", "", subject)
    )

g1_info_tnbc <- read_excel("~/Documents/Irenaeus/ArcherDX/G1_manifests/WU CH sample pull list 13Feb23.xlsx", sheet = 1)
g1_info_tnbc <- g1_info_tnbc %>%
    select(`Originating ID`, `Site #`, `Subject ID`, whereincycle, Trilaciclib, Cycles_Received, Cycle_1, Cycle_7, Post, Post_60, Months_From, Age_Start, Gender) %>%
    mutate(DeID = gsub("-", "", `Subject ID`)) %>%
    dplyr::rename(
        SampleDeID = `Originating ID`,
        site = `Site #`
    )
g1_info_tnbc <- g1_info_tnbc %>%
    mutate(
        Age = case_when(
            whereincycle == "pre" ~ Age_Start,
            whereincycle == "during" ~ Age_Start + (4.8/12),
            whereincycle == "post" & Post_60 == "no" ~ Age_Start + (Months_From/12),
            whereincycle == "post" & Post_60 == "yes" ~ Age_Start + ((Months_From-2)/12),   # Since post_60 is 60 days after post (2 months) post should be ~postpost - 2 months
            whereincycle == "postpost" ~ Age_Start + (Months_From/12)
        )
    )

g1_crc_bisi_sample_info <- fread("/Users/irenaeuschan/Documents/Irenaeus/CDK4:6/g1_bisi_keylist_edit.csv")
