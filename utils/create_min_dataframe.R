# Code used to generate the minimum required dataframes for the CH project

library(data.table)
library(tidyverse)

raw_data_file_path <- "/Users/chani/Library/CloudStorage/Box-Box/Irenaeus/CDK46_CH/data/RawData/"

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
    ) #%>% fwrite("~/Documents/irenaeus/CDK46_CH/data/g1_sclc_df.minimum.csv")

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
    dplyr::select(
        deid_key, DeID, key, Gene, Trilaciclib,
        all_time_points, prepreVAF, preVAF, duringVAF, postVAF, postpostVAF,
        rescued, missing, called, at_limit,
        preAge, duringAge, postAge, postpostAge, in_cohort,
        DDR, Splice, DTA, Gene_Class,
        change_in_days_during, growth_rate_during, change_in_VAF, compare_group_during
    ) %>% mutate(
        key = gsub(" |>", ":", key),
        compare_group_during = case_when(
            compare_group_during == "PreDuring" ~ "C1D1vsC7D1",
            compare_group_during == "PrePost" ~ "C1D1vsMaintenance"
        )    
    ) %>%
    dplyr::rename(
        growth_rate = growth_rate_during,
        change_in_days = change_in_days_during,
        compare_group = compare_group_during
    ) %>% fwrite("~/Documents/irenaeus/CDK46_CH/data/tnbc_timepoints.minimum.new.csv")
crc_tp %>%
    filter(fromwhere == "G1 Therapeutics") %>%
    mutate(
        duringAge = NA,
        postpostAge = NA,
    ) %>%
    dplyr::select(
        deid_key, DeID, key, Gene, Trilaciclib,
        all_time_points, prepreVAF, preVAF, duringVAF, postVAF, postpostVAF,
        rescued, missing, called, at_limit,
        preAge, duringAge, postAge, postpostAge, in_cohort,
        DDR, Splice, DTA, Gene_Class,
        change_in_days, growth_rate, change_in_VAF, compare_group
    ) %>% fwrite("~/Documents/irenaeus/CDK46_CH/data/crc_timepoints.minimum.csv")
#fwrite(g1_tp, "~/Documents/irenaeus/CDK46_CH/data/g1_timepoints.csv")

keylist <- fread("/Users/irenaeuschan/Documents/Irenaeus/ArcherDX/data/keylistall6.csv")
table_parp <- fread("/Users/irenaeuschan/Documents/Irenaeus/ArcherDX/data/table1.csv")
g1_csv <- fread("/Users/irenaeuschan/Documents/Irenaeus/CDK4:6/OldData/g1.csv")
carlos_table <- fread("/Users/irenaeuschan/Documents/Irenaeus/ArcherDX/data/carlos_table1.csv")
library(readxl)
g1_info_extra <- read_excel("~/Documents/Irenaeus/ArcherDX/G1_manifests/WU CH sample pull list 13Feb23.xlsx", sheet = 2)
g1_info_tnbc <- read_excel("~/Documents/Irenaeus/ArcherDX/G1_manifests/WU CH sample pull list 13Feb23.xlsx", sheet = 1)
g1_crc_bisi_sample_info <- fread("/Users/irenaeuschan/Documents/Irenaeus/CDK4:6/g1_bisi_keylist_edit.csv")

format_date_vector <- function(datetime_vector) {
    sapply(datetime_vector, function(x) {
        tryCatch(
            {
                if (is.na(x)) {
                    return(NA)
                }
                # Split on space and take first part
                date_part <- strsplit(x, " ")[[1]][1]
                # Split date and pad with zeros
                date_components <- strsplit(date_part, "/")[[1]]
                month <- sprintf("%02d", as.numeric(date_components[1]))
                day <- sprintf("%02d", as.numeric(date_components[2]))
                year <- date_components[3]
                # Return formatted date
                paste(year, month, day, sep = "-")
            },
            error = function(e) NA
        )
    })
}

get_random_date_from_deid <- function(deid) {
    # Use DEID as seed
    set.seed(as.numeric(deid))
    # Generate random date
    start_date <- as.Date("2000-01-01")
    end_date <- as.Date("2022-12-31")
    start_date + sample(0:(end_date - start_date), 1)
}

rbind(
    carlos_table %>%
    mutate(
        TimePoint = ifelse(whichdraw == 1, "TP1", "TP2"),
        SampleCollectionDate = case_when(
            whichdraw == 1 ~ get_random_date_from_deid(gsub("PA", "", DEID)),
            whichdraw == 2 ~ get_random_date_from_deid(gsub("PA", "", DEID)) + CycleLength
        )
    ) %>%
    dplyr::select(DEID, SampleDeID, SampleCollectionDate, TimePoint),
    g1_csv %>%
        mutate(
            TimePoint = ifelse(whichdraw == 1, "C1D1", "C5D1"),
            SampleDeID = sapply(strsplit(SampleDeID, "_"), `[`, 2)
        ) %>%
        dplyr::select(DEID, SampleDeID, SampleCollectionDate, TimePoint),
    g1_info_extra %>%
        mutate(
            DEID = gsub("-", "", `Subject #`),
            SampleDeID = as.character(`Originating ID`),
            TimePoint = case_when(
                grepl("Cycle 1", `Time Point`) ~ "C1D1",
                grepl("C1D1", `Time Point`) ~ "C1D1",
                grepl("Cycle 5", `Time Point`) ~ "C5D1",
                grepl("C5D1", `Time Point`) ~ "C5D1",
                grepl("D90", `Time Point`) ~ "Post D90",
                grepl("Day 90", `Time Point`) ~ "Post D90"
            ),
            SampleCollectionDate = as.Date(format_date_vector(`Received In`))
        ) %>%
        dplyr::select(DEID, SampleDeID, SampleCollectionDate, TimePoint),
    g1_info_tnbc %>%
        mutate(
            DEID = gsub("-", "", `Subject ID`),
            SampleDeID = `Originating ID`,
            TimePoint = case_when(
                whereincycle == "pre" ~ "C1D1",
                whereincycle == "during" ~ "C7D1",
                whereincycle == "post" ~ "Post",
                whereincycle == "postpost" ~ "Post D60"
            ),
            SampleCollectionDate = as.Date(format_date_vector(`Received In`))
        ) %>%
        dplyr::select(DEID, SampleDeID, SampleCollectionDate, TimePoint),
    g1_crc_bisi_sample_info %>%
        dplyr::filter(fromwhere == "G1 Therapeutics - CRC") %>%
        mutate(
            DEID = patientID,
            SampleDeID = sampleID,
            TimePoint = case_when(
                whereincycle == "MAINTENANCE C1D1" ~ "C1D1 Maintenance",
                whereincycle == "C1D1" ~ "C1D1",
                whereincycle == "C5D1" ~ "C5D1"
            ),
            SampleCollectionDate = as.Date(format(as.Date(CollectionDate, "%d-%b-%y"), "%Y-%m-%d"))
        ) %>%
        dplyr::select(DEID, SampleDeID, SampleCollectionDate, TimePoint)
) #%>% fwrite("~/Documents/irenaeus/CDK46_CH/data/g1_sra_extra_columns.csv")

untreated_demo <- carlos_table %>%
    dplyr::select(
        DEID, SampleDeID, whichdraw, 
        SEX, RACE, Ethnicity, TOBAC30, Age_at_BloodDraw
    ) %>%
    dplyr::rename(
        DeID = DEID,
        Age = Age_at_BloodDraw,
        Sex = SEX,
        Race = RACE,
        SmokingStatus = TOBAC30
    ) %>%
    mutate(
        whichdraw = ifelse(whichdraw == 1, "PreTx", "PostTx"),
        SampleDeID = as.character(SampleDeID),
        SmokingStatus = case_when(
            SmokingStatus == 1 ~ "Current",
            SmokingStatus == 0 ~ "Never",
            is.na(SmokingStatus) ~ "Unknown"
        ),
        Trilaciclib = "Untreated",
        Cohort = "Untreated"
    )

g1_kelly <- g1_csv %>%
    dplyr::select(
        DEID, SampleDeID, whichdraw, `BioInventory Registration Number`, `BioInventory Project Name`, 
        `Sample Type`, `Protocol #`, `Site Number`, `Time Point`
    ) %>%
    dplyr::rename(
        DeID = DEID,
        `BioInventory Group Name` = `Sample Type`,
        `Site #` = `Site Number`
    ) %>%
    mutate(
        `Originating ID` = SampleDeID,
        Division = "",
        whichdraw = ifelse(whichdraw == 1, "PreTx", "PostTx"),
        Cohort = "SCLC"
    )
g1_kelly$`Subject #` <- paste0(g1_kelly$`Site #`, "-", substr(g1_kelly$DeID, nchar(g1_kelly$`Site #`) + 1, nchar(g1_kelly$DeID)))

g1_kelly <- rbind(
    g1_kelly,
    g1_info_extra %>% 
        mutate(
            DeID = gsub("-", "", `Subject #`),
            SampleDeID = as.character(`Originating ID`),
            whichdraw = case_when(
                grepl("Cycle 1", `Time Point`) ~ "PreTx",
                grepl("C1D1", `Time Point`) ~ "PreTx",
                grepl("Cycle 5", `Time Point`) ~ "PostTx",
                grepl("C5D1", `Time Point`) ~ "PostTx",
                grepl("D90", `Time Point`) ~ "PostPostMaintenance",
                grepl("Day 90", `Time Point`) ~ "PostPostMaintenance"
            )
        ) %>%
        dplyr::select(
            DeID, SampleDeID, whichdraw, `BioInventory Registration Number`, `BioInventory Project Name`,
            `BioInventory Group Name`, `Protocol #`, 
            `Site #`, `Time Point`, `Subject #`, `Originating ID`, Division
        ) %>%
        mutate(Cohort = "SCLC")
)

g1_kelly <- rbind(
    g1_kelly,
    g1_info_tnbc %>%
        mutate(
            DeID = gsub("-", "", `Subject ID`),
            SampleDeID = `Originating ID`,
            whichdraw = case_when(
                whereincycle == "pre" ~ "PreTx",
                whereincycle == "during" ~ "MidTx",
                whereincycle == "post" ~ "PostTx",
                whereincycle == "postpost" ~ "PostPostMaintenance"
            )
        ) %>%
        dplyr::select(
            DeID, SampleDeID, whichdraw, `BioInventory Registration Number`, `BioInventory Project Name`,
            `BioInventory Group Name`, `Protocol #`,
            `Site #`, `Time Point`, `Subject ID`, `Originating ID`, Division
        ) %>%
        dplyr::rename(
            `Subject #` = `Subject ID`
        ) %>%
        mutate(Cohort = "TNBC")
)

g1_kelly <- rbind(
    g1_kelly,
    g1_crc_bisi_sample_info %>%
        filter(fromwhere == "G1 Therapeutics - CRC") %>%
        mutate(
            `Protocol #` = "G1T28-207",
            `BioInventory Project Name` = "G1T28-207 CRC",
            DeID = patientID,
            `Subject #` = DeID,
            `Site #` = Site,
            `Time Point` = whereincycle,
            SampleDeID = sampleID,
                whichdraw = case_when(
                whereincycle == "MAINTENANCE C1D1" ~ "PostTx",
                whereincycle == "C1D1" ~ "PreTx",
                whereincycle == "C5D1" ~ "MidTx"
            ),
            `BioInventory Registration Number` = SampleDeID,
            `BioInventory Group Name` = "",
            `Originating ID` = "",
            Division = ""
        ) %>%
        dplyr::select(
            DeID, SampleDeID, whichdraw, `BioInventory Registration Number`, `BioInventory Project Name`,
            `BioInventory Group Name`, `Protocol #`,
            `Site #`, `Time Point`, `Subject #`, `Originating ID`, Division
        ) %>%
        mutate(Cohort = "CRC")
)

g1_kelly %>% 
    filter(Cohort != "Untreated") %>% 
    fwrite("G1_Samples_Request.csv")
    filter(ifelse(Cohort == "TNBC", whichdraw == "PreTx" | whichdraw == "MidTx", whichdraw == "PreTx" | whichdraw == "PostTx")) %>%
    filter(DeID %in% c(g1_kelly %>% 
                        filter(Cohort != "Untreated") %>% 
                        filter(ifelse(Cohort == "TNBC", whichdraw == "PreTx" | whichdraw == "MidTx", whichdraw == "PreTx" | whichdraw == "PostTx")) %>%
                        group_by(Cohort, DeID) %>%
                        summarise(n = n()) %>%
                        filter(n > 1) %>%
                        pull(DeID))
    ) %>% fwrite("G1_Samples_Request.csv")

g1_kelly %>% 
    filter(Cohort != "Untreated") %>% 
    filter(ifelse(Cohort == "TNBC", whichdraw == "PreTx" | whichdraw == "MidTx", whichdraw == "PreTx" | whichdraw == "PostTx")) %>%
    filter(DeID %in% c(g1_kelly %>% 
                        filter(Cohort != "Untreated") %>% 
                        filter(ifelse(Cohort == "TNBC", whichdraw == "PreTx" | whichdraw == "MidTx", whichdraw == "PreTx" | whichdraw == "PostTx")) %>%
                        group_by(Cohort, DeID) %>%
                        summarise(n = n()) %>%
                        filter(n > 1) %>%
                        pull(DeID))
    ) %>% filter(Cohort == "CRC") %>% dplyr::select(DeID) %>% unique()


sclc_demo <- rbind(
g1_csv %>%
    mutate(
        DEID = as.character(DEID),
        Trilaciclib = ifelse(grepl("Placebo", ARM), "Placebo", "Trilaciclib")
    ) %>%
    dplyr::select(
        OurID, DEID, whichdraw,
        SEX, RACE, ETHNIC, Age, Trilaciclib
    ) %>%
    dplyr::rename(
        DeID = DEID,
        Sex = SEX,
        Race = RACE,
        Ethnicity = ETHNIC,
        SampleDeID = OurID
    ) %>%
    mutate(
        Sex = ifelse(Sex == "F", "Female", "Male"),
        Race = case_when(
            Race == "AMERICAN INDIAN OR ALASKA NATIVE" ~ "American Indian or Alaska Native",
            Race == "WHITE" ~ "White",
            Race == "ASIAN" ~ "Asian"
        ),
        Ethnicity = case_when(
            Ethnicity == "HISPANIC OR LATINO" ~ "Hispanic or Latino",
            Ethnicity == "NOT HISPANIC OR LATINO" ~ "Not Hispanic or Latino",
            Ethnicity == "UNKNOWN" ~ "Unknown"
        ),
        SmokingStatus = "Unknown",
        Cohort = "SCLC",
        whichdraw = ifelse(whichdraw == 1, "PreTx", "PostTx")
    ),
g1_info_extra %>% 
    dplyr::select(`Originating ID`, `Subject #`, `Time Point`, AGE, Trilaciclib) %>%
    dplyr::rename(
        SampleDeID = `Originating ID`,
        subject = `Subject #`,
        timepoint = `Time Point`,
        Age = AGE
    ) %>%
    mutate(
        SampleDeID = as.character(SampleDeID),
        whichdraw = case_when(
            grepl("Cycle 1", timepoint) ~ "PreTx",
            grepl("C1D1", timepoint) ~ "PreTx",
            grepl("Cycle 5", timepoint) ~ "PostTx",
            grepl("C5D1", timepoint) ~ "PostTx",
            grepl("D90", timepoint) ~ "PostPostMaintenance",
            grepl("Day 90", timepoint) ~ "PostPostMaintenance"
        ),
        DeID = gsub("-", "", subject),
        Sex = NA_character_,
        Race = NA_character_,
        Ethnicity = NA_character_,
        SmokingStatus = NA_character_,
        Trilaciclib = ifelse(Trilaciclib == 0, "Placebo", "Trilaciclib"),
        Cohort = "SCLC"
    ) %>%
    dplyr::select(-c(timepoint, subject))
)

tnbc_demo <- g1_info_tnbc %>%
    dplyr::select(`Originating ID`, `Subject ID`, whereincycle, Trilaciclib, Cycles_Received, Cycle_1, Cycle_7, Post, Post_60, Months_From, Age_Start, Gender) %>%
    mutate(DeID = gsub("-", "", `Subject ID`)) %>%
    dplyr::rename(
        SampleDeID = `Originating ID`,
        Sex = Gender
    ) %>%
    mutate(
        Sex = ifelse(Sex == "F", "Female", "Male"),
        Age = case_when(
            whereincycle == "pre" ~ Age_Start,
            whereincycle == "during" ~ Age_Start + (4.8/12),
            whereincycle == "post" & Post_60 == "no" ~ Age_Start + (Months_From/12),
            whereincycle == "post" & Post_60 == "yes" ~ Age_Start + ((Months_From-2)/12),   # Since post_60 is 60 days after post (2 months) post should be ~postpost - 2 months
            whereincycle == "postpost" ~ Age_Start + (Months_From/12)
        ),
        whichdraw = case_when(
            whereincycle == "pre" ~ "PreTx",
            whereincycle == "during" ~ "MidTx",
            whereincycle == "post" ~ "PostTx",
            whereincycle == "postpost" ~ "PostPostMaintenance"
        ),
        Trilaciclib = ifelse(Trilaciclib == "trila", "Trilaciclib", "Placebo"),
        Race = NA_character_,
        Ethnicity = NA_character_,
        SmokingStatus = NA_character_,
        Cohort = "TNBC"
    ) %>%
    dplyr::select(-c(`Subject ID`, Cycles_Received, Cycle_1, Cycle_7, Post, Post_60, Months_From, Age_Start, whereincycle))

crc_demo <- g1_crc_bisi_sample_info %>%
    filter(fromwhere == "G1 Therapeutics - CRC") %>%
    dplyr::select(
        sampleID, patientID, whereincycle, Trilaciclib
    ) %>%
    dplyr::rename(
        SampleDeID = sampleID,
        DeID = patientID
    ) %>%
    mutate(
        whichdraw = case_when(
            whereincycle == "MAINTENANCE C1D1" ~ "PostTx",
            whereincycle == "C1D1" ~ "PreTx",
            whereincycle == "C5D1" ~ "MidTx"
        ),
        Age = NA_real_,
        Sex = NA_character_,
        Race = NA_character_,
        Ethnicity = NA_character_,
        SmokingStatus = NA_character_,
        Cohort = "CRC"
    ) %>%
    dplyr::select(-whereincycle)

#rbind(untreated_demo, sclc_demo, tnbc_demo, crc_demo) %>% fwrite("~/Documents/irenaeus/CDK46_CH/data/g1_demo.csv")

g1_all_kelly <- rbind(
sclc_tp %>% 
    mutate(
        Cohort = "SCLC",
        duringAge = NA,
        compare_group_during = NA,
        change_in_days_during = NA,
        growth_rate_during = NA
    ) %>%
    dplyr::select(
        deid_key, key, DeID, Gene,
        all_time_points, prepreVAF, preVAF, duringVAF, postVAF, postpostVAF,
        prepreVAF_ignore_noise, preVAF_ignore_noise, duringVAF_ignore_noise, postVAF_ignore_noise, postpostVAF_ignore_noise,
        rescued, missing, called, at_limit, preAge, duringAge, postAge, in_cohort,
        DDR, Splice, DTA, Gene_Class, change_in_days, growth_rate, change_in_VAF, Trilaciclib,
        compare_group, compare_group_during, change_in_days_during, growth_rate_during
  ) %>% mutate(
        key = gsub(" |>", ":", key)
  ),
tnbc_tp %>%
    filter(Trilaciclib != "Untreated") %>%
    mutate(
        Cohort = "TNBC"
    ) %>%
    dplyr::select(
        deid_key, key, DeID, Gene,
        all_time_points, prepreVAF, preVAF, duringVAF, postVAF, postpostVAF,
        prepreVAF_ignore_noise, preVAF_ignore_noise, duringVAF_ignore_noise, postVAF_ignore_noise, postpostVAF_ignore_noise,
        rescued, missing, called, at_limit, preAge, duringAge, postAge, in_cohort,
        DDR, Splice, DTA, Gene_Class, change_in_days, growth_rate, change_in_VAF, Trilaciclib,
        compare_group, compare_group_during, change_in_days_during, growth_rate_during
  ),
crc_tp %>%
    filter(Trilaciclib != "Untreated") %>%
    mutate(
        Cohort = "CRC"
    ) %>%
    dplyr::select(
        deid_key, key, DeID, Gene,
        all_time_points, prepreVAF, preVAF, duringVAF, postVAF, postpostVAF,
        prepreVAF_ignore_noise, preVAF_ignore_noise, duringVAF_ignore_noise, postVAF_ignore_noise, postpostVAF_ignore_noise,
        rescued, missing, called, at_limit, preAge, duringAge, postAge, in_cohort,
        DDR, Splice, DTA, Gene_Class, change_in_days, growth_rate, change_in_VAF, Trilaciclib,
        compare_group, compare_group_during, change_in_days_during, growth_rate_during
  )
)

g1_all_kelly %>%
    left_join(
        g1_df %>%
            dplyr::select(
                key, DeID, gene_aachange, oncoKB, oncoKB_reviewed, starts_with("clinvar"),
                n.loci.vep, source.totals.loci, n.loci.truncating.vep, source.totals.loci.truncating, n.HGVSp, source.totals.p, n.HGVSc, source.totals.c, CosmicCount, heme_cosmic_count, myeloid_cosmic_count
            )
    ) %>% fwrite("~/Documents/irenaeus/CDK46_CH/data/g1_all_kelly.timepoints.csv")
   

#g1_demo <- fread("~/Documents/irenaeus/CDK46_CH/data/g1_demo.csv")
#g1_demog <- fread("~/Documents/irenaeus/CDK46_CH/data/G1_demog.csv")

rbind(
g1_demo %>%
    mutate(
        SmokingStatus = ifelse(Cohort != "Untreated", NA_character_, SmokingStatus)
    ) %>%
    left_join(
        g1_demog %>%
            mutate(
                Cohort = case_when(
                    Protocol == "G1T28-02" ~ "SCLC",
                    Protocol == "G1T28-04" ~ "TNBC",
                    Protocol == "G1T28-05" ~ "SCLC",
                    Protocol == "G1T28-207" ~ "CRC"
                ),
                DeID = ifelse(
                    Cohort == "CRC",
                    Subject,
                    gsub("-", "", Subject)
                )
            )
    ) %>%
    filter(Cohort != "Untreated") %>%
    mutate(
        Sex = ifelse(Sex == "", SEX, Sex),
        Sex = case_when(
            Sex == "F" ~ "Female",
            Sex == "M" ~ "Male",
            TRUE ~ Sex
        ),
        Race = ifelse(Race == "", RACE, Race),
        Race = case_when(
            Race == "AMERICAN INDIAN OR ALASKA NATIVE" ~ "American Indian or Alaska Native",
            Race == "WHITE" ~ "White",
            Race == "ASIAN" ~ "Asian",
            Race == "BLACK OR AFRICAN AMERICAN" ~ "Black",
            Race == "OTHER" ~ "Other",
            Race == "UNKNOWN" ~ "Unknown",
            Race == "NOT REPORTED" ~ "Unknown",
            TRUE ~ Race
        ),
        Ethnicity = ifelse(Ethnicity == "", ETHNIC, Ethnicity),
        Ethnicity = case_when(
            Ethnicity == "HISPANIC OR LATINO" ~ "Hispanic or Latino",
            Ethnicity == "NOT HISPANIC OR LATINO" ~ "Not Hispanic or Latino",
            Ethnicity == "UNKNOWN" ~ "Unknown",
            Ethnicity == "NOT REPORTED" ~ "Unknown",
            TRUE ~ Ethnicity
        ),
        SmokingStatus = ifelse(is.na(SmokingStatus), SMKSTAT, SmokingStatus),
        SmokingStatus = case_when(
            SmokingStatus == "Current Smoker" ~ "Current",
            SmokingStatus == "Never Smoked" ~ "Never",
            SmokingStatus == "Former Smoker" ~ "Former",
            SmokingStatus == "" ~ "Unknown",
            TRUE ~ SmokingStatus
        ),
        Age = ifelse(is.na(Age), AGE, Age)
    ) %>%
    dplyr::select(-Subject, -AGE, -AGEU, -SEX, -RACE, -ETHNIC, -SMKSTAT, -Protocol),
    g1_demo %>% filter(Cohort == "Untreated")
) %>% fwrite("~/Documents/irenaeus/CDK46_CH/data/g1_demo.csv")

g1_demo <- fread("~/Documents/irenaeus/CDK46_CH/data/g1_demo.csv")
g1_lab <- fread("~/Documents/irenaeus/CDK46_CH/data/g1_labresults.csv")

g1_lab <- g1_lab %>%
    mutate(
        whichdraw = case_when(
            protocol == "G1T28-02" & visit == "Cycle 1 Day 1" ~ "PreTx",
            protocol == "G1T28-02" & visit == "Cycle 5 Day 1" ~ "PostTx",
            protocol == "G1T28-02" & visit == "Post Treatment" ~ "PostTx",

            protocol == "G1T28-05" & visit == "Induction Day 01 Cycle 01" ~ "PreTx",
            protocol == "G1T28-05" & visit == "Maintenance Cycle 05" ~ "PostTx",
            protocol == "G1T28-05" & visit == "90 Days Post-treatment visit" ~ "PostPostMaintenance",
            protocol == "G1T28-05" & visit == "Maintenance Cycle 03" ~ "PostTx",
            protocol == "G1T28-05" & visit == "Maintenance Cycle 04" ~ "PostTx",
            protocol == "G1T28-05" & visit == "Maintenance Cycle 01" ~ "PostTx - Remove",
            protocol == "G1T28-05" & visit == "Screening" ~ "PreTx",

            protocol == "G1T28-04" & visit == "Cycle 1 Day 1" ~ "PreTx",
            protocol == "G1T28-04" & visit == "Cycle 7 Day 1" ~ "MidTx",
            protocol == "G1T28-04" & visit == "Post-Treatment Visit 1" ~ "PostTx",
            protocol == "G1T28-04" & visit == "Post-Treatment Visit 2 (+60 days)" ~ "PostMaintenance",
            protocol == "G1T28-04" & visit == "UNSCHEDULED 1315.03" ~ "PostTx",

            protocol == "G1T28-207" & visit == "Induction Cycle 1 Day 1" ~ "PreTx",
            protocol == "G1T28-207" & visit == "Maintenance Cycle 1 Day 1" ~ "PostTx",
            protocol == "G1T28-207" & visit == "Induction Cycle 5 Day 1" ~ "MidTx",
            protocol == "G1T28-207" & visit == "Induction Cycle 6 Day 1" ~ "MidTx",
            protocol == "G1T28-207" & visit == "Last Induction Cycle Day 15" ~ "PostTx",
            protocol == "G1T28-207" & visit == "Screening" ~ "PreTx",

            TRUE ~ "MISSING"
        ),
        Cohort = case_when(
            protocol == "G1T28-02" ~ "SCLC",
            protocol == "G1T28-05" ~ "SCLC",
            protocol == "G1T28-04" ~ "TNBC",
            protocol == "G1T28-207" ~ "CRC"
        )
    )

g1_lab %>%
    pivot_wider(
        id_cols = c(protocol, subject, whichdraw),
        names_from = Testcode,
        values_from = Labvalue
    )
    
g1_lab %>%
    filter(
        case_when(
            Cohort == "SCLC"
        )
    )
    pivot_wider(
        id_cols = c(protocol, subject),
        names_from = c(Testcode, whichdraw),
        values_from = Labvalue
    )
    mutate(
        change_in_WBC = case_when(
            protocol == "G1T28-02" | protocol == "G1T28-05" ~ WBC_PostTx - WBC_PreTx,
            protocol == "G1T28-04" ~ WBC_MidTx - WBC_PreTx,
            protocol == "G1T28-207" ~ WBC_PostTx - WBC_PreTx
        ),
        change_in_NEUT = case_when(
            protocol == "G1T28-02" | protocol == "G1T28-05" ~ NEUT_PostTx - NEUT_PreTx,
            protocol == "G1T28-04" ~ NEUT_MidTx - NEUT_PreTx,
            protocol == "G1T28-207" ~ NEUT_PostTx - NEUT_PreTx
        )
    ) %>%
    filter(!is.na(change_in_WBC) | !is.na(change_in_NEUT))


g1_df %>%
    left_join(
        g1_lab_ %>% dplyr::select(subject, whichdraw, WBC, NEUT),
        by = c("DeID" = "subject", "whichdraw" = "whichdraw")
    ) %>%
    filter(Gene %in% c("TP53", "PPM1D", "CHEK2")) %>%
    filter(Cohort != "Untreated") %>%
    filter(Gene == "CHEK2") %>%
    glm(
        average_AF ~ WBC + NEUT,
        data = .,
        family = gaussian
    ) %>%
    sjPlot::get_model_data(type = "est")

g1_df %>%
    left_join(
        g1_lab_ %>% dplyr::select(subject, whichdraw, WBC, NEUT),
        by = c("DeID" = "subject", "whichdraw" = "whichdraw")
    ) %>%
    filter(Gene %in% c("TP53", "PPM1D", "CHEK2")) %>%
    filter(Cohort != "Untreated") %>%
    ggplot(
        aes(
            x = NEUT,
            y = average_AF,
            group = Gene,
            color = Gene
        )
    ) +
    geom_point() +
    geom_smooth(method = "lm") +
    theme_bw()

g1_tp %>% 
    left_join(
        g1_lab_ %>% dplyr::select(subject, change_in_WBC, change_in_NEUT),
        by = c("DeID" = "subject")
    ) %>%
    filter(Gene %in% c("TP53", "PPM1D", "CHEK2")) %>%
    filter(Cohort != "Untreated") %>%
    filter(growth_rate != "Inf" & growth_rate != "-Inf") %>%
    ggplot(
        aes(
            x = change_in_WBC,
            y = change_in_VAF,
            group = Gene,
            color = Gene
        )
    ) + 
    geom_point() +
    geom_smooth(method = "lm") +
    theme_bw()
    
g1_tp %>% 
    left_join(
        g1_lab_ %>% dplyr::select(subject, change_in_WBC, change_in_NEUT),
        by = c("DeID" = "subject")
    ) %>%
    filter(Gene %in% c("TP53", "PPM1D", "CHEK2")) %>%
    filter(Cohort != "Untreated") %>%
    filter(growth_rate != "Inf" & growth_rate != "-Inf") %>%
    lm(
        formula = growth_rate ~ change_in_WBC + change_in_NEUT
    ) %>% sjPlot::get_model_data(type = "est")



# g1_demo %>% 
#     filter(Cohort != "Untreated") %>% 
#     filter(ifelse(Cohort == "TNBC", whichdraw == "PreTx" | whichdraw == "MidTx", whichdraw == "PreTx" | whichdraw == "PostTx")) %>%
#     filter(DeID %in% c(g1_demo %>% 
#                         filter(Cohort != "Untreated") %>% 
#                         filter(ifelse(Cohort == "TNBC", whichdraw == "PreTx" | whichdraw == "MidTx", whichdraw == "PreTx" | whichdraw == "PostTx")) %>%
#                         group_by(Cohort, DeID) %>%
#                         summarise(n = n()) %>%
#                         filter(n > 1) %>%
#                         pull(DeID))
#     ) %>%
#     left_join(
#         sample_mapping_all_g1
#     ) %>% fwrite("G1_Samples_Request.csv")

sample_mapping_all_g1 <- rbind(
    tnbc_samplemap %>% 
            dplyr::select(`ESP ID`, `Library Name`) %>%
            mutate(`Library Name` = as.character(`Library Name`)) %>%
            left_join(
                g1_info_tnbc, #%>%
                #     dplyr::select(SampleDeID, site, DeID, whereincycle, Trilaciclib),
                by = c("Library Name" = "SampleDeID")
            ) %>%
            mutate(
                `Sample Name` = `Library Name`,
                SampleDeID = `Sample Name`
            ) %>%
            dplyr::select(`ESP ID`, `Sample Name`, DeID, SampleDeID, site, whereincycle),
    g1_extra_samplemap %>% 
    dplyr::select(`GTAC ESP ID`, `Sample Name`) %>%
    filter(
        !grepl("U10", `Sample Name`),
        !grepl("MCIA", `Sample Name`),
        grepl("-", `Sample Name`)
    ) %>%
    mutate(
        `Sample Name` = as.character(`Sample Name`),
        SampleDeID = str_split(`Sample Name`, "_") %>% sapply("[", 1),
        DeID = str_split(`Sample Name`, "_") %>% sapply("[", 2) %>% gsub("-", "", .),
        `ESP ID` = `GTAC ESP ID`
    ) %>%
    left_join(
        g1_info_extra %>%
            dplyr::select(sample_id, site, DeID, timepoint) %>%
            mutate(
                whereincycle = case_when(
                    grepl("Cycle 1", timepoint) ~ "pre",
                    grepl("C1D1", timepoint) ~ "pre",
                    grepl("Cycle 5", timepoint) ~ "post",
                    grepl("C5D1", timepoint) ~ "post",
                    grepl("D90", timepoint) ~ "postpost",
                    grepl("Day 90", timepoint) ~ "postpost"
                )
            ),
        by = c("SampleDeID" = "sample_id", "DeID" = "DeID")
    ) %>%
    dplyr::select(`ESP ID`, `Sample Name`, DeID, SampleDeID, site, whereincycle),
g1_crc_bisi_sample_info %>%
dplyr::select(MGI_ID, sampleID, patientID, Site, whereincycle) %>%
dplyr::rename(`Sample Name` = sampleID, `ESP ID` = MGI_ID, DeID = patientID, site = Site) %>%
mutate(SampleDeID = `Sample Name`)
) #%>% fwrite("G1_Samples_Request.csv")