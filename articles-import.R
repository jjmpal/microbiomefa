import_filter_data <- function(species = FALSE, genus = FALSE, vars, dropna = TRUE) {
    stopifnot(xor(species, genus))
    fphyloseq <- case_when(species ~ "data/phfinrisk_species_all_drop50k_2018-11-16.RDs",
                           genus ~ "data/phfinrisk_genus_all_drop50k_2018-11-16.RDs")
    pseq.full <- readRDS(fphyloseq)
    pseq.meta <- filter.phenotype.data(pseq = pseq.full, vars = c(vars), dropna = dropna)
    phyloseq::sample_data(pseq.full) <- phyloseq::sample_data(pseq.meta)
    return(pseq.full)
}

filter.phenotype.data <-
    function(pseq,
             vars,
             fextra = "pheno/2015_60_Salomaa_Jain_dataFR02_FU17_2020-11-15.txt",
             extra = c("Barcode",
                       "AF_AGE",
                       "AF_AGEDIFF",
                       "INCIDENT_AF",
                       "PREVAL_AF",
                       "PREVAL_HFAIL_STRICT",
                       "ALKI2_FR02"),
             allowna = c("NA."),
             dropna = TRUE) {
        meta(pseq) %>%
            tibble::rownames_to_column(var = "rowname") %>%
            select_if(function(x) all(!x %in% extra)) %>%
            left_join(., read_tsv(fextra, col_types = cols(
                                           .default = col_double(),
                                           CDT = col_character(),
                                           Sample_ID = col_character(),
                                           Barcode = col_character(),
                                           WESID = col_character(),
                                           WGSID = col_character(),
                                           BATCH = col_character(),
                                           FID = col_character(),
                                           CRP = col_character(),
                                           K_TPKS = col_character(),
                                           K_VKS = col_character(),
                                           K_M1 = col_character(),
                                           K_M2 = col_character(),
                                           APOE_BATCH = col_character())) %>%
                      select(one_of(extra)),
                      by = c("rowname" = "Barcode")) %>%
            dplyr::mutate(Sample_ID = rowname,
                          BP_TREAT = factor(ifelse(BL_USE_RX_C03  == 1 |
                                                   BL_USE_RX_C07 == 1 |
                                                   BL_USE_RX_C08 == 1  |
                                                   BL_USE_RX_C09 == 1, 1, 0)),
                          SEX = factor(ifelse(MEN == "Female", 1, 0)),
                          HEARTFAILURE = as.factor(PREVAL_HFAIL_STRICT),
                          ETOH = log(ALKI2_FR02 + 1)) %>%
            dplyr::mutate_at(vars(starts_with("PREVAL_"), starts_with("INCIDENT_")), as.factor) %>%
            tibble::column_to_rownames(var = "rowname") %>%
            select(one_of(vars)) %>%
            { if (dropna) tidyr::drop_na(., vars %difference% allowna) else . }
}
