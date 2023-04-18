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
          add_HFC_to_dataframe %>%
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
          { if ("Q57X" %in% vars) mutate(., Q57X = factor(Q57X, ordered = FALSE)) else . } %>%
          tibble::column_to_rownames(var = "rowname") %>%
            select(one_of(vars)) %>%
            { if (dropna) tidyr::drop_na(., vars %difference% allowna) else . }
}

add_HFC_to_dataframe <- function(x) {
  food_categories <- c('Ruis-/nakkileipa' = 'KY100_1',
                       'Hiiva-/graham-/sekaleipa' = 'KY100_2',
                       'Ranskanleipa' = 'KY100_3',
                       'Puurot' = 'KY100_5',
                       'Mysli' = 'KY100_6',
                       'Riisi/makaroni' = 'KY100_7',
                       'Tuoreet vihannekset/juurekset' = 'KY100_14',
                       'Keitetyt kasvikset/palkokasvit' = 'KY100_15',
                       'Kasvisruoka' = 'KY100_16',
                       'Hedelmat' = 'KY100_17',
                       'Tuoreet/pakastetut marjat' = 'KY100_18',
                       'Hedelma-/marjatuoremehu' = 'KY100_19',
                       'Kala/kalaruoka' = 'KY100_20',
                       'Broileri/kanaruoka' = 'KY100_21',
                       'Liharuoka' = 'KY100_22',
                       'Makkarat/nakit' = 'KY100_23',
                       'Leikkelemakkarat' = 'KY100_24',
                       'Lihaleikkeleet' = 'KY100_25',
                       'Viili/jogurtti' = 'KY100_8',
                       'Vaharasvaiset juustot' = 'KY100_9',
                       'Piirakat/pasteijat' = 'FR02_101_1',
                       'Salaattikastike/Oljy' = 'FR02_101_3',
                       'Pahkinat' = 'FR02_101_7',
                       'Siemenet' = 'FR02_101_8',
                       'Maitolaatu' = 'FR02_112',
                       'Ruoanvalmistusrasva' = 'FR02_104',
                       'Pizza' = 'FR02_101_4',
                       'Hampurilaiset' = 'FR02_101_5',
                       'Suolaiset naposteltavat' = 'KY100_31',
                       'Sokeroidut, virvoitusjuomat' = 'KY100_29',
                       'Pikaruoka' = 'KY100_37',
                       'Jaatelot/vanukkaat/rahkat' = 'FR02_101_2',
                       'Suklaa' = 'KY100_27',
                       'Makeiset' = 'KY100_28',
                       'Makea, kahvileipa' = 'KY100_4',
                       'Muut juustot' = 'KY100_10',
                       'Kananmuna' = 'KY100_26')

  dietary_frequency <- list(
    rarely_daily =  c("Pizza",
                      "Hampurilaiset",
                      "Suolaiset naposteltavat",
                      "Jaatelot/vanukkaat/rahkat"),
    often_daily = c("Tuoreet vihannekset/juurekset",
                    "Hedelmat",
                    "Vaharasvaiset juustot",
                    "Muut juustot",
                    "Ruis-/nakkileipa",
                    "Hiiva-/graham-/sekaleipa"),
    red_meat =  c("Liharuoka",
                  "Makkarat/nakit",
                  "Leikkelemakkarat",
                  "Lihaleikkeleet"),
    others = c("Ranskanleipa",
               "Puurot",
               "Mysli",
               "Riisi/makaroni",
               "Keitetyt kasvikset/palkokasvit",
               "Kasvisruoka",
               "Tuoreet/pakastetut marjat",
               "Hedelma-/marjatuoremehu",
               "Kala/kalaruoka",
               "Broileri/kanaruoka",
               "Viili/jogurtti",
               "Piirakat/pasteijat",
               "Salaattikastike/Oljy",
               "Pahkinat",
               "Siemenet",
               "Maitolaatu",
               "Ruoanvalmistusrasva",
               "Sokeroidut, virvoitusjuomat",
               "Pikaruoka",
               "Suklaa",
               "Makeiset",
               "Makea, kahvileipa",
               "Kananmuna"))

  food_groups_for_diet_scores <- c(BREAD = c("Ruis-/nakkileipa",
                                             "Hiiva-/graham-/sekaleipa"),
                                   VEGETABLES = c("Tuoreet vihannekset/juurekset",
                                                  "Keitetyt kasvikset/palkokasvit",
                                                  "Kasvisruoka"),
                                   FRUIT = c("Hedelmat"),
                                   BERRIES = c("Tuoreet/pakastetut marjat"),
                                   JUICES = c("Hedelma-/marjatuoremehu"),
                                   FISH = c("Kala/kalaruoka"),
                                   POULTRY = c("Broileri/kanaruoka"),
                                   LOW_FAT_CHEESE = c("Vaharasvaiset juustot"),
                                   DRESSINGS_OILS = c("Salaattikastike/Oljy"),
                                   NUTS_SEEDS = c("Pahkinat", "Siemenet"))


  food_score_factory <- function(high = 30, reverse = FALSE) {
    y = c(0.5, 1.5, 4.3, 8.6, 21.5, high) %>%
      { if(reverse) rev(.) else .}
    function(x) {
      case_when(x == 1 ~ y[1], x == 2 ~ y[2], x == 3 ~ y[3], x == 4 ~ y[4], x == 5 ~ y[5], x == 6 ~ y[6], TRUE ~ 0)
    }
  }

  df_HFC <- x %>%
    dplyr::mutate_at(vars(food_categories), .funs = as.factor) %>%
    dplyr::mutate_at(vars(food_categories[dietary_frequency$rarely_daily]), food_score_factory()) %>%
    dplyr::mutate_at(vars(food_categories[dietary_frequency$often_daily]), food_score_factory(high = 60)) %>%
    dplyr::mutate_at(vars(food_categories[dietary_frequency$red_meat]), food_score_factory(high = 45, reverse = TRUE)) %>%
    dplyr::mutate_at(vars(food_categories[dietary_frequency$others]), food_score_factory(high = 45)) %>%
    dplyr::mutate(HFC := eval(parse(text = paste(food_categories[food_groups_for_diet_scores], collapse = "+")))) %>%
    dplyr::select(Barcode, HFC)

  dplyr::left_join(x, df_HFC, by = "Barcode")
}
