getfactorvariables <- function(df, vars) {
    classes <- lapply(c2l(vars), function(x) class(df[[x]]))
    classes[classes == "factor"] %>% names
}

characteristics.names <- function(onlyvars = FALSE) {
    if (onlyvars == TRUE)
        return(c("BL_AGE",
                 "SEX",
                 "BMI",
                 "SYSTM",
                 "PREVAL_AF",
                 "INCIDENT_AF",
                 "PREVAL_DIAB",
                 "CURR_SMOKE",
                 "BP_TREAT",
                 "KOL"))
    list("BL_AGE" = "Age, y (SD)",
         "SEX" = "Female, N (%)",
         "BMI" = "BMI, kg/m² (SD)",
         "SYSTM" = "Systolic BP, mmHg (SD)",
         "PREVAL_DIAB" = "Diabetes mellitus, N (%)", 
         "PREVAL_AF" = "Prevalent AF, N (%)",
         "INCIDENT_AF" = "Incident AF, N (%)", 
         "CURR_SMOKE" = "Current smoker, N (%)",
         "BP_TREAT" = "Antihypertensive medication, N (%)",
         "KOL" = "  Total cholesterol, mmol/l (SD)")
}

characteristicsTable <- function(dset, strata) {
    nostrata <- missing(strata)
    characteristics <- tableone::CreateTableOne(
                                     strata = strata,
                                     data = dset,
                                     vars = characteristics.names(TRUE),
                                     factorVars = getfactorvariables(dset,
                                                                     characteristics.names(TRUE)),
                                     test = !missing(strata))
    print(characteristics,
                      exact = "stage",
                      quote = FALSE,
                      noSpaces = TRUE,
                      printToggle = FALSE,
                      digits = 1,
                      pDigits = 3,
                      contDigits=1)  %>%
        as.data.frame %>%
        tibble::rownames_to_column(var = "rowname") %>%
        format(justify = "left", trim = TRUE) %>%
        mutate(rowname = characteristics.names()[gsub("^ *([A-Za-z_0-9]+).*", "\\1", rowname)]) %>%
        { if (!nostrata) select(., -test) else . }
}

characteristics <- function(dset, tableone.names, tableone.factors, extras = list()) {
    title <- "Characteristics"
    overall <- paste0("Cases, n=", dim(dset)[1])
    tableobject <- tableone::CreateTableOne(data = dset, vars = names(tableone.names), factorVars = tableone.factors)
    tablecsv <- print(tableobject,
                      exact = "stage",
                      quote = FALSE,
                      noSpaces = TRUE,
                      printToggle = FALSE,
                      digits = 1,
                      pDigits = 3,
                      contDigits=1)

    tableone.fullnames <- c(tableone.names, extras)
    tablecsv %>%
        as.data.frame %>%
        tibble::rownames_to_column(var = "rowname") %>%
        dplyr::filter(row_number() > 1) %>%
        format(justify = "left", trim = TRUE) %>%
        rowwise() %>%
        mutate(id = gsub("^ *([A-Za-z_0-9]+).*", "\\1", rowname)) %>%
        mutate(present = id %in%  names(tableone.fullnames)) %>%
        mutate(rowname = ifelse(present == TRUE, tableone.fullnames[[id]], rowname)) %>%
        select(rowname, Overall)
}

typologyformatter <- function(data, font = 12, typology, left = c(1), hleft = c(1)) {
  flex <- flextable(data = data) %>%
    flextable::theme_booktabs() %>%
    flextable::border(border = fp_border(width=0), part="body") %>%
    flextable::border(border = fp_border(width=0), part="header") %>%
    flextable::border(part="header", border.bottom = fp_border(width=1))

  if (!missing(typology)) {
      flex <- flex %>%
          set_header_df(mapping = typology, key = "col_keys") %>%
          merge_h(part = "header") %>%
          flextable::border(part="header", border.bottom = fp_border(width=1))
      if (missing(hleft)) {
          hleft <- c(2)
      }
  }

  flex %>%
      flextable::border(i = nrow(data), part="body", border.bottom = fp_border(width=1)) %>%
      flextable::bold(bold = FALSE, part = "header") %>%
      flextable::bold(bold = FALSE, part = "body") %>%
      flextable::fontsize(size = font, part = "header") %>%
      flextable::fontsize(size = font, part = "body") %>%
      flextable::align(align = "center", part = "all") %>%
      flextable::align(align = "left", part = "header", j = left, i = hleft) %>%
      flextable::align(align = "left", part = "body", j = left)
}

myspread <- function(ret, list = c2l("lfc_se", "p.value"), term  = "Feature", key = "Name") {
    lapply(list, function(x) ret %>%
                             select(term, key, x) %>%
                             spread(key, x) %>%
                             rename_at(vars(-term), funs(paste0(., "_", x)))) %>%
        Reduce(function(...) full_join(..., by = term), .) %>%
        dplyr::select(term, noquote(order(colnames(.))))
}