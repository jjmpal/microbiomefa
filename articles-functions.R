`%difference%` <- function(a, b) {
    a[!a %in% b]
}

`%intersect%` <- function(a, b) {
    intersect(a, b)
}


`%union%` <- function(a, b) {
    c(a, b)
}

replace.brackets <- function (genus) {
    gsub('(.*) \\(+(.+)\\)', '\\1.\\2', genus)
}

renametaxa <- function (names) {
    suffix <- ifelse(grepl("BacteriaPlasmid", names), "*", "")
    base <- ifelse(grepl("_", names),
           gsub('(.).*_(.*) \\(+(.+)\\)', '\\1. \\2', names),
           gsub('(.*) \\(+(.+)\\)', '\\1', names))
    paste0(base, suffix)
}

mycbind <- function(list) {
    max.length <- c(list) %>% lapply(length) %>% unlist %>% max
    lapply(list, function(x) {
        length(x) <- max.length
        x}) %>%
        do.call(cbind, .)
}

list.partition <- function(list, parts = 3) {
    partition <- floor(parts*seq(0, length(list)-1)/length(list))
    split(list, partition)
}

rename.genus  <- function (genus, markword = "Plasmid", mark = "*") {
    name <- gsub('(.*) \\(.*\\)', '\\1', genus)
    star <- ifelse(grepl(markword, genus), mark, "")
    paste0(gsub("_", " ", name), star)
}

myscale <- function(x){
  (x - mean(x, na.rm=TRUE)) / sd(x, na.rm=TRUE)
}


meta.merge.alphadiversity <- function(pseq, index = "shannon") {
    alphadiversity  <- microbiome::alpha(pseq, index = index)
    base::merge(meta(pseq), alphadiversity, by=0, all=TRUE) %>%
        dplyr::rename(Sample_ID = Row.names) %>%
        dplyr::mutate(diversity_shannon = scale(diversity_shannon))
}

calculate.beta.matrix <- function(pseq) {
    compositional.profile <- microbiome::transform(pseq, 'compositional')
    otu <- microbiome::abundances(compositional.profile)
    bray.dist.m <- vegan::vegdist(t(otu), method="bray")
    dist.matrix <- as.matrix(bray.dist.m)
    attr(dist.matrix, "method") <- "bray"
    dist.matrix
}

calculateglm <- function(dset,
                         responses,
                         min_n_for_continuous = 10,
                         covariates = c(),
                         filterstr = ".") {
    stopifnot(!missing(dset), !missing(responses))
    glmlist <- lapply(responses, function(response) {
        is.logistic <- length(unique(pull(dset, response))) < min_n_for_continuous
        fo.family <- ifelse(is.logistic, stats::binomial, stats::gaussian)
        fo <- sprintf("%s ~ %s", response, paste0(covariates, collapse = "+"))
        stats::glm(formula = as.formula(fo), family = fo.family, data = dset) %>%
            broom::tidy(conf.int = TRUE, conf.level = 0.95, exponentiate = is.logistic) %>%
                dplyr::filter(grepl(filterstr, term)) %>%
                dplyr::mutate(response = response, fo = fo)
    })
    data.table::rbindlist(glmlist, id = "model_name") %>%
        as.data.frame
}


 merge.pheno.abu <- function(pseq,
                             core = FALSE,
                             core.detection = 0.001,
                             core.prevalence = 0.01) {
     pseq.comp <- microbiome::transform(pseq, "compositional")
     if (core) {
         pseq.comp  <- microbiome::core(pseq.comp, detection = core.detection, prevalence = core.prevalence)
    }
     pheno <- microbiome::meta(phyloseq::sample_data(pseq.comp))
     sample.ids <- rownames(pheno)
     abud.rel <- abundances(pseq.comp)
     rownames(abud.rel) <- sapply(rownames(abud.rel), replace.brackets)
     cbind(pheno, t(abud.rel))
}


legend_extractor <- function (tmp) {
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  tmp$grobs[[leg]]
}

ylabl_extractor <- function(tmp) {
  yl <- which(grepl('axis.title.y.left', sapply(tmp$grobs, function(x) x$name)))
  tmp$grobs[[yl]]
}

xlabb_extractor <- function(tmp) {
  xl <- which(grepl('axis.title.x.bottom', sapply(tmp$grobs, function(x) x$name)))
  tmp$grobs[[xl]]
}

vectortolist <- function(c) {
  l <- as.list(c)
  names(l) <- l
  l
}

getdescriptions <- function() {
    tribble(~Covariate, ~Category, ~Name, ~Desc,
            "PREVAL_AF", "Physical", "Atrial fibrillation", "Prevalent atrial fibrillation")
}

mygrep <- function(..., word, ignorecase = TRUE, complement = FALSE) {
    c(...)[xor(grepl(word, c(...), ignore.case = ignorecase), (complement == TRUE))]
}


myin <- function(x,y, complement = FALSE) {
    x[xor(x %in% y, complement)]
}

cc <- function(...) {
    c2l(c(...))
}


c2l <- function(...) {
    l <- as.list(c(...))
    names(l) <- c(...)
    l
}

pub.p <- function(p) {
    p <- as.numeric(p)
    case_when(p < 0.01 ~ sprintf("%.3e", p),
              TRUE ~ sprintf("%.3f", p))
}


firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}

mykable <- function(x, ...) {
  capture.output(x <- print(x))
  knitr::kable(x, ...)
}

myinstall.packages <- function(...) {
    list.of.packages <- c(...)
    new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
    if (length(new.packages) == 0) { return(TRUE) }
    for (package in new.packages) {
        message(sprintf("Installing: %s", package))
        myinstall.packages(gtools::getDependencies(package))
        install.packages(package)
    }
}


mydropna <- function(...) {
    c(...)[!is.na(c(...))]
}


diversities <- function(pseq, vars, responses, betadiversity) {
    stopifnot(!missing(pseq), !missing(vars), !missing(responses), !missing(betadiversity))
    lapply(c2l(names(vars)), function(var) {
        alphadiversity <- calculateglm(meta.merge.alphadiversity(pseq),
                                       responses = responses,
                                       covariates = c("diversity_shannon", vars[[var]]),
                                       filterstr = "shannon") %>%
            select(response, alpha.p = p.value, alpha.effect = estimate,
                   alpha.low = conf.low, alpha.high = conf.high)
        betadiversity <- betadiversity[[var]] %>%
            map_df(., ~as.data.frame(.x) %>%
                          dplyr::filter(term %in% responses) %>%
                          select(response, beta.R2 = R2, beta.p=`Pr(>F)`))
        list(beta = betadiversity, alpha = alphadiversity)
        full_join(alphadiversity, betadiversity, by = "response") %>%
            merge(getdescriptions() %>% select(Covariate, Name), by.x = "response", by.y = "Covariate")
    })
}

diversities.tidy <- function(diversity) {
    map_df(diversity, ~as.data.frame(.x), .id = "covariates") %>%
        mutate(alpha = sprintf("%.2f (%.2f-%.2f)", alpha.effect, alpha.low, alpha.high),
               alpha.p = pub.p(alpha.p),
               beta.R2 = sprintf("%.3f%%", 100*beta.R2),
               beta.p = pub.p(beta.p)) %>%
        select(Name, covariates, alpha, alpha.p, beta.R2, beta.p) 
}

calculate.betadiversity <- function(dset,
                                    matrix,
                                    vars,
                                    responses,
                                    npermutations = 999,
                                    maxcores = 20) {
    stopifnot(!missing(dset), !missing(matrix), !missing(vars), !missing(responses))
    lapply(vars, function(var) {
        mclapply(responses, function(response) {
            fo <- sprintf("matrix ~ %s + %s", paste(var, collapse = " + "), response)
            ad <- adonis(formula = as.formula(fo), data = dset, permutations = npermutations)
            ad$aov.tab %>%
                as.data.frame %>%
                tibble::rownames_to_column(var = "term") %>%
                dplyr::mutate(response = response)
        }, mc.cores = min(maxcores, length(responses)))
    })
}

prune_lactobacillus <- function(pseq, transform = "compositional", drop = TRUE) {
    microbiome::transform(pseq, transform = transform) %>%
        prune_taxa(c("g_Lactobacillus"), .) %>%
        psmelt %>%
        spread(OTU, Abundance) %>%
        dplyr::mutate(lacto.prevalent = ifelse(g_Lactobacillus > 0.1/100, 1, 0),
                      duna.present = ifelse(is.na(NA.), 0, 1)) %>%
        { if (drop == TRUE) dplyr::filter(., duna.present == 1) else . }
}

pseq_prevalence <- function(pseq, taxa = "Lactobacillus (Bacteria)", limit = 0.1/100, fun = mean) {
    pseq %>%
        microbiome::transform(transform = "compositional") %>%
        abundances %>%
        as.data.frame %>%
        .[taxa, ] %>%
        { ifelse(. > limit, 1, 0) } %>%
        fun
}

mytableone <- function(dset, variables) {
    { if (!missing(variables)) dplyr::select(dset, variables) else dset } %>%
        select(-Sample) %>% 
        gather(key,value) %>%
        group_by(key) %>%
        summarise_all(funs(N = n(),
                           Mean = mean(as.numeric(value), na.rm = TRUE),
                           ones = sum(value == 1, na.rm = TRUE),
                           sd = sd(as.numeric(value), na.rm = TRUE),
                           NAs = sum(is.na(value))))  %>%
        mutate_if(is.numeric, round, 3)
}

deseq.formula <- function(..., fo = "~ %s") {
    as.formula(sprintf(fo, paste0(c(...), collapse = "+")))
}

deseq.list <- function(var.CL, var.CL.min) {
    l <- lapply(c2l(var.CL), function(x) {
        unique(c(var.CL.min, x))
    })
    l[!duplicated(l)]
}

coretaxa <- function(pseq, detection = 0.1/100, prevalence = 1/100) {
    microbiome::transform(pseq, "compositional") %>%
        core(detection = detection, prevalence = prevalence) %>%
        taxa
}

pseqsubset <- function(pseq, coretaxa = NULL, transform = NULL, subset = TRUE) {
    participants <- pseq %>% meta %>% subset(PREVAL_AF == 0) %>% rownames
    pseq %>%
        { if (subset) prune_samples(participants, .) else . } %>%
        { if (!is.null(transform)) microbiome::transform(., transform) else . } %>%
        { if (!is.null(coretaxa)) prune_taxa(coretaxa, .) else . }
}

myDESeq <- function(pseq, vars, coretaxa, coreterm = ".", saltsubset = FALSE) {
    mydds.pseq <- pseqsubset(pseq, coretaxa = coretaxa, saltsubset = saltsubset)
    mydds.data <- phyloseq_to_deseq2(mydds.pseq, deseq.formula(vars))
    mydds.size <- estimateSizeFactors(mydds.data)
    mydds.dispersion <- taxa(mydds.pseq) %>%
        mygrep(word = coreterm) %>%
        mydds.size[.,] %>%
        estimateDispersions(fitType="parametric")
    nbinomWaldTest(mydds.dispersion)
}

loop.cox <- function(dset,
                     response,
                     loops,
                     covariates = c()) {
    dset.fixed <- dset %>%
        filter_at(vars(!!sym(paste0("PREVAL_", response))), ~ . == 0) %>%
        mutate_at(vars(!!sym(paste0("INCIDENT_", response))), ~as.numeric(as.character(.)))
    lapply(c2l(loops), function(loop, df) {
        fo <- sprintf("Surv(%s_AGEDIFF, INCIDENT_%s) ~ %s + %s",
                      response,
                      response,
                      loop,
                      paste(covariates, collapse = " + "))
        message(sprintf("%s for N=%s", fo, nrow(dset.fixed)))
        ret <- survival::coxph(as.formula(fo), ties = "breslow", data = df)    
        ret$call <- as.formula(fo)
        ret
    }, df = dset.fixed)
}

loop.results <- function(...,
                         filterstr = "Bacteria",
                         fdrbygroup = FALSE,
                         exponentiate = FALSE) {
    purrr::map_df(c(...), ~as.data.frame(tidy(.x, exponentiate = exponentiate))) %>%
        dplyr::filter(grepl(filterstr, term)) %>%
        dplyr::mutate(conf.low = estimate - qnorm(0.975) * std.error,
                      conf.high = estimate + qnorm(0.975) * std.error) %>%
        { if (fdrbygroup) group_by(., response) %>%
                              mutate(qval = p.adjust(p.value, method="BH")) %>%
                              ungroup
          else mutate(., qval = p.adjust(p.value, method="BH")) }
}

matdf <- function(matrix) {
    as.data.frame(matrix) %>% tibble::rownames_to_column(var = "rowname")
}

pseqlongform <- function(pseq, prunetaxa = NULL, subset = NULL, transform = "clr") {
    participants <- pseq %>% meta %>% subset(PREVAL_AF == 0) %>% rownames
    
    pseq.transprune <- pseq %>%
        { if (!is.null(prunetaxa)) prune_taxa(prunetaxa, .) else .} %>%
        { if (!is.null(subset)) prune_samples(participants, .) else .} %>%
        microbiome::transform(., transform = transform)

    pseq.abundances <- microbiome::abundances(pseq.transprune) %>%
        as_tibble(rownames = "taxa") %>%
        mutate(taxa = replace.brackets(taxa)) %>%
        gather(Sample_ID, abundance, -taxa) %>%
        spread(taxa, abundance)
    
    pseq.meta <- microbiome::meta(pseq.transprune) %>%
        as_tibble(rownames = "Sample_ID")
    
    full_join(pseq.meta, pseq.abundances, by = "Sample_ID")
}

mykeggget <- function(list, index = "DEFINITION") {
    lapply(list, function(x) paste(KEGGREST::keggGet(x)[[1]][[index]], collapse = ", ")) %>%
        unlist %>%
        gsub(" \\[.*\\]", "", .)
}

getmissing <- function(...) {
    rbind(...) %>% filter(dropmissing == TRUE) %>% pull(var)
}
