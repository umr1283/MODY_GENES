## Do multigene analysis in HNP ##

#### head ----
options(stringsAsFactors = FALSE)

script_name <- "06_Multigene_analysis_HNP"
project_directory <- "./"
working_directory <- gsub("/PROJECT/", "/DATATMP/", project_directory)
output_directory <- paste(working_directory, script_name, sep = "/")
dir.create(path = output_directory, recursive = TRUE, showWarnings = FALSE, mode = "0777")

### Load packages
library(parallel)
library(scales) ## scale
library(readr) ## read_tsv and cols
library(writexl)
library(dplyr) ## %>% 
library(purrr) # map 
library(tidyr) ## pivot
library(tidyverse)

## input 
input_directory <- paste(working_directory, "03_Replic_HealthyNevada", sep = "/")

#### Functions
SourceDir <- paste0(project_directory, "/Scripts/utils")
source(paste0(SourceDir, "/Zoom_Genotype.R"))
source(paste0(SourceDir, "/mist.R"))
source(paste0(SourceDir, "/Summary_clinique.R"))
source(paste0(SourceDir, "/Effectif_mutation_CC.R"))
source(paste0(SourceDir, "/DO_QC.R"))
Effectif_adapted <- function(data, y_name) {
  if (ncol(data) == 1) {
    Effectifs_CC <- data.frame(matrix(
      data = NA, 
      nrow = 1, 
      ncol = 9, 
      dimnames = list(
        NULL, 
        c(
          "Variant", 
          map_chr(.x = transpose(expand.grid(
            stringsAsFactors = FALSE,
            c("genotype_Missing", "genotype_0", "genotype_1", "genotype_2"),
            c("CTRL", "CASE")
          )), .f = ~paste(., collapse = "_"))
        )
      )
    ))[-1, ]
    return(Effectifs_CC)
  }
  
  counts_table <- map(.x = c(0, 1), .f = function(.x) {
    line_sample <- which(data[, y_name] %in% .x)
    data <- data[line_sample, ]
    Mat_eff <- Zoom_Genotype(select(data,-y_name)) %>%
      select(Variant, Missing_genotype, G_0, G_1, G_2)
    
    Mat_percent <- Mat_eff %>% 
      transpose() %>% 
      map_df(.f = function(.l) {
        .l <- unlist(.l)
        out <- unlist(.l[-grep("Variant", names(.l))])
        as_tibble(matrix(
          data = paste0(round(
            as.numeric(out) / sum(as.numeric(out), na.rm = TRUE) * 100,
          digits = 4), "%"), 
          byrow = TRUE, 
          nrow = 1, dimnames = list(NULL, names(out))
        )) %>% 
          mutate(Variant = .l["Variant"])
      }) %>% 
      rename(
       Missing_genotype_percent = Missing_genotype,
       genotype_0_percent = G_0,
       genotype_1_percent = G_1,
       genotype_2_percent = G_2
      )
    
    percent_Eff <- inner_join(x = Mat_percent, y = Mat_eff, by = "Variant") %>% 
      mutate(
        genotype_Missing = paste0(Missing_genotype_percent, " (", Missing_genotype, ")"),
        genotype_0 = paste0(genotype_0_percent, " (", G_0, ")"),
        genotype_1 = paste0(genotype_1_percent, " (", G_1, ")"),
        genotype_2 = paste0(genotype_2_percent, " (", G_2, ")")
      ) %>% 
      select(Variant, genotype_Missing, genotype_0, genotype_1, genotype_2)
    
    new_colnames <- switch(
      EXPR = as.character(.x),
      "1" = colnames(percent_Eff)[-1] <- paste0(colnames(percent_Eff)[-1], "_CASE"),
      "0" = colnames(percent_Eff)[-1] <- paste0(colnames(percent_Eff)[-1], "_CTRL")
    )
    colnames(percent_Eff)[-1] <- new_colnames
    return(percent_Eff)
  })
  
  Effectifs_CC <- inner_join(x = counts_table[[1]], y = counts_table[[2]], by = "Variant")
  return(Effectifs_CC)
}
Summary_clinique_adapted <- function(data, x, na_symbol = "NA") {
  continuous_trait <- function(x, digits = 3) {
    x_mean <- mean(x, na.rm = TRUE)
    x_sd <- sd(x, na.rm = TRUE)
    if (is.na(x_mean)) {
      txt <- na_symbol
    } else {
      txt <- paste0(
        signif(x_mean, digits),
        " (", signif(x_sd, digits), ")"
      )
    }
    return(txt)
  }

  mean_sd_n <- function(x) {
    if (all(is.na(x))) {
      na_symbol ## "\\-"
    } else {
      paste(
        paste(
          format(mean(x, na.rm = TRUE), digits = 2, nsmall = 2, drop0trailing = FALSE),
          format(sd(x, na.rm = TRUE), digits = 2, nsmall = 2, drop0trailing = FALSE),
          sep = "±" # "\\$\\pm\\$"
        ),
        paste0(
          "(n=",
          sum(!is.na(x)),
          ")"
        )
      )
    }
  }
  ## drop0trailing = FALSE keep 0 in right

  out <- data %>%
    mutate(Status = get(x)) %>%
    filter(!is.na(Status)) %>%
    group_by(Status) %>%
    summarise(
      N = format(n(), big.mark = ",", scientific = FALSE),
      AGE = mean_sd_n(age_now),
      BMI = mean_sd_n(median_BMI),
      SEX = paste0(
        "M:", format(sum(!is_female), big.mark = ",", scientific = FALSE),
        " / F:", format(sum(is_female), big.mark = ",", scientific = FALSE)
      )
    ) %>%
    mutate(Status = c("0" = "Controls", "1" = "Cases")[as.character(Status)]) %>%
    as.data.frame() %>%
    column_to_rownames(var = "Status") %>%
    t() %>%
    as.data.frame() %>%
    rownames_to_column(var = x)

  if (ncol(out) == 1) {
    # out <- out %>% mutate(Data = na_symbol)
    out <- NULL
  }

  return(out)
}


#### READ DATA ever QCed ----
alldta <- read_rds(path = paste0(input_directory, "/DATA_QCed.rds"))
Genes_Variants <- read_rds(path = paste0(input_directory, "/Genes_Variants.rds"))
genes_of_interest <- names(Genes_Variants)
length(genes_of_interest)

#### def group ----
## 6/ ## 8 gene of interest
## only considered "8 gènes analysés" 
Genes_Variants03 <- readr::read_rds(
  path = paste0(
    gsub(script_name, "03_Replic_HealthyNevada", output_directory), "/Genes_Variants.rds")
)
##  8 genes seen in replication HNP 
Group6 <- sort(names(Genes_Variants03))
Group6
length(Group6)
Group6_final <- Group6

## 7/ 7 gene of interest (sans than group6 without HNF1B)
Genes_Variants02 <- readr::read_rds(
  path = paste0(
    gsub(script_name, "02_Replic_UKBioBank", output_directory), "/Genes_Variants.rds")
)
##  7 genes seen in replication UK 
Group7 <- sort(names(Genes_Variants02))
Group7
length(Group7)
Group7_final <- Group7

y_name <- "T2D_status"
binary <- TRUE
covar_names <- covariates_names <- c("age_now", "is_female", "median_BMI", "PC1", "PC2", "PC3", "PC4", "PC5") 


for (igroup in c("Group6", "Group7")) {
  message(">>> >>> ", igroup)
  
  Genes_Variants <- Genes_Variants03 ## genes and variants dispo in HNP
  if(igroup=="Group7"){
    ## keep only genes in common
    Genes_Variants <- Genes_Variants[intersect(x = names(Genes_Variants), names(Genes_Variants02))]
  }
  message(paste(names(Genes_Variants), collapse = ", "))
  
  variants <- unlist(Genes_Variants, use.names = FALSE)
  G_rare <- alldta[, variants, drop = FALSE]

  #### multigene mist ----

  message(paste("[MULTI-GENE]", "After all QCs,", ncol(G_rare), "rare variants will be analysed."))
  
  info_filtred <- Genes_Variants
  for (i in 1:length(Genes_Variants)) {
    tmp <- as.data.frame(Genes_Variants[[names(Genes_Variants)[i]]]) 
    names(tmp) <- "VARIANTS"
    tmp$GENE <- names(Genes_Variants)[i]
    info_filtred[[names(Genes_Variants)[i]]] <- tmp 
  }
  info_filtred <- bind_rows(info_filtred)
  
  Z <- model.matrix(~GENE -1, data = info_filtred)
  rownames(Z) <- info_filtred$VARIANTS
  G_rare <- G_rare[, info_filtred$VARIANTS] 
  
  out <- try(
    mist(
      y = alldta[, y_name],
      X = alldta[, covar_names, drop = FALSE],
      G = G_rare,
      Z = Z,
      method = "liu",
      model = "binary"
    ),
    silent = TRUE
  )
  
  if (!is(out, "try-error")) {
    stat_rare <- as.data.frame(out$out_rare) 
    rownames(stat_rare) <- gsub("^GZGENE", "", rownames(stat_rare))
    stat_rare <- cbind.data.frame("SubClusters" = rownames(stat_rare), stat_rare)
    stat_MiST <- as.data.frame(out$out_MiST)
    stat_MiST$ErrorReason <- attr(out, "condition")$message
  } else {
    message(paste("[MULTI-GENE]", paste0(attr(out, "condition")$message)))
    stat_rare <- data.frame(
      matrix(
            NA, 
            nrow = 1,
            ncol = 7,
            dimnmes = list(NULL, c("SubClusters", "Pi_hat", "CI_2.5", "CI_97.5", "SE", "Pval", "OR"))
          )
    )
    stat_MiST <- data.frame(
      matrix(
        data = NA, 
        nrow = 1, 
        ncol = 6,
        dimnames = list(NULL, c("S.pi", "p.value.S.pi", "S.tau", "p.value.S.tau", "p.value.overall", "ErrorReason"))
      )
    )
    stat_MiST$ErrorReason <- attr(out, "condition")$message
  }
  stat_MiST$nb_rare_var <- ncol(G_rare)
  stat_rare$SubClusters <- gsub("^GZGene", "", stat_rare$SubClusters)
  
  stat_MiST$trait = y_name
  stat_MiST$covariates = paste(covar_names, collapse = ";")
  stat_MiST$sample_size = nrow(alldta)
  
  saveRDS(object = stat_MiST, paste0(output_directory, "/stat_MiST_",igroup,".RDS"))
  saveRDS(object = stat_rare, paste0(output_directory, "/stat_rare_",igroup,".RDS"))
  message("MiST analysis done and writen !")
  
  
  #### then ---
  ## clinical information
  clinicalData <- Summary_clinique_adapted(
    data = alldta[,c("T2D_status", covar_names)], 
    x = "T2D_status", 
    na_symbol = "NA"
  )
  
  ## Annot
  myEffectif <- Effectif_adapted(
    data = alldta[, c("T2D_status", variants)],
    y_name = "T2D_status"
  )
  myZoom <- Zoom_Genotype(MatriceG = alldta[, variants]) %>% 
    select(-c("ProperAlleleREF", "VarClean", "chr_pos"))
  MyAnnotation <- full_join(x = myZoom, y = myEffectif, by = "Variant") %>% 
    left_join(x = ., y = rename(info_filtred, "Variant" = "VARIANTS"), by = "Variant")  %>% 
    select(GENE, everything())
  
  
  writexl::write_xlsx(
    x = list("MyAnnotation" = MyAnnotation, "clinicalData" = clinicalData), 
    path = paste0(output_directory, "/multigene_HNP_",igroup,".xlsx"), 
    col_names = TRUE
  )
  message("xlsx writen ! ")

}



  
