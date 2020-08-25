## UK bio bank Replication on MODY ##
##  run 09/03/2020 ## 


#### head ----
options(stringsAsFactors = FALSE)
rm(list=ls())
script_name <- "02_Replic_UKBioBank"
project_directory <- "./"
SourceDir <- paste0(project_directory, "/Scripts/utils")
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
source(paste0(SourceDir, "/MiST.R"))
source(paste0(SourceDir, "/Summary_clinique.R"))

#### Read raw data ----
# input <- "/UKBioBank/MgD_UKB35987indiv_PLP_EXP085_ht320.tsv"
input <- "/UKBioBank/1583699503320_UKB_carriers_8genes_EXP-085-v3_ht221_20200306.tsv"
dta <- read_tsv(
  file = input, 
  col_names = TRUE, 
  col_types = cols(
    # .default = col_double(),
    'is_female' = col_logical(),
    'genes_mutated' = col_character(),
    'hgvsc_present' = col_character(),
    'hgvsp_present' = col_character(),
    # 'is_female_1' = col_logical(), 
    'is_female'= col_logical(),
    'UKB_T2D_case'= col_double(),
    'UKB_T2D_control'= col_double(),
    'number_of_P_LP_variants'= col_double(),
    'genes_mutated'= col_character(), 
    'hgvsc_present'= col_character(), 
    'hgvsp_present'= col_character(), 
    'n_P_LP'= col_double(),
    'SEX'= col_double(),
    'age_attended_enrollment'= col_double(),
    'diabetes_diagnosed_by_doctor_2443'= col_double(),
    'age_diabetes_diagnosed_2976'= col_double(),
    'age_diabetes_related_eye_disease_5901'= col_double(),
    'started_insulin_within_one_year_diabetes_diagnosis_2986'= col_double(),
    'glucose_30740'= col_double(), 
    'number_of_icd10'= col_double(),
    'number_of_E11'= col_double(),
    'number_of_E10'= col_double(),
    'number_of_R73'= col_double(),
    'number_of_O24'= col_double(),
    'T2D_self_report'= col_double(),
    'T1D_self_report'= col_double(),
    'diabetes_self_report'= col_double(),
    'gestational_diabetes_self_report'= col_double(),
    'diabetes_father'= col_double(),
    'diabetes_mother'= col_double(),
    'diabetes_siblings'= col_double(),
    'insulin_field6153'= col_double(),
    'insulin_field6177'= col_double(),
    'BMI_21001'= col_double(), 
    'alcohol_intake_1558'= col_double(),
    'smoking_status_20116'= col_double(),
    'smoking_current_1239'= col_double(),
    'smoking_past_1249'= col_double(),
    'glycated_hemoglobin_30750'= col_double(), 
    'PC1'= col_double(), 
    'PC2'= col_double(), 
    'PC3'= col_double(), 
    'PC4'= col_double(), 
    'PC5'= col_double(), 
    'PC6'= col_double(), 
    'PC7'= col_double(), 
    'PC8'= col_double(), 
    'PC9'= col_double(), 
    'PC10'= col_double(), 
    'age_glucose_test'= col_double(),
    'age_glycated_hemoglobin_test'= col_double()
  )
) %>% 
  tibble::rownames_to_column(., "ID_LINE")
dim(dta)
covariates_names <- c("age_attended_enrollment", "is_female", "BMI_21001", "PC1", "PC2", "PC3", "PC4", "PC5")
var_of_interest <- c("ID_LINE", "UKB_T2D_case", "UKB_T2D_control", covariates_names, "hgvsc_present", "number_of_P_LP_variants")

if (!file.exists(paste0(output_directory, "/DATA_QCed.rds"))) {

  dta$genes_list <- map(.x = dta$genes_mutated, .f = function(.x){
    myres <- unlist(strsplit(x = .x, split = ","))
    myres
    myres2 <- myres
    myres2 <- gsub("\"", "", myres2, fixed=TRUE)
    myres2 <- gsub("[", "", myres2, fixed=TRUE)
    myres2 <- gsub("]", "", myres2, fixed=TRUE)
    return(myres2)
  })
  genes_of_interest <- dta$genes_list %>% unlist %>% unique 
  genes_of_interest <- genes_of_interest[!genes_of_interest%in%""]

  #### Format ----
  var_of_interest <- c(var_of_interest, "genes_list")
  dta <- select(dta, var_of_interest) %>% 
    as.data.frame()

  Genes_Variants <- vector(mode = "list", length = length(genes_of_interest))
  names(Genes_Variants) <- genes_of_interest
  for (i_line in 1:nrow(dta)) {
    current_genes <- unlist(dta$genes_list[i_line])
    if (all(!current_genes%in%"")) {
      current_variants <- strsplit(x = dta$hgvsc_present[i_line], split = ",") %>% 
        unlist(.) %>% 
        gsub("\"", "", ., fixed=TRUE) %>% 
        gsub("[", "", ., fixed=TRUE) %>% 
        gsub("]", "", ., fixed=TRUE)  
      names(current_variants) <- current_genes
      for (i_gene in current_genes) {
        Genes_Variants[i_gene][[1]] <- unique(c(unlist(Genes_Variants[i_gene]), current_variants[i_gene]))
      }
    }
  }
  
  G <- map_df(.x = dta$hgvsc_present, .f = function(.x) {
    current_variants <- strsplit(x = .x, split = ",") %>% 
      unlist(.) %>% 
      gsub("\"", "", ., fixed=TRUE) %>% 
      gsub("[", "", ., fixed=TRUE) %>% 
      gsub("]", "", ., fixed=TRUE)  
      
    myG <- table(current_variants) %>% 
      as.data.frame %>% 
      mutate(current_variants = ifelse(current_variants=="", "NO_Variant", as.character(current_variants))) %>% 
      tidyr::pivot_wider(data = ., names_from = current_variants, values_from = Freq)
    return(myG)
  })
  
  G[is.na(G)] <- 0 ## Here in G all NA means no mutation on those gene are observed <=> genotype 0
  G$NO_Variant <- NULL ## no more need this column, remove it 
  
  ALL_DATA <- bind_cols(dta, G) ## Pheno and Geno ... 
   
  ALL_DATA <- ALL_DATA %>% 
    mutate(
      T2D_status = map2_dbl(.x = UKB_T2D_case, .y = UKB_T2D_control, .f = function(.x, .y){
        status <- NA
        if (.x == 1) {status <- 1}
        if (.y == 1) {status <- 0}
        return(status)
      })
    )
  
  write_rds(x = ALL_DATA, path = paste0(output_directory, "/ALL_DATA.rds"))
  
  #### QC ----
  
  ## QC1 ## Filter on phenotype.
  ## remove samples if missing data in covariates
  DATA_QCed <- ALL_DATA %>% 
    filter_at(.vars = c("T2D_status", covariates_names), .vars_predicate = ~!is.na(.))
  
  QC1_missingValues <- dta$ID_LINE[!dta$ID_LINE %in% DATA_QCed$ID_LINE] 
  
  write_rds(x = QC1_missingValues, path = paste0(output_directory, "/QC1_missingValues.rds"))
  write_rds(x = Genes_Variants, path = paste0(output_directory, "/Genes_Variants.rds"))
  write_rds(x = DATA_QCed, path = paste0(output_directory, "/DATA_QCed.rds"))
  
} else {
  DATA_QCed <- read_rds(path = paste0(output_directory, "/DATA_QCed.rds"))
  Genes_Variants <- read_rds(path = paste0(output_directory, "/Genes_Variants.rds"))
  genes_of_interest <- names(Genes_Variants)
}

#### mist Analyses ----
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
    Mat_eff <- Summary_clinique(select(data,-y_name)) %>%
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
      AGE = mean_sd_n(age_attended_enrollment),
      BMI = mean_sd_n(BMI_21001),
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

all_res <- map_df(.x = genes_of_interest, .f = function(gene_i) {
  print(gene_i)
  y <- DATA_QCed[, "T2D_status"]
  X <- as.matrix(DATA_QCed[, covariates_names])
  G <- as.matrix(DATA_QCed[, unlist(Genes_Variants[gene_i]), drop = FALSE])
  Z <- matrix(data = rep(1, ncol(G)), ncol = 1)
  
  ## mist
  resmist <- mist(y = y, X = X, G = G, Z = Z, method = "liu", model = "binary")
  
  ## clinical information
  clinicalData <- Summary_clinique_adapted(
    data = DATA_QCed[,c("T2D_status", covariates_names)], 
    x = "T2D_status", 
    na_symbol = "NA"
  )
  ## Annot
  myEffectif <- Effectif_adapted(
    data = DATA_QCed[, c("T2D_status", unlist(Genes_Variants[gene_i]))],
    y_name = "T2D_status"
  )
  myZoom <- Summary_clinique(MatriceG = G) %>% 
    select(-c("ProperAlleleREF", "VarClean", "chr_pos"))
  MyAnnotation <- full_join(x = myZoom, y = myEffectif, by = "Variant") %>% 
    mutate(Gene = gene_i) %>% 
    select(Gene, everything())
  
  writexl::write_xlsx(
    x = list("MyAnnotation" = MyAnnotation, "clinicalData" = clinicalData), 
    path = paste0(output_directory, "/", gene_i, ".xlsx"), 
    col_names = TRUE
  )
  message("xlsx writen ! ")
  
  myline <- cbind.data.frame(
    "Gene_Name" = gene_i, 
    "N_samples" = length(y), 
    resmist$out_rare, 
    "nb_rare_var" = ncol(G),
    resmist$out_MiST
  )
  return(myline)
})


writexl::write_xlsx(
  x = all_res %>% 
    arrange(p.value.overall), 
  path = paste0(output_directory, "/MiST_results.xlsx"), 
  col_names = TRUE
)

message("MiST analysis done and writen !")
