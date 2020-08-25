#!/usr/bin/env Rscript

PARAM_remise <- TRUE
message("is remise ? ", PARAM_remise)
# PARAM_remise <- FALSE 
all_samples_size <- c(seq(from = 500, to = 6000, by = 500), 9999)
simu_sample_size <- all_samples_size

#############################################################
#### MAKE ANALYSES calcul de la puissance - samples size #### 
#############################################################

# "ACMG-CUSTOM"  
## CUSTOM analyses of ** 8 actionable ** GENES DU DIAB MONOGENIQUE

options(stringsAsFactors = FALSE)
args <- commandArgs(trailingOnly = TRUE)
options(warn = 2) ## Turn warning into error => can be trapped

project_directory <- "./"

#### Load packages ---------------------------------------------------------------------------------
suppressPackageStartupMessages(
  invisible(
    lapply(
      c(
        "datasets", "utils", "grDevices", "graphics", "stats", "methods",
        "parallel", "scales", "readxl", "writexl", "tidyverse", "mongolite"
      ),
      library,
      character.only = TRUE
    )
  )
)

RawDir <- "/sequencing_data/"
PhenotypeDir <- "/phenotypes_data/"
Phenotype_file <- paste0(PhenotypeDir, "Phenotyes_latest.xlsx")
SourceDir <- paste0(project_directory, "/Scripts/utils")
InputDir <- paste0(project_directory, "/Data/01_MAKE_DATASET/")


#### Functions -------------------------------------------------------------------------------------
source(paste0(SourceDir, "/MiST.R"))
source(paste0(SourceDir, "/Parameters-class.R"))
source(paste0(SourceDir, "/QC-class.R"))
source(paste0(SourceDir, "/Summary_clinique.R"))
source(paste0(SourceDir, "/Zoom_genotype.R"))
source(paste0(SourceDir, "/Effectif_mutation_YCC.R"))
source(paste0(SourceDir, "/Effectif_mutation_Yquanti.R"))
source(paste0(SourceDir, "/Define_QC_Object.R"))
source(paste0(SourceDir, "/DO_QC.R"))


f <- function(objectParams) {
  slot_names <- slotNames(objectParams)
  out <- as.list(rep(NA, length(slot_names))) %>%
    `names<-`(slot_names)
  for (islot in slot_names) {
    out[[islot]] <- objectParams[islot]
  }
  return(as.data.frame(out))
}
f2 <- function(objectQC) {
  out <- data.frame(matrix(objectQC["Excluded"], byrow = TRUE, ncol = 1)) %>%
    `names<-`(paste0("QC_", objectQC["Step"], "_", objectQC["Type"]))
  return(out)
}
f3 <- function(listf2) {
  nrowQCtable <- max(unlist(lapply(listf2, nrow)))
  out <- data.frame(matrix(nrow = nrowQCtable, ncol = length(listf2)))
  names(out) <- unlist(lapply(listf2, names))
  for (el in listf2) {
    if (nrow(el) < nrowQCtable) {
      out[, names(el)] <- c(el[, names(el)], rep(NA, nrowQCtable - nrow(el)))
    } else {
      out[, names(el)] <- el
    }
  }
  return(out)
}
f4 <- function(objectQC) {
  listout <- list(
    "Why" = objectQC["Why"],
    "Type" = objectQC["Type"],
    "N_origin" = objectQC["N_origin"],
    "Excluded" = objectQC["Excluded"],
    "N_excluded" = objectQC["N_excluded"],
    "Comment" = objectQC["Comment"]
  )
  return(listout)
}

#### DEF -------------------------------------------------------------------------------------------
LARGEdefSNP <- c(
  "initiator_codon_variant", "synonymous_variant", "stop_retained_variant", "stop_lost", "stop_gained",
  "missense_variant", "incomplete_terminal_codon_variant", "coding_sequence_variant", "splice_acceptor_variant",
  "splice_donor_variant", "splice_region_variant", "regulatory_region_variant"
)
STRICTdefSNP <- c(
  "initiator_codon_variant", "stop_retained_variant", "stop_lost", "stop_gained", "missense_variant",
  "splice_donor_variant", "splice_acceptor_variant"
)
LARGEdefINDEL <- c("all")
STRICTdefINDEL <- c("frameshift_variant", "inframe_insertion", "splice_acceptor_variant", "splice_donor_variant")

#### SET PARAMETERS --------------------------------------------------------------------------------

analyse_LARGE_STRICT <- "STRICT"   
## TO perform  the "ACMG-CUSTMOM" analyses, start form STRICT QC and removed other variant in QC11 based on CUSTOM FILE ! 
y_name <- "CC_T2D61" 
binary <- TRUE
SNP_INDEL <- "BOTH"
UseCustomFile <- TRUE
PopFilter <- "ALL"
stopifnot(PopFilter%in%c("EUR", "ALL"))

for(samples_size_pourcent in simu_sample_size) { 
  ## 9999 will means the whole sample size
  
  if(PARAM_remise) {
    OutDir <- paste0(".//Data/08_Simu_T2Dmono/simu_avecRemise/") 
  } else {
    OutDir <- paste0(".//Data/08_Simu_T2Dmono/simu_sansRemise/")
  }
  
  CustomFile <- ".//Docs/CUSTOM_files/custom_03102019.xlsx"
  
  
  #### Selection of genes   -------------------------------------------------------------
  
  
  GroupeGene <- "CUSTOM_T2Dmono_CCDiabete" ## ==> But filter only 8 actionable genes
  
  liste_TR <- read.table(
    file = paste0(RawDir, "Docs/Latest_info_gene.csv"),
    header = TRUE,
    sep = ";"
  ) %>%
    filter(get(GroupeGene) == 1) %>%
    select(gene_name, Transcript) %>%
    filter(!gene_name %in% c("BOLA2", "BOLA2B", "DTX2P1-UPK3BP1-PMS2P11", "RMST"))
  
  GroupeGene <- "8_actionable_genes_T2Dmono"
  genes_actionable <- unlist(strsplit(x = "HNF1A, GCK, GATA4, KCNJ11, HNF1B, ABCC8, HNF4A, GATA6", split = ", "))
  liste_TR <- liste_TR %>% 
    filter(gene_name %in% genes_actionable)
  
  
  message(paste("[MONO-GENE]", "<===>"))
  message(paste("[DATE-TIME]",Sys.time()))
  message(paste("[MONO-GENE]", "Group selected:", GroupeGene))
  message(paste("[MONO-GENE]", "Analysis:", analyse_LARGE_STRICT, " ~> ACMG-CUSTMOM"))
  message(paste("[MONO-GENE]", "Pheno:", y_name))
  message(paste("[MONO-GENE]", "Binary:", binary))
  message(paste("[MONO-GENE]", SNP_INDEL, "(about SNP and INDEL) will be analysed."))
  message(paste("[MONO-GENE]", PopFilter, "(subset to do on Pop_pred?)."))
  
  message(paste("[MONO-GENE]", samples_size_pourcent , " over the sample size is concered for the study power."))
  
  
  #### BEGIN OF ANALYSIS -----------------------------------------------------------------------------
  
  for (T2Ddef in c(5.6, 5.8, 6, 6.1, 7)) {
    
    if(T2Ddef == 7 ) {CustomFile <- ".//Docs/CUSTOM_files/custom_21072020.xlsx"}
    
    #### LOOP ON ALL GENE-TRANSCRIPT SELECTED in GroupeGene ------------------------------------------
    for (iline in 1:nrow(liste_TR)) { 
      
      #### ALL parameters ----------------------------------------------------------------------------
      Params <- new.MODY_Parameters()
      Params["FolderOut"] <- OutDir
      Params["GroupeGene"] <- GroupeGene
      Params["Gene_Symbol"] <- liste_TR$gene_name[iline]
      Params["ENST"] <- liste_TR$Transcript[iline]
      Params["analyse_LARGE_STRICT"] <- analyse_LARGE_STRICT
      Params["threshold"] <- 0.01
      Params["ExludeSynonymous"] <- FALSE
      Params["cluster"] <- "none"
      Params["y_name"] <- y_name
      Params["binary"] <- binary
      Params["SNP_INDEL"] <- SNP_INDEL
      # Params
      
      message(paste(
        "[MONO-GENE]", ">>> Go on", Params["Gene_Symbol"], Params["ENST"],
        "Version", Params["analyse_LARGE_STRICT"]
      ))
      
      ## if exist, the Dataset ready to analyses is here:
      dataset_file <- paste0(
        InputDir,
        Params["Gene_Symbol"], "/", Params["Gene_Symbol"], "_", Params["ENST"], ".RData"
      )
      
      if (!file.exists(dataset_file)) {
        message(paste("[MONO-GENE]", "Check in Log_gene DB what happened.")) 
        readr::write_file(
          x = paste(
            ">>>>\n\n", 
            "check `", dataset_file, "`\n", 
            "! not existing file !\n", 
            "\n\n"
          ), 
          path = paste0(Params["FolderOut"], "README.TXT"), append = TRUE)
      } else {
        ## Analysis can be done ##
        
        message(paste("[MONO-GENE]", "<Gene & Transcript ever seen, can load Phenotypes and RData..."))
        
        #### Phenotypes
        Phenotypes <- read_excel(
          path = Phenotype_file,
          col_names = TRUE,
          na = "NA",
          guess_max = 9365
        ) %>%
          as.data.frame()
        rownames(Phenotypes) <- paste(Phenotypes$RUN, Phenotypes$ID, sep = "_")

        ####//!\\ RE definie the T2D status with the threshold T2Ddef ##
        Phenotypes[["CC_T2D_simu"]] <- ifelse(
          test = Phenotypes[["CC_T2D61"]] %in% 1,
          yes = 1,
          no = ifelse(
            test = (
              Phenotypes[["selectInd"]]%in%"adult" & 
                matrixStats::rowMedians(as.matrix(Phenotypes[, grep("^FG_D", colnames(Phenotypes))]), na.rm = TRUE)<T2Ddef & 
                rowSums(Phenotypes[, grep("^AGE_D", colnames(Phenotypes))]>=40, na.rm = TRUE)>=1
            ) & (
              !is.na(
                Phenotypes[["selectInd"]]%in%"adult" & 
                  matrixStats::rowMedians(as.matrix(Phenotypes[, grep("^FG_D", colnames(Phenotypes))]), na.rm = TRUE)<T2Ddef &
                  rowSums(Phenotypes[, grep("^AGE_D", colnames(Phenotypes))]>=40, na.rm = TRUE)>=1
              )
            ), 
            yes = 0, 
            no = NA
          )
        )
        
        with(Phenotypes, table(CC_T2D56, CC_T2D_simu, useNA = "ifany"))
        y_name <- "CC_T2D_simu"
        Params["y_name"] <-  "CC_T2D_simu"
        
        #### RData ready to analyses
        load(dataset_file)
        
        ## Select covariates ##
        if (Params["y_name"] %in% c("CC_T2D61", "CC_T2D56", "CC_T2D_simu")) {
          covar_names <- c("AGE", "SEX", "BMI", "PC1", "PC2", "PC3", "PC4", "PC5") 
        }
        if (Params["y_name"] %in% c("CC_OBESITE")) {
          covar_names <- c("AGE", "SEX", "INFO_DIAB", "PC1", "PC2", "PC3", "PC4", "PC5")
        }
        if (Params["y_name"] %in% c("CC_Obesity_child")) {
          covar_names <- c("SEX", "AGE", "PC1", "PC2", "PC3", "PC4", "PC5")
        }
        if (Params["y_name"] %in% c("FG")) {
          covar_names <- c("SEX", "AGE", "BMI", "PC1", "PC2", "PC3", "PC4", "PC5")
        }
        
        #### TEST IF ANALYSIS EVER DONE --------------------------------------------------------------
        {
          
          #### RUN ANALYSIS --------------------------------------------------------------------------
          message(paste("[MONO-GENE]", "Go to run analysis...")) 
          
          #### KEEP ONLY SAMPLES INVOLVED IN THE STUDY ####
          message(paste("[MONO-GENE]", nrow(Phenotypes), "samples at the begining."))
          QC_PHENO <- Phenotypes %>%
            select("IID", Params["y_name"]) %>%
            `colnames<-`(c("IID", "y_name")) %>%
            filter(!is.na(y_name))
          
          ## Remove sample
          Phenotypes <- Phenotypes %>%
            filter(IID %in% QC_PHENO$IID)
          
          all_genotype_data <- all_genotype_data[Phenotypes$IID, , drop = FALSE]
          all_genotype_data$IID <- rownames(all_genotype_data)
          all_genotype_data <- all_genotype_data %>%
            filter(IID %in% QC_PHENO$IID)
          
          ## if option PopFilter == "EUR" ## 
          if (PopFilter %in% "EUR") {
            message("[MONO-GENE]", "Only EUR samples will be kept...")
            Phenotypes <- Phenotypes %>%
              filter(Pop_pred %in% "EUR")
            all_genotype_data <- all_genotype_data %>%
              filter(IID %in% Phenotypes$IID)
          }
          
          
          message(paste("[MONO-GENE]", nrow(Phenotypes), "samples involved in this study."))
          
          #### If custom file used, remove all frequent variant before QC ####
          if (UseCustomFile) {
            
            QC00 <- new.MODY_QC(
              Step = 0,
              Why = "Custom file used, so all frequent variants are removed before QC.",
              Type = "Variant"
            )
            QC00["N_origin"] <- ncol(all_genotype_data) - 1
            QC00["Excluded"] <- Detail_Geno[Detail_Geno$MAF > Params["threshold"], "Variant_Old"]
            
            Detail_Geno <- filter(.data = Detail_Geno, MAF <= Params["threshold"])
            
            all_genotype_data <- all_genotype_data %>% 
              select(c("IID", Detail_Geno$Variant_Old))
            
            all_annotation_data <- filter(all_annotation_data, Variant %in% Detail_Geno$Variant_Old)
            
            all_quality_data <- filter(all_quality_data, chr_pos %in% unique(all_annotation_data$chr_pos))
            
            QC00["N_excluded"] <- length(QC00["Excluded"])
            QC00["Comment"] <- paste0(
              "After QC", QC00["Step"], ", ",
              QC00["N_excluded"], " ",
              tolower(QC00["Type"]), ifelse(test = QC00["N_excluded"]>1, yes = "s", no = ""),
              " have been excluded leading to ",
              QC00["N_origin"] - QC00["N_excluded"], " in the dataset."
            )
            
            message(paste(
              "[MONO-GENE]", QC00["N_excluded"], 
              "frequent variants removed because custom analyses will use only rare variants."
            )
            )
          }
          message(paste("[MONO-GENE]", ncol(all_genotype_data) - 1, "variants ready to be analysed."))
          
          #### QC steps ------------------------------------------------------------------------------
          ## define QC object
          QC <- QC10 
          
          
          #### PARAMETER SNP_INDEL (BOTH , ...) ####
          if (SNP_INDEL == "BOTH") {
            
            #### DO 10 QC ####
            res_QC <- DO_QC(
              Params = Params, 
              covar_names = covar_names,
              LARGEdefSNP = LARGEdefSNP, 
              LARGEdefINDEL = LARGEdefINDEL,
              STRICTdefSNP = STRICTdefSNP, 
              STRICTdefINDEL = STRICTdefINDEL,
              CustomFile = CustomFile,
              QC = QC, 
              Phenotypes = Phenotypes,
              all_genotype_data = all_genotype_data, 
              all_annotation_data = all_annotation_data, 
              Detail_Geno = Detail_Geno, 
              all_quality_data = all_quality_data,
              qc_type = "classic", 
              GroupeGene = GroupeGene
            )
            QC <- res_QC$QC
            Phenotypes <- res_QC$Phenotypes
            all_genotype_data <- res_QC$all_genotype_data
            all_annotation_data <- res_QC$all_annotation_data
            Detail_Geno <- res_QC$Detail_Geno
            all_quality_data <- res_QC$all_quality_data
            
            QCdone <- TRUE
            message(paste("[MONO-GENE]", "QC10 Done"))
          } else {
            if (SNP_INDEL %in% c("SNP", "INDEL")) {
              
              #### QC0 FILTER ON SNP_INDEL ####
              message(paste("[MONO-GENE]", "QC0 FILTER ONLY", SNP_INDEL, "kept"))
              
              QC0 <- new.MODY_QC(
                Step = 0.1,
                Why = paste("Only", SNP_INDEL, "are selected in this analysis"),
                Type = "Variant"
              )
              QC0["N_origin"] <- ncol(all_genotype_data) - 1
              QC0["Excluded"] <- all_annotation_data$Variant[all_annotation_data$Type_Variant != SNP_INDEL]
              
              all_annotation_data <- all_annotation_data %>%
                filter(Type_Variant %in% SNP_INDEL)
              ## maj geno, quality and Detail_Geno
              all_genotype_data <- all_genotype_data %>%
                select(c("IID", names(all_genotype_data) %in% all_annotation_data$Variant))
              
              Detail_Geno <- Detail_Geno %>%
                filter(Variant_Old %in% names(all_genotype_data))
              
              all_quality_data <- all_quality_data %>%
                filter(chr_pos %in% unique(all_annotation_data$chr_pos))
              
              QC0["N_excluded"] <- QC0["N_origin"] - (ncol(all_genotype_data) - 1)
              QC0["Comment"] <- paste0(
                "After QC", QC0["Step"], ", ",
                QC0["N_excluded"]," ",
                tolower(QC0["Type"]), ifelse(test = QC0["N_excluded"]>1, yes = "s", no = ""),
                " have been excluded leading to ",
                QC0["N_origin"] - QC0["N_excluded"], " in the dataset."
              )
              
              if (ncol(all_genotype_data) > 0) {
                #### DO 10 QC ####
                res_QC <- DO_QC(
                  Params = Params, 
                  covar_names = covar_names,
                  LARGEdefSNP = LARGEdefSNP, 
                  LARGEdefINDEL = LARGEdefINDEL,
                  STRICTdefSNP = STRICTdefSNP, 
                  STRICTdefINDEL = STRICTdefINDEL,
                  CustomFile = CustomFile,
                  QC = QC, 
                  Phenotypes = Phenotypes,
                  all_genotype_data = all_genotype_data, 
                  all_annotation_data = all_annotation_data, 
                  Detail_Geno = Detail_Geno, 
                  all_quality_data = all_quality_data,
                  qc_type = "classic"
                )
                QC <- res_QC$QC
                Phenotypes <- res_QC$Phenotypes
                all_genotype_data <- res_QC$all_genotype_data
                all_annotation_data <- res_QC$all_annotation_data
                Detail_Geno <- res_QC$Detail_Geno
                all_quality_data <- res_QC$all_quality_data
                
                QCdone <- TRUE
                message(paste("[MONO-GENE]", "QC10 Done"))
              } else {
                QCdone <- FALSE
                message(paste("[MONO-GENE]", "QC10 IMPOSSIBLE"))
                all_genotype_data <- data.frame(matrix(nrow = 0, ncol = 0))
                Phenotypes <- data.frame(matrix(nrow = 0, ncol = 0))
              }
            } else {
              ## wrong value in SNP_INDEL
              message(paste("[MONO-GENE]", "INVALID value for SNP_INDEL parameter. MUST BE 'BOTH', 'SNP' or 'INDEL'."))
              stop(paste("[MONO-GENE]", "INVALID value for SNP_INDEL parameter. MUST BE 'BOTH', 'SNP' or 'INDEL'."))
            }
          }
          
          if (QCdone) {
            qc11done <- FALSE
            
            #### QC11 ####
            if (file.exists(CustomFile) & Params["analyse_LARGE_STRICT"] != "LARGE") {
              ## Read
              InfoCUSTOM <- read_excel(path = CustomFile, sheet = 1, col_names = TRUE, col_types = "text")
              
              if (nrow(InfoCUSTOM) > 0) {
                InfoCUSTOM <- InfoCUSTOM %>%
                  dplyr::filter(Gene == Params["Gene_Symbol"]) %>% 
                  dplyr::rename("VARIANT" = "CUSTOM")
                
                ## EXCLUD 0
                QC11 <- new.MODY_QC(
                  Step = 11,
                  Why = "Manual curation",
                  Type = "Variant"
                )
                QC11["N_origin"] <- ncol(all_genotype_data) - 1
                QC11["Excluded"] <- names(all_genotype_data)[!(names(all_genotype_data) %in% InfoCUSTOM$VARIANT)]
                QC11["Excluded"] <- QC11["Excluded"][-grep("IID", QC11["Excluded"])]
                ## maj geno, quality and Detail_Geno
                all_genotype_data <- all_genotype_data %>%
                  select(names(all_genotype_data)[!(names(all_genotype_data) %in% QC11["Excluded"])])
                Detail_Geno <- Detail_Geno %>%
                  filter(Variant %in% names(all_genotype_data))
                all_quality_data <- all_quality_data %>%
                  filter(chr_pos %in% unique(Detail_Geno$chr_pos))
                
                QC11["N_excluded"] <- QC11["N_origin"] - (ncol(all_genotype_data) - 1)
                QC11["Comment"] <- paste0(
                  "After QC", QC11["Step"], ", ",
                  QC11["N_excluded"]," ",
                  tolower(QC11["Type"]), ifelse(test = QC11["N_excluded"]>1, yes = "s", no = ""),
                  " have been excluded leading to ",
                  QC11["N_origin"] - QC11["N_excluded"], " in the dataset."
                )
                
                qc11done <- TRUE
              
                stopifnot(all(Detail_Geno$Variant %in% InfoCUSTOM$VARIANT))
                
              } else {
                ## do not have this analyses to run ... 
                #  exclud everything
                QC11 <- new.MODY_QC(
                  Step = 11,
                  Why = "Manual curation",
                  Type = "Variant"
                )
                QC11["N_origin"] <- ncol(all_genotype_data) - 1
                QC11["Excluded"] <- names(all_genotype_data)
                QC11["Excluded"] <- QC11["Excluded"][-grep("IID", QC11["Excluded"])]
                ## maj geno, quality and Detail_Geno
                # all_annotation_data  ## no more used
                all_genotype_data <- all_genotype_data %>%
                  select(names(all_genotype_data)[!(names(all_genotype_data) %in% QC11["Excluded"])])
                Detail_Geno <- Detail_Geno %>%
                  filter(Variant %in% names(all_genotype_data))
                all_quality_data <- all_quality_data %>%
                  filter(chr_pos %in% unique(Detail_Geno$chr_pos))
                
                QC11["N_excluded"] <- QC11["N_origin"] - (ncol(all_genotype_data) - 1)
                QC11["Comment"] <- paste0(
                  "After QC", QC11["Step"], ", ",
                  QC11["N_excluded"]," ",
                  tolower(QC11["Type"]), ifelse(test = QC11["N_excluded"]>1, yes = "s", no = ""),
                  " have been excluded leading to ",
                  QC11["N_origin"] - QC11["N_excluded"], " in the dataset."
                )
                
                qc11done <- TRUE
              }
            } 
            
            if (qc11done) {
              message(paste("[MONO-GENE]", "QC11: qcmanuel done."))
            } else {
              message(paste("[MONO-GENE]", "NO QC11."))
            }
          }
          
          if (ncol(all_genotype_data) == 1) {
            all_genotype_data <- data.frame(matrix(nrow = 0, ncol = 0))
            Phenotypes <- data.frame(matrix(nrow = 0, ncol = 0))
          }
          
          #### MAKE VARY THE N_size by re-sampling #### 
          
          back_up_Phenotypes <- Phenotypes
          back_up_all_genotype_data <- all_genotype_data
          back_up_Detail_Geno <- Detail_Geno
          
          trash <- mclapply(X = 1:100, mc.cores = 25, mc.preschedule = FALSE, FUN = function(simu_i) {
            if(samples_size_pourcent>nrow(back_up_Phenotypes)) {samples_size_pourcent <- nrow(back_up_Phenotypes)}
            
            message("Re sampling : ", samples_size_pourcent)
            message("Simulation : ", simu_i)
            
            get_percent_CaseControl <- sum(back_up_Phenotypes[Params["y_name"]]==1, na.rm = TRUE) / 
              (sum(back_up_Phenotypes[Params["y_name"]]==0, na.rm = TRUE) + sum(back_up_Phenotypes[Params["y_name"]]==1, na.rm = TRUE))
            n_Cases <- round(samples_size_pourcent * get_percent_CaseControl)
            n_Controls <- round(samples_size_pourcent - n_Cases)
            message(n_Cases, " cases vs ", n_Controls, " controls")
            # Echantion
            cases <- back_up_Phenotypes %>% 
              filter(get(y_name) %in% 1)
            samples_cases <- sample(x = cases$IID, size = n_Cases, replace = PARAM_remise)
            cases <- left_join(x = tibble(IID = samples_cases), y = cases, by = "IID")
            controls <- back_up_Phenotypes %>% 
              filter(get(y_name) %in% 0) 
            samples_ctrls <- sample(x = controls$IID, size = n_Controls, replace = PARAM_remise)
            controls <- left_join(x = tibble(IID = samples_ctrls), y = controls, by = "IID")
            Phenotypes <- bind_rows(cases, controls)
            message(nrow(Phenotypes), " are detected in Phenotypes")
            # keep Nsize in detail analyses and update object ! 
            all_dta <- left_join(x = Phenotypes, y = back_up_all_genotype_data, by = "IID")
            all_genotype_data <- as.data.frame(all_dta[, names(back_up_all_genotype_data)])
            Phenotypes <- all_dta[, names(Phenotypes)]
            
            Detail_Geno <- back_up_Detail_Geno
            Detail_Geno$samples_size_pourcent = samples_size_pourcent
            message("done")
            
            #### clinicalData --------------------------------------------------------------------------
            if (ncol(Phenotypes) != 0 & nrow(Phenotypes) != 0) {
              if (Params["binary"]) {
                clinicalData <- Summary_clinique(data = Phenotypes, x = Params["y_name"], na_symbol = "NA")
              } else {
                clinicalData <- summary(Phenotypes[, Params["y_name"]]) %>%
                  t() %>%
                  as.data.frame() %>%
                  select(-Var1) %>%
                  `colnames<-`(c(Params["y_name"], "Values")) %>%
                  mutate(
                    Values = round(x = Values, digit = 4)
                  )
                clinicalData[, 1] <- as.character(clinicalData[, 1])
                clinicalData <- rbind.data.frame(clinicalData, c("N", nrow(Phenotypes)))
              }
              clinicalData <- list(clinicalData)
              message(paste("[MONO-GENE]", "clinicalData Done"))
            } else {
              clinicalData <- NULL
              message(paste("[MONO-GENE]", "/!\\ No phenotypes"))
            }
            
            #### Analysis using MiST -------------------------------------------------------------------
            message(paste("[MONO-GENE]", "Analysis using MiST"))
            variantRARE <- Detail_Geno$Variant[Detail_Geno$MAF < Params["threshold"]]
            G_rare <- all_genotype_data[, names(all_genotype_data) %in% variantRARE, drop = FALSE]
            message(paste("[MONO-GENE]", "G_rare Done"))
            if (ncol(G_rare) > 1) { ## rare analysis possible
              message(paste("[MONO-GENE]", "After all QCs,", ncol(G_rare), "rare variants will be analysed."))
              
              # param cluster = "none" ## or "TYPE"
              if (Params["cluster"] == "TYPE") {
                info_filtred <- Detail_Geno[Detail_Geno$Variant %in% names(G_rare), ]
                if (length(unique(info_filtred$Consequence)) > 1) {
                  Z <- model.matrix(~ Consequence - 1, data = info_filtred) ## Z mat indicatrice
                } else { ## 1 cluster because 1 Type of consequence
                  Z <- matrix(rep(1, ncol(G_rare)), ncol = 1)
                }
                message(paste("[MONO-GENE]", "WARNING cluster=TYPE selected, out_mist must be adapted..."))
                
              } else { ## cluster = "none" => 1 cluster
                Z <- matrix(rep(1, ncol(G_rare)), ncol = 1)
              }
              
              if (Params["binary"]) {
                
                has_issues <- try(tools::assertCondition(out_mist <- mist(
                  y = Phenotypes[[Params["y_name"]]],  ## must be a vector !
                  X = Phenotypes[, covar_names, drop = FALSE], 
                  G = G_rare, 
                  Z = Z, 
                  model = "binary"
                )), silent = TRUE)
                has_issues <- has_issues[
                  sapply(has_issues, function(el) {length(intersect(class(el), c("warning", "error")))!=0})
                  ]
                if (length(has_issues)==0) {
                  # out_mist <- mist(y, X, G, Z, model)
                  out_mist <- mist_print(out_mist)
                  out <- tibble::tibble(
                    trait = Params["y_name"], 
                    covariates = paste(covar_names, collapse = ";"),
                    sample_size = nrow(Phenotypes),
                    error = "None",
                    mist_estimate = list(out_mist$estimate),
                    mist_statistics = list(out_mist$statistic)
                  ) 
                } else {
                  out <- tibble::tibble(
                    trait = Params["y_name"], 
                    covariates = paste(covar_names, collapse = ";"),
                    sample_size = nrow(Phenotypes),
                    error = paste(
                      unique(sapply(
                        X = has_issues, 
                        FUN = function(el) {paste0("[", class(el)[2], "] ", el$message)}
                      )),
                      collapse = ";\n"
                    ),
                    mist_estimate = ifelse(
                      test = any(
                        unique(sapply(
                          X = has_issues, 
                          FUN = function(el) {class(el)[2]}
                        )) != "warning") ,
                      yes = list(NA),
                      no = list(mist_print(out_mist)$estimate) ## if only warning keep estimates
                    ),
                    mist_statistics = ifelse(
                      test = any(
                        unique(sapply(
                          X = has_issues, 
                          FUN = function(el) {class(el)[2]}
                        )) != "warning") ,
                      yes = list(NA),
                      no = list(mist_print(out_mist)$statistic) ## if only warning keep statistic
                    )
                  ) 
                }
                stat_MiSTrare <- tidyr::unnest(data = out, cols = c(mist_estimate, mist_statistics)) 
                
              } else {
                ## analysis quanti
                
                has_issues <- try(tools::assertCondition(out_mist <- mist(
                  y = Phenotypes[[Params["y_name"]]],  ## must be a vector ! 
                  X = Phenotypes[, covar_names, drop = FALSE], 
                  G = G_rare, 
                  Z = Z, 
                  model = "continuous"
                )), silent = TRUE)
                has_issues <- has_issues[
                  sapply(has_issues, function(el) {length(intersect(class(el), c("warning", "error")))!=0})
                  ]
                if (length(has_issues)==0) {
                  # out_mist <- mist(y, X, G, Z, model)
                  out_mist <- mist_print(out_mist)
                  out <- tibble::tibble(
                    trait = Params["y_name"], 
                    covariates = paste(covar_names, collapse = ";"),
                    sample_size = nrow(Phenotypes),
                    error = "None",
                    mist_estimate = list(out_mist$estimate),
                    mist_statistics = list(out_mist$statistic)
                  ) 
                } else {
                  out <- tibble::tibble(
                    trait = Params["y_name"], 
                    covariates = paste(covar_names, collapse = ";"),
                    sample_size = nrow(Phenotypes),
                    error = paste(
                      unique(sapply(
                        X = has_issues, 
                        FUN = function(el) {paste0("[", class(el)[2], "] ", el$message)}
                      )),
                      collapse = ";\n"
                    ),
                    mist_estimate = ifelse(
                      test = any(
                        unique(sapply(
                          X = has_issues, 
                          FUN = function(el) {class(el)[2]}
                        )) != "warning") ,
                      yes = list(NA),
                      no = list(mist_print(out_mist)$estimate)
                    ),
                    mist_statistics = ifelse(
                      test = any(
                        unique(sapply(
                          X = has_issues, 
                          FUN = function(el) {class(el)[2]}
                        )) != "warning") ,
                      yes = list(NA),
                      no = list(mist_print(out_mist)$statistic)
                    )
                  )
                }
                stat_MiSTrare <- tidyr::unnest(data = out, cols = c(mist_estimate, mist_statistics)) 
                stat_MiSTrare$SubClusters <- ifelse(stat_MiSTrare$SubClusters%in%"M", "None", stat_MiSTrare$SubClusters)
                ## fix notation about fit in linear mist
                
                ## Add in clinicalData, mean Trait ~ rare variant count --
                
                tmp <- cbind.data.frame(Y = Phenotypes[, Params["y_name"]], G_rare)
                tmp$nb_mut = apply(X = tmp[,-1], MARGIN = 1, FUN = sum)
                
                meantrait_countmut <- tmp %>% 
                  group_by(nb_mut) %>% 
                  summarise(
                    n = n(), 
                    min_trait = min(Y),
                    q25_trait = quantile(Y, probs = 0.25),
                    median_trait = median(Y), 
                    q75_trait = quantile(Y, probs = 0.75), 
                    max_trait = max(Y),
                    mean_trait = mean(Y), 
                    sd_trait = sd(Y)
                  )
                
                bx_plot <- ggplot(data = tmp, mapping = aes(x = as.factor(nb_mut), y = Y)) + 
                  ggbeeswarm::geom_quasirandom(aes(color = as.factor(nb_mut)), width = 0.25, shape = 21) +
                  geom_violin(fill = "transparent") +
                  geom_boxplot(outlier.shape = NA, fill = "transparent", width = 0.25) + 
                  theme(legend.position = "none", axis.ticks = element_blank()) + 
                  scale_x_discrete(
                    breaks = meantrait_countmut$nb_mut,
                    labels = paste0(meantrait_countmut$nb_mut, "\n(N=", meantrait_countmut$n, ")")
                  ) + 
                  labs(
                    title = paste0("Boxplot of ", Params["y_name"], " per number of carried risk alleles"),
                    subtitle = paste0("MiST p.value.overall = ", round(x = stat_MiSTrare$p.value.overall, digits = 4)),
                    ## add stat_MiSTrare$p.value.overall in the fig
                    x = paste0("Number of carried risk alleles\namong the cluster of ", ncol(G_rare)," rare variants"), 
                    y = Params["y_name"], 
                    caption = "N, the number of individuals carrying the allele(s)."
                  ) 
                
                clinicalData <- append(x = clinicalData, values = list(bx_plot, meantrait_countmut))
              }
            } else {
              message(paste("[MONO-GENE]", "After all QCs, not enough rare variants."))
              ## After all QC, not enough rare variants... Rare analyses impossible
              if (Params["binary"]) {
                stat_MiSTrare <- data.frame(matrix(rep(NA, 14), byrow = TRUE, nrow = 1))
                colnames(stat_MiSTrare) <- c("trait", "covariates", "sample_size", "error", "SubClusters", 
                                             "Pi_hat", "CI_2.5", "CI_97.5", "OR", "S.pi", "p.value.S.pi", 
                                             "S.tau", "p.value.S.tau", "p.value.overall")
                stat_MiSTrare$error <- "[Dev msg] MiST Not applicable : After all QCs, not enough rare variants."
              } else {
                stat_MiSTrare <- data.frame(matrix(rep(NA, 13), byrow = TRUE, nrow = 1))
                colnames(stat_MiSTrare) <- c("trait", "covariates", "sample_size", "error", "SubClusters", 
                                             "Pi_hat", "CI_2.5", "CI_97.5", "S.pi", "p.value.S.pi", 
                                             "S.tau", "p.value.S.tau", "p.value.overall")
                stat_MiSTrare$error <- "[Dev msg] MiST Not applicable : After all QCs, not enough rare variants."
              }
            }
            
            ## Output
            # stat_MiSTrare
            stat_MiSTrare$nb_rare_var <- ncol(G_rare)
            stat_MiSTrare
            
            
            #### Effectif mutation  --------------------------------------------------------------------
            if (nrow(Phenotypes)>0) {
              if (Params["binary"]) {
                N_cases = sum(Phenotypes[ Params["y_name"]]==1, na.rm = TRUE)
                N_controles = sum(Phenotypes[ Params["y_name"]]==0, na.rm = TRUE)
                
                Effectif_mutation <- Effectif_mutation_YCC(
                  Phenotypes = Phenotypes,
                  G_mat = all_genotype_data,
                  y_name = Params["y_name"]
                ) %>%
                  select(-starts_with("genotype_Missing"))
              } else {
                Effectif_mutation <- Effectif_mutation_Yquanti(
                  G_mat = all_genotype_data[, -grep("IID", names(all_genotype_data)), drop = FALSE]
                )
              }
            } else {
              if (Params["binary"]) {
                N_cases = 0
                N_controles = 0
                Effectif_mutation <- structure(
                  list(
                    Variant = character(0), genotype_0_CTRL = character(0),
                    genotype_1_CTRL = character(0), genotype_2_CTRL = character(0),
                    genotype_0_CASE = character(0), genotype_1_CASE = character(0),
                    genotype_2_CASE = character(0)
                  ),
                  row.names = integer(0), class = c(
                    "tbl_df",
                    "tbl", "data.frame"
                  )
                )
              } else {
                Effectif_mutation <- structure(
                  list(
                    Variant = character(0), genotype_0 = character(0),
                    genotype_1 = character(0), genotype_2 = character(0)
                  ),
                  row.names = integer(0), class = c(
                    "tbl_df",
                    "tbl", "data.frame"
                  )
                )
              }
            }
            
            #### Output reshape ------------------------------------------------------------------------
            
            ## reshape QCobj
            if (QCdone) {
              xxqc <- f(Params) %>%
                as_tibble() %>%
                select(-GroupeGene)
              
              if (UseCustomFile) { xxqc$QC00 <- list(f4(QC00)) }
              if (SNP_INDEL %in% c("SNP", "INDEL")) { xxqc$QC0 <- list(f4(QC0)) }
              xxqc$QC1 <- list(f4(QC[[1]]))
              xxqc$QC2 <- list(f4(QC[[2]]))
              xxqc$QC3 <- list(f4(QC[[3]]))
              xxqc$QC4 <- list(f4(QC[[4]]))
              xxqc$QC5 <- list(f4(QC[[5]]))
              xxqc$QC6 <- list(f4(QC[[6]]))
              xxqc$QC7 <- list(f4(QC[[7]]))
              xxqc$QC8 <- list(f4(QC[[8]]))
              xxqc$QC9 <- list(f4(QC[[9]]))
              xxqc$QC10 <- list(f4(QC[[10]]))
              if (qc11done) { xxqc$QC11 <- list(f4(QC11)) }
              
            } else {
              xxqc <- f(Params) %>%
                as_tibble() %>%
                select(-GroupeGene)
              if (UseCustomFile) { xxqc$QC00 <- list(f4(QC00)) }
              if (SNP_INDEL %in% c("SNP", "INDEL")) { xxqc$QC0 <- list(f4(QC0)) }
            }
            
            if (nrow(Detail_Geno) == 0) {
              Annotation_qced <- left_join(
                x = Detail_Geno,
                y = all_quality_data %>%
                  select(-Type_Variant) %>%
                  slice(0),
                by = c("chr_pos")
              ) %>% 
                left_join(
                  x = .,
                  y = Effectif_mutation,
                  by = "Variant"
                ) %>% 
                mutate(Samples_with_mutation = NA)
            } else {
              Annotation_qced <- left_join(
                x = Detail_Geno,
                y = all_quality_data,
                by = c("chr_pos", "Type_Variant")
              ) %>%
                left_join(
                  x = .,
                  y = Effectif_mutation,
                  by = "Variant"
                )
            }
            
            res_rare <- stat_MiSTrare %>% 
              select(trait, covariates, sample_size, everything()) %>% 
              as.data.frame() 
            
            if(Params["binary"]) {
              res_rare <- res_rare %>% 
                mutate(
                  N_cases = N_cases, 
                  N_controles = N_controles 
                )
            }
            
            ## Data analysed at the end
            if (nrow(Annotation_qced) > 0) {
              Annotation_qced <- Annotation_qced %>%
                mutate(
                  Samples_with_mutation = map_chr(
                    .x = Variant,
                    .f = function(x) {
                      all_genotype_data[, c(x, "IID"), drop = FALSE] %>%
                        filter(get(x) != 0) %>%
                        `[[`("IID") %>%
                        paste(., collapse = ";")
                    }
                  )
                ) %>%
                select(Variant, Variant_Old, everything())
              ## For frequent variant Samples_with_mutation is too long for excel (limit of 32,767 characters)
              Annotation_qced[Annotation_qced$MAF > 0.05, "Samples_with_mutation"] <- "..."
              
              Annotation_qced <- Annotation_qced %>% 
                mutate(list_ind_mut = NULL, INFO = NULL)
            }
            
            ####  SAVE it -------------------------------------------------------------------
            message(">>> SAVE <<<")
            
            xxr <- f(Params) %>%
              as_tibble() %>%
              select(-c("GroupeGene", "FolderOut"))
            xxr$res_rare <- list(res_rare)
            xxr$qc11_done <- list(qc11done)
            
            write_tsv(
              x = xxr %>% unnest(cols = c(res_rare, qc11_done)), 
              path = paste0(
                Params["FolderOut"],
                simu_i, "_",
                "T2D", T2Ddef, "_",
                samples_size_pourcent, "_",
                Params["Gene_Symbol"],
                "_xxr", ".tsv"
              ), 
              append = FALSE,
              col_names = TRUE
            )
            
            #### WRITE ALL RES to send by xlsx ---------------------------------------------------------
            
            listf2 <- lapply(X = QC, FUN = f2)
            if (qc11done | UseCustomFile) { listf2 <- lapply(X = c(QC00, QC, QC11), FUN = f2) }
            
            QC_table <- f3(listf2)
            QC_table <- QC_table[, -grep("_noExcluded", names(QC_table))]
            if (nrow(Annotation_qced) > 0) {
              write_xlsx(
                x = list(Annotation_qced, QC_table),
                path = paste0(
                  Params["FolderOut"],
                  simu_i, "_",
                  "T2D", T2Ddef, "_",
                  samples_size_pourcent, "_", 
                  Params["Gene_Symbol"], "_", 
                  Params["ENST"], "_", 
                  Params["y_name"], ".xlsx"
                ),
                col_names = TRUE
              )
              message(paste("[MONO-GENE]", "annotation_qced and QC_table xlsx"))
            }
            
            ## End of the analysis
          }) ## end of simu_i
        }
      }
      
      ## next iline
    }
    
    
    message(paste("[DATE-TIME]",Sys.time()))
  } ## end loop about T2D definition threshold
  
} ## end of loop about samples_size_pourcent

message("end simu")
