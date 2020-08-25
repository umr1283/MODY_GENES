#!/usr/bin/env Rscript

# "ACMG-CUSTOM"  
## CUSTOM analyses of 8 GENES DU DIAB MONOGENIQUE

options(stringsAsFactors = FALSE)
# args <- commandArgs(trailingOnly = TRUE)
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
InputDir <- paste0(project_directory, "/Data/00_MAKE_DATASET/")


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

GroupeGene <- "CUSTOM_T2Dmono_CCDiabete"
analyse_LARGE_STRICT <- "STRICT"   
## TO perform  the "ACMG-CUSTMOM" analyses, start form STRICT QC and removed other variant in QC11 based on CUSTOM FILE ! 
y_name <- "CC_T2D7" 
binary <- TRUE
SNP_INDEL <- "BOTH"
UseCustomFile <- TRUE
PopFilter <- "ALL"
stopifnot(PopFilter%in%c("EUR", "ALL"))

OutDir <- ".//Data/09_ACMG_CUSTOM_T2D7/" 
dir.create(
  path = OutDir,
  recursive = TRUE,
  showWarnings = FALSE,
  mode = "0777"
)
CustomFile <- ".//Docs/CUSTOM_files/custom_03102019.xlsx"

InfoCUSTOM <- read_excel(path = CustomFile, sheet = 1, col_names = TRUE, col_types = "text")

## update custom file for T2D7 analyses ## 20/07/2020
updateCustom <- ".//Docs/CUSTOM_files/Mut-prediabetique_20072020.xlsx"
updateCustomdf <- read_excel(path = updateCustom, sheet = 1, col_names = TRUE, col_types = "text") %>% 
  dplyr::filter(ACMG %in% "1") %>%  
  dplyr::select(Gene_Symbol, Variant) %>% 
  dplyr::rename("CUSTOM" = "Variant", "Gene" = "Gene_Symbol")
newcustomfile <- ".//Docs/CUSTOM_files/custom_21072020.xlsx"
writexl::write_xlsx(x = unique(bind_rows(InfoCUSTOM, updateCustomdf)), path = newcustomfile, col_names = TRUE)
file.copy(from = newcustomfile, to = paste0(OutDir, basename(newcustomfile)))
CustomFile <- newcustomfile


#### Selection of genes and connexion  -------------------------------------------------------------

liste_TR <- read.table(
  file = paste0(RawDir, "Docs/Latest_info_gene.csv"),
  header = TRUE,
  sep = ";"
) %>%
  filter(get(GroupeGene) == 1) %>%
  select(gene_name, Transcript) %>%
  filter(!gene_name %in% c("BOLA2", "BOLA2B", "DTX2P1-UPK3BP1-PMS2P11", "RMST"))



message(paste("[MONO-GENE]", "<===>"))
message(paste("[DATE-TIME]",Sys.time()))
message(paste("[MONO-GENE]", "Group selected:", GroupeGene))
message(paste("[MONO-GENE]", "Analysis:", analyse_LARGE_STRICT, " ~> ACMG-CUSTMOM"))
message(paste("[MONO-GENE]", "Pheno:", y_name))
message(paste("[MONO-GENE]", "Binary:", binary))
message(paste("[MONO-GENE]", SNP_INDEL, "(about SNP and INDEL) will be analysed."))
message(paste("[MONO-GENE]", PopFilter, "(subset to do on Pop_pred?)."))


#### BEGIN OF ANALYSIS -----------------------------------------------------------------------------

 

# for (analyse_LARGE_STRICT in c("STRICT", "ACMG", "CUSTOM")) 
{

  
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
    
    message(paste(
      "[MONO-GENE]", ">>> Go on", Params["Gene_Symbol"], Params["ENST"],
      "Version", Params["analyse_LARGE_STRICT"]
    ))
    
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
      
      dir.create(
        path = paste0(Params["FolderOut"], Params["Gene_Symbol"]),
        recursive = TRUE,
        showWarnings = FALSE,
        mode = "0777"
      )
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
  
      
      #### RData ready to analyses
      load(dataset_file)
      
      ## Select covariates 
      if (Params["y_name"] %in% c( "CC_T2D61", "CC_T2D56", "CC_T2D7")) {
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
          ##  Curation manuel (to apply for analyses STRICT and CUSTOM)
          
          if (file.exists(CustomFile) & Params["analyse_LARGE_STRICT"] != "LARGE") {
            ## Read
            InfoCUSTOM <- read_excel(path = CustomFile, sheet = 1, col_names = TRUE, col_types = "text")
            
            InfoCUSTOM <- InfoCUSTOM %>%
              filter(Gene == Params["Gene_Symbol"]) ## --here adapt on 03102019 custom
              
            if (nrow(InfoCUSTOM) > 0) {
              InfoCUSTOM <- InfoCUSTOM %>%
                dplyr::rename(VARIANT = "CUSTOM")
              
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
              
              stopifnot(all(Detail_Geno$Variant %in% InfoCUSTOM$VARIANT))
              # stopifnot(nrow(InfoCUSTOM)==nrow(Detail_Geno))
            
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
        
        #### clinicalData --------------------------------------------------------------------------
        if (ncol(Phenotypes) != 0 & nrow(Phenotypes) != 0) {
          if (Params["binary"]) {
            clinicalData <- Summary_clinique(data = Phenotypes, x = Params["y_name"], na_symbol = "NA")
          } else {
            clinicalData <- summary(Phenotypes[, Params["y_name"]]) %>%
              t() %>%
              as.data.frame() %>%
              select(-Var1) %>%
              # rename(Values = Freq) %>%
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
              y = Phenotypes[, Params["y_name"]], 
              X = Phenotypes[, covar_names, drop = FALSE], 
              G = G_rare, 
              Z = Z, 
              model = "binary"
            )), silent = TRUE)
            has_issues <- has_issues[
              sapply(has_issues, function(el) {length(intersect(class(el), c("warning", "error")))!=0})
            ]
            if (length(has_issues)==0) {
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
              y = Phenotypes[, Params["y_name"]], 
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
        
        #### Analysis using Logistic regression ----------------------------------------------------
        
        variantFREQ <- Detail_Geno$Variant[Detail_Geno$MAF >= Params["threshold"]]
        G_freq <- all_genotype_data[, names(all_genotype_data) %in% variantFREQ, drop = FALSE]
        
        if (Params["binary"]) {
          cnames <- c("Y_name", "Covar_names", "Variant", "Beta_Variant", "CI_2.5", "CI_97.5", "Std_Error", "P_val", "OR", "ErrorReason")
        } else {
          cnames <- c("Y_name", "Covar_names", "Variant", "Beta_Variant", "CI_2.5", "CI_97.5", "Std_Error", "P_val", "ErrorReason")
        }
        
        if (Params["binary"]) {
          if (ncol(G_freq) > 0) {
            colsummaryY <- c("Distribution.per.variants", paste(Params["y_name"], "=", levels(as.factor(Phenotypes[, Params["y_name"]]))))
            
            message(paste(
              "[MONO-GENE]",
              "After all QCs,",
              ifelse(
                test = ncol(G_freq) == 1,
                yes = "1 frequent variant is analysed.",
                no = paste0(ncol(G_freq), " frequent variants are analysed.")
              )
            ))
            
            stat_freq_detail <- lapply(seq_len(ncol(G_freq)), function(s) {
              # print(s)
              if ("SampletoKeep" %in% ls()) {
                rm(SampletoKeep)
              }
              SampletoKeep <- rownames(G_freq[(G_freq[, s, drop = FALSE] >= 0), , drop = FALSE]) ##  NA coded -1 so keep only lines with geno >= 0
              
              Y <- Phenotypes[SampletoKeep, Params["y_name"]]
              X <- Phenotypes[SampletoKeep, covar_names, drop = FALSE]
              G_filtred <- G_freq[SampletoKeep, s]
              
              # n
              headcount <- cbind.data.frame(
                names(G_freq)[s],
                matrix(table(as.factor(Y)), byrow = TRUE, nrow = 1)
              )
              colnames(headcount) <- colsummaryY
              
              if (var(G_filtred) == 0) {
                ## no test possible
                ligne <- cbind.data.frame(
                  y_name, paste(covar_names, collapse = " + "),
                  names(G_freq)[s], NA, NA, NA, NA, NA, NA, "var = 0"
                )
              } else {
                prep_reg <- paste(" Y ~ ", paste(paste("X$", covar_names, sep = ""), collapse = " + "), " + G_filtred")
                if ("res" %in% ls()) {
                  rm(res)
                }
                res <- try(glm(as.formula(prep_reg), family = binomial), silent = TRUE)
                
                if (!is(res, "try-error")) {
                  b <- res$coefficients["G_filtred"]
                  if (is.na(b)) {
                    ligne <- cbind.data.frame(Params["y_name"], paste(covar_names, collapse = " + "), names(G_freq)[s], NA, NA, NA, NA, NA, NA, "is.na(b)==TRUE")
                  } else {
                    err_std <- summary(res)$coefficients["G_filtred", "Std. Error"]
                    p <- summary(res)$coefficients["G_filtred", "Pr(>|z|)"]
                    or <- exp(coef(res)["G_filtred"])
                    CI <- try(confint(res)["G_filtred", ], silent = TRUE)
                    if (!is(CI, "try-error")) {
                      CI <- matrix(CI, byrow = TRUE, nrow = 1)
                      ligne <- cbind.data.frame(Params["y_name"], paste(covar_names, collapse = " + "), names(G_freq)[s], b, CI, err_std, p, or, NA)
                    } else {
                      ligne <- cbind.data.frame(Params["y_name"], paste(covar_names, collapse = " + "), names(G_freq)[s], b, NA, NA, err_std, p, or, attr(CI, "condition")$message)
                    }
                  }
                } else {
                  ligne <- cbind.data.frame(Params["y_name"], paste(covar_names, collapse = "+"), names(G_freq)[s], NA, NA, NA, NA, NA, NA, attr(res, "condition")$message)
                }
              }
              colnames(ligne) <- cnames
              return(list(ligne = ligne, headcount = headcount))
            })
            
            stat_freq <- do.call("rbind", lapply(stat_freq_detail, "[[", "ligne"))
            stat_freq_n <- do.call("rbind", lapply(stat_freq_detail, "[[", "headcount"))
          } else {
            message(paste("[MONO-GENE]", "After all QCs, no more frequent variant."))
            ## after all QC, no more frequent variant
            stat_freq <- data.frame(matrix(rep(NA, length(cnames)), byrow = TRUE, nrow = 1))[-1, ]
            colnames(stat_freq) <- cnames
            colsummaryY <- c("Distribution.per.variants", paste(Params["y_name"], "=", levels(as.factor(c(0, 1)))))
            stat_freq_n <- data.frame(matrix(rep(NA, length(colsummaryY)), byrow = TRUE, nrow = 1))[-1, ]
            colnames(stat_freq_n) <- colsummaryY
          }
        } else {
          ## analysis freq quanti
          if (ncol(G_freq) > 0) {
            colsummaryY <- c("Distribution.per.variants", "SampleSize")
            
            message(paste(
              "[MONO-GENE]",
              "After all QCs,",
              ifelse(
                test = ncol(G_freq) == 1,
                yes = "1 frequent variants is analysed.",
                no = paste0(ncol(G_freq), " frequent variants are analysed.")
              )
            ))
            
            stat_freq_detail <- lapply(seq_len(ncol(G_freq)), function(s) {
              if ("SampletoKeep" %in% ls()) {
                rm(SampletoKeep)
              }
              SampletoKeep <- rownames(G_freq[(G_freq[, s, drop = FALSE] >= 0), , drop = FALSE]) ##  NA coded -1 so keep only lines with geno >= 0
              
              Y <- Phenotypes[SampletoKeep, Params["y_name"]]
              X <- Phenotypes[SampletoKeep, covar_names, drop = FALSE]
              G_filtred <- G_freq[SampletoKeep, s]
              
              headcount <- cbind.data.frame(
                names(G_freq)[s],
                length(G_filtred)
              )
              colnames(headcount) <- colsummaryY
              
              if (var(G_filtred) == 0) {
                ## no test possible
                ligne <- cbind.data.frame(
                  y_name, paste(covar_names, collapse = " + "),
                  names(G_freq)[s], NA, NA, NA, NA, NA, "var = 0"
                )
              } else {
                prep_reg <- paste(
                  "Y ~", paste(paste0("X$", covar_names), collapse = " + "), "+ G_filtred"
                )
                if ("res" %in% ls()) {
                  rm(res)
                }
                res <- try(lm(as.formula(prep_reg)), silent = TRUE)
                
                if (!is(res, "try-error")) {
                  b <- res$coefficients["G_filtred"]
                  if (is.na(b)) {
                    ligne <- cbind.data.frame(Params["y_name"], paste(covar_names, collapse = " + "), names(G_freq)[s], NA, NA, NA, NA, NA, "is.na(b)==TRUE")
                  } else {
                    err_std <- summary(res)$coefficients["G_filtred", "Std. Error"]
                    p <- summary(res)$coefficients["G_filtred", "Pr(>|t|)"]
                    CI <- try(confint(res)["G_filtred", ], silent = TRUE)
                    if (!is(CI, "try-error")) {
                      CI <- matrix(CI, byrow = TRUE, nrow = 1)
                      ligne <- cbind.data.frame(Params["y_name"], paste(covar_names, collapse = " + "), names(G_freq)[s], b, CI, err_std, p, NA)
                    } else {
                      ligne <- cbind.data.frame(Params["y_name"], paste(covar_names, collapse = " + "), names(G_freq)[s], b, NA, NA, err_std, p, attr(CI, "condition")$message)
                    }
                  }
                } else {
                  ligne <- cbind.data.frame(Params["y_name"], paste(covar_names, collapse = "+"), names(G_freq)[s], NA, NA, NA, NA, NA, attr(res, "condition")$message)
                }
              }
              colnames(ligne) <- cnames
              return(list(ligne = ligne, headcount = headcount))
            })
            
            stat_freq <- do.call("rbind", lapply(stat_freq_detail, "[[", "ligne"))
            stat_freq_n <- do.call("rbind", lapply(stat_freq_detail, "[[", "headcount"))
          } else {
            message(paste("[MONO-GENE]", "After all QCs, no more frequent variant."))
            ## after all QC, no more frequent variant
            stat_freq <- data.frame(matrix(rep(NA, length(cnames)), byrow = TRUE, nrow = 1))[-1, ]
            colnames(stat_freq) <- cnames
            colsummaryY <- c("Distribution.per.variants", "SampleSize")
            stat_freq_n <- data.frame(matrix(rep(NA, length(colsummaryY)), byrow = TRUE, nrow = 1))[-1, ]
            colnames(stat_freq_n) <- colsummaryY
          }
        }
        
        #### Effectif mutation  --------------------------------------------------------------------
        if (nrow(Phenotypes)>0) {
          if (Params["binary"]) {
            Effectif_mutation <- Effectif_mutation_YCC(
              Phenotypes = Phenotypes,
              G_mat = all_genotype_data,
              y_name = Params["y_name"]
            ) %>%
              select(-starts_with("genotype_Missing"))
            ## subset column missing_genotype because there is no samples with missing values in the final analyses
          } else {
            Effectif_mutation <- Effectif_mutation_Yquanti(
              G_mat = all_genotype_data[, -grep("IID", names(all_genotype_data)), drop = FALSE]
            )
          }
        } else {
          if (Params["binary"]) {
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
          select(trait, covariates, sample_size, SubClusters, everything()) %>% 
          as.data.frame()
        
        res_freq <- inner_join(stat_freq, stat_freq_n, by = c("Variant" = "Distribution.per.variants"))
        if (nrow(stat_freq) > 0) {
          if (Params["binary"]) {
            res_freq$SampleSize <- paste0("Ctrls:", comma(res_freq[, 11]), " / Cases:", comma(res_freq[, 12]))
            res_freq <- res_freq[, -c(11, 12)] ## remove column CC_Diabete = 0 and CC_Diabete = 1
          }
        } else {
          res_freq <- inner_join(
            res_freq,
            data.frame(matrix(rep(NA, 2), byrow = TRUE, nrow = 1, dimnames = list(NULL, c("Variant", "SampleSize"))))[-1, ],
            by = "Variant"
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
        
        #### save it -------------------------------------------------------------------
		
        message(">>> SAVE RDS <<<")
        xxLOG <- f(Params) %>%
          as_tibble() %>%
          select(-c("GroupeGene")) %>%
          mutate(
            n_var_qced = nrow(Annotation_qced),
            n_freq_res = nrow(res_freq),
            n_MiST_res = nrow(res_rare)
          )
        if (!is.null(clinicalData)) {
          if(length(clinicalData)==1){
            xxLOG$clinicalData <- list(clinicalData) 
          } else {
            ## cas with quanti trait and ggplot to save
            ggsave(
              filename = paste0(
                Params["FolderOut"],
                Params["Gene_Symbol"], "/",
                Params["Gene_Symbol"], "_", 
                Params["ENST"], "_", 
                Params["y_name"], "_boxplot.png"
              ), 
              plot = clinicalData[2][[1]], 
              dpi = 100
            )
            xxLOG$clinicalData <- list(clinicalData[-2]) 
            
          }
        }
        
        saveRDS(
          object = xxLOG, 
          file = paste0(
            Params["FolderOut"], Params["Gene_Symbol"], "/",
            Params["Gene_Symbol"], 
            "_xxLog", ".RDS"
          )
        )

        saveRDS(
          object = xxqc, 
          file = paste0(
            Params["FolderOut"], Params["Gene_Symbol"], "/",
            Params["Gene_Symbol"], 
            "_xxqc", ".RDS"
          )
        )
        
    
        xxr <- f(Params) %>%
          as_tibble() %>%
          select(-c("GroupeGene", "FolderOut"))
        xxr$res_rare <- list(res_rare)
        xxr$qc11_done <- list(qc11done)
        
        saveRDS(
          object = xxr, 
          file = paste0(
            Params["FolderOut"], Params["Gene_Symbol"], "/",
            Params["Gene_Symbol"], 
            "_xxr", ".RDS"
          )
        )
        
        if (nrow(res_freq) > 0) {
          xxf <- res_freq %>%
            as_tibble() %>%
            cbind(
              .,
              f(Params) %>%
                as_tibble() %>%
                select(-c(GroupeGene, FolderOut))
            )
          
          saveRDS(
            object = xxf, 
            file = paste0(
              Params["FolderOut"], Params["Gene_Symbol"], "/",
              Params["Gene_Symbol"], 
              "_xxf", ".RDS"
            )
          )
          
          message(paste("[MONO-GENE]", "insert RES_FREQ"))
        } else {
          message(paste("[MONO-GENE]", "No insert RES_FREQ"))
        }
        
        
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
              Params["Gene_Symbol"], "/",
              Params["Gene_Symbol"], "_", 
              Params["ENST"], "_", 
              Params["y_name"], ".xlsx"
            ),
            col_names = TRUE
          )
          message(paste("[MONO-GENE]", "annotation_qced and QC_table xlsx"))
        }
        
        ## End of the analysis
      }
    }
    
    ## next iline
  }

  message(paste("[DATE-TIME]",Sys.time()))
}


message("END 38")
