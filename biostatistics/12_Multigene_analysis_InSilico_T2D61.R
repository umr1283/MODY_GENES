#### MULTI GENE ANALYSIS ####

options(stringsAsFactors = FALSE)
VERSION <- "ACMGCUSTOM"
y_name <- "CC_T2D61"
binary <- TRUE
SNP_INDEL <- "BOTH"
GroupeGene <- "CUSTOM_T2Dmono_CCDiabete"

script_name <- "12_Multigene_analysis_InSilico_T2D61"
project_directory <- "./"
working_directory <- gsub("/PROJECT/", "/DATATMP/", project_directory)
output_directory <- paste(working_directory, script_name, sep = "/")

dir.create(path = output_directory, recursive = TRUE, showWarnings = FALSE, mode = "0777")

## Select covariates
if (y_name %in% c("CC_T2D61", "CC_T2D56", "CC_T2D7")) {
  covar_names <- c("AGE", "SEX", "BMI", "PC1", "PC2", "PC3", "PC4", "PC5") 
}
if (y_name %in% c("CC_OBESITE")) {
  covar_names <- c("AGE", "SEX", "INFO_DIAB", "PC1", "PC2", "PC3", "PC4", "PC5")
}
if (y_name %in% c("CC_Obesity_child")) {
  covar_names <- c("SEX", "AGE", "PC1", "PC2", "PC3", "PC4", "PC5")
}
if (y_name %in% c("FG")) {
  covar_names <- c("SEX", "AGE", "PC1", "PC2", "PC3", "PC4", "PC5")
}

SourceDir <- paste0(project_directory, "/Scripts/utils")
InputDir <- paste0(project_directory, "/Data/00_MAKE_DATASET/")
CustomDir <- paste0(project_directory, "/Docs/CUSTOM_files/")
CustomFile <- ""

RawDir <- "/sequencing_data/"
PhenotypeDir <- "/phenotypes_data/"
Phenotype_file <- paste0(PhenotypeDir, "Phenotyes_latest.xlsx")

outputPath <- paste0(output_directory ,"/", VERSION)
if (file.exists(outputPath)) {
  unlink(x = outputPath, recursive = TRUE) 
  dir.create(outputPath)
} else {
  dir.create(outputPath)
}


### Load packages
suppressPackageStartupMessages(
  lapply(
    X = c(
      "datasets", "utils", "grDevices", "graphics", "stats", "methods",
      "parallel", "scales", "readxl", "writexl", "tidyverse"
    ),
    FUN = library,
    character.only = TRUE
  )
)

#### Functions
source(paste0(SourceDir, "/MiST.R"))
source(paste0(SourceDir, "/Parameters-class.R"))
source(paste0(SourceDir, "/QC-class.R"))
source(paste0(SourceDir, "/Summary_clinique.R"))
source(paste0(SourceDir, "/Zoom_genotype.R"))
source(paste0(SourceDir, "/Effectif_mutation_YCC.R"))
source(paste0(SourceDir, "/Effectif_mutation_Yquanti.R"))
source(paste0(SourceDir, "/Define_QC_Object.R"))
source(paste0(SourceDir, "/DO_QC.R"))

### LARGE or STRICT analysis
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


#### Create group of genes ####
liste_TR <- read.table(
  file = paste0(RawDir, "Docs/Latest_info_gene.csv"),
  header = TRUE,
  sep = ";"
) %>% 
  filter(get(GroupeGene)==1) %>% 
  select(gene_name, Transcript) %>% 
  filter(!gene_name %in% c("BOLA2", "BOLA2B", "DTX2P1-UPK3BP1-PMS2P11", "RMST")) 

#### MUST EXIST ONLY ONE ENST BY GENE ####
if (length(unique(liste_TR$gene_name)) != length(unique(liste_TR$Transcript))) {
  stop(paste("[MULTI-GENE]", "MUST EXIST ONLY ONE ENST BY GENE, else, this script need to be adapted (about files names, the loading data part, etc..."))
}


## 1/ tous les gènes
Group1 <- liste_TR$gene_name
Group1_final <- Group1
## 2/  les gènes "autosomiques dominants".
Group2 <- "ABCC8, APPL1, BLK, CEL, GATA4, GATA6, GCK, HNF1A, HNF1B, HNF4A, INS, KCNJ11, KLF11, NEUROD1, PAX4, PDX1, RFX6, STAT3, WFS1" %>% 
  strsplit(., ", ") %>%
  unlist()
Group2_final <- Group2
## 3/ les mêmes gènes en enlevant les gènes "bordeline (pour lesquels la communauté à de gros doutes)" 
## ainsi que KCNJ11 et ABCC8 (car ils peuvent aussi engendrer une hypoglycémie).
Group3 <- "APPL1, GATA4, GATA6, GCK, HNF1A, HNF1B, HNF4A, INS, NEUROD1, PDX1, RFX6, STAT3, WFS1" %>%
  strsplit(., ", ") %>%
  unlist()
Group3_final <- Group3
## 4/ les gènes "autosomiques recessifs".
Group4 <- "DNAJC3, FOXP3, GLIS3, IER3IP1, MNX1, NEUROG3, NKX2-2, PAX6, PCBD1, PPP1R15B, PTF1A, SLC19A2, SLC2A2, TRMT10A" %>%
  strsplit(., ", ") %>%
  unlist()
Group4_final <- Group4

## 5/  recalculer l'association avec le DT2 (STRICT + ACMG) rassemblant les gènes "officiels" du MODY [analyse multi-genes CC_T2D61]
## que sont : GCK, HNF1A, HNF4A, PDX1, CEL, PAX4, INS, NEUROD1 et HNF1B (et laisser les autres gènes de côté) 
Group5 <- "GCK, HNF1A, HNF4A, PDX1, CEL, PAX4, INS, NEUROD1, HNF1B" %>%
  strsplit(., ", ") %>%
  unlist()
Group5_final <- Group5

## 6/ 8 gene of interest
Genes_Variants03 <- readr::read_rds(
  path = paste0(
    gsub(script_name, "03_Replic_HealthyNevada", output_directory), "/Genes_Variants.rds")
)
##  8 genes seen in replication HNP 
Group6 <- sort(names(Genes_Variants03))
Group6
length(Group6)
Group6_final <- Group6

## 7/  7 gene of interest (sans than group6 without HNF1B)
Genes_Variants02 <- readr::read_rds(
  path = paste0(
    gsub(script_name, "02_Replic_UKBioBank", output_directory), "/Genes_Variants.rds")
)
##  7 genes seen in replication UK 
Group7 <- sort(names(Genes_Variants02))
Group7
length(Group7)
Group7_final <- Group7

#### SAVE QCed DATA for all genes in this group ####
allgenes <- unique(c(Group1, Group2, Group3, Group4, Group5, Group6, Group7)) 


#### init ####
everything_finish <- tibble()


#### Do the analyses, on Groups (Super-Clusters) ####
## classique Do the 3 analyses
if (VERSION == "CLASSIQUE" ) {
 analyses_to_do <- c("STRICT", "ACMG", "CUSTOM")  
}
if (VERSION == "ACMGCUSTOM") {
  analyses_to_do <- c("STRICT") 
}

for (analyse_LARGE_STRICT in analyses_to_do) {
  message(paste("[MULTI-GENE]", ">>>", analyse_LARGE_STRICT))

  #### prepare all genotype (QC gene by gene) ####
  for (gene in allgenes) {
    message(paste("[MULTI-GENE]", ">>> gene ongoing:", gene))

    #### ALL parameters ####
    Params <- new.MODY_Parameters()

    Params["FolderOut"] <- "NULL"
    Params["GroupeGene"] <- GroupeGene
    Params["Gene_Symbol"] <- gene
    Params["ENST"] <- liste_TR[grep(gene, liste_TR$gene_name), "Transcript"]
    Params["analyse_LARGE_STRICT"] <- analyse_LARGE_STRICT
    Params["threshold"] <- 0.01
    Params["ExludeSynonymous"] <- FALSE
    Params["cluster"] <- "none"
    Params["y_name"] <- y_name
    Params["binary"] <- binary
    Params["SNP_INDEL"] <- SNP_INDEL
    
    #### Phenotypes ####
    Phenotypes <- read_excel(
      path = Phenotype_file,
      col_names = TRUE,
      na = "NA",
      guess_max = 9365
    ) %>%
      as.data.frame()
    rownames(Phenotypes) <- paste(Phenotypes$RUN, Phenotypes$ID, sep = "_")
  
    #### LOAD RAW DATA from RData ####
    dataset_file <- paste0(
      InputDir,
      Params["Gene_Symbol"], "/", Params["Gene_Symbol"], "_", Params["ENST"], ".RData"
    )
    if (!file.exists(dataset_file)) {
      stop(paste('[MULTI-GENE]', 'The data does not exists in:\n      ', dataset_file, ''))
    }
    load(dataset_file)
    all_genotype_data <- all_genotype_data[Phenotypes$IID, ]
    
    #### Perform QCs ####
    if (ncol(all_genotype_data) == 0) {
      stop(paste("[MULTI-GENE]", "No genotype data. Consider removing", Params["Gene_Symbol"], "_", Params["ENST"], "from the ongoing group."))
    }
    #### KEEP ONLY SAMPLES INVOLVED IN THE STUDY ####
    message(paste("[MULTI-GENE]", nrow(Phenotypes), "samples at the begining."))

    all_genotype_data$IID <- rownames(all_genotype_data)
 
    QC_PHENO <- Phenotypes %>%
      select("IID", Params["y_name"]) %>%
      `colnames<-`(c("IID", "y_name")) %>%
      filter(!is.na(y_name))

    ## Remove sample
    Phenotypes <- Phenotypes %>%
      filter(IID %in% QC_PHENO$IID)
    
    all_genotype_data <- all_genotype_data %>%
      filter(IID %in% QC_PHENO$IID)

    message(paste("[MULTI-GENE]", nrow(Phenotypes), "samples involved in this study."))
    
    ## QC00 : remove all frequent variant
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
    
    #### QC steps ####
    ## define QC object
    QC <- QC5

    #### DO 1 to 5 QC gene by gene ####
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
      qc_type = "first5", 
      GroupeGene = GroupeGene
    )
    QC <- res_QC$QC
    Phenotypes <- res_QC$Phenotypes
    all_genotype_data <- res_QC$all_genotype_data
    all_annotation_data <- res_QC$all_annotation_data
    Detail_Geno <- res_QC$Detail_Geno
    all_quality_data <- res_QC$all_quality_data

    if (ncol(all_genotype_data) == 1 | nrow(Phenotypes) == 0) {
      message(paste("[MULTI-GENE]", paste(Params["Gene_Symbol"], "is removed from all group with this analyses.. No rds saved !!! !!! ")))

      if (Params["Gene_Symbol"] %in% Group1_final) {
        Group1_final <- Group1_final[-grep(Params["Gene_Symbol"], Group1_final)]
      }
      if (Params["Gene_Symbol"] %in% Group2_final) {
        Group2_final <- Group2_final[-grep(Params["Gene_Symbol"], Group2_final)]
      }
      if (Params["Gene_Symbol"] %in% Group3_final) {
        Group3_final <- Group3_final[-grep(Params["Gene_Symbol"], Group3_final)]
      }
      if (Params["Gene_Symbol"] %in% Group4_final) {
        Group4_final <- Group4_final[-grep(Params["Gene_Symbol"], Group4_final)]
      }
      if (Params["Gene_Symbol"] %in% Group5_final) {
        Group5_final <- Group5_final[-grep(Params["Gene_Symbol"], Group5_final)]
      }
      if (Params["Gene_Symbol"] %in% Group6_final) {
        Group6_final <- Group6_final[-grep(Params["Gene_Symbol"], Group6_final)]
      }
      if (Params["Gene_Symbol"] %in% Group7_final) {
        Group7_final <- Group7_final[-grep(Params["Gene_Symbol"], Group7_final)]
      }
    } else {
      #### write RDS 1 gene qc'ed 1 to 5 ####
      write_rds(
        x = all_genotype_data,
        path = paste0(outputPath, "/", Params["Gene_Symbol"], "_all_genotype_data.RDS")
      )
      write_rds(
        x = Detail_Geno,
        path = paste0(outputPath, "/", Params["Gene_Symbol"], "_Detail_Geno.RDS")
      )
      write_rds(
        x = all_annotation_data,
        path = paste0(outputPath, "/", Params["Gene_Symbol"], "_all_annotation_data.RDS")
      )
      message(paste("[MULTI-GENE]", "3 RDS saved. END of gene QC "))
    }
  }
  rm(Params)
  rm(gene)

  #### all genes are qc'ed 1 to 5
  #### FOR all groups, for this analysis, read all rds of genes and finish the QC 6 to 11 ####
  message(paste("[MULTI-GENE]", "read all rds of genes."))

  everything_old <- tibble(
    Analysis = rep(analyse_LARGE_STRICT, 7),  
    group = c("Group1", "Group2", "Group3", "Group4", "Group5", "Group6", "Group7"),
    genes = list(Group1, Group2, Group3, Group4, Group5, Group6, Group7),
    genes_final = list(Group1_final, Group2_final, Group3_final, Group4_final, Group5_final, Group6_final, Group7_final)
  ) %>% 
    mutate(
      Genotypes = map(.x = genes_final, .f = function(igenes) {
        map(.x = igenes, .f = function(igene) {
          file_gene <- paste0(outputPath, "/", igene, "_all_genotype_data.RDS")
          read_rds(path = file_gene)
        }) %>%
          Reduce(f = function(x, y) {inner_join(x, y, by = "IID")}, x = .) 
      }),
      Detail_Genos = map(.x = genes_final, .f = function(igenes) {
        map_df(.x = igenes, .f = function(igene) {
          file_gene <- paste0(outputPath, "/", igene, "_Detail_Geno.RDS")
          mutate(.data = read_rds(path = file_gene), Gene = igene)
        })
      }),
      All_annos = map(.x = genes_final, .f = function(igenes) {
        map_df(.x = igenes, .f = function(igene) {
          file_gene <- paste0(outputPath, "/", igene, "_all_annotation_data.RDS")
          mutate(.data = read_rds(path = file_gene), Gene = igene)
        })
      }),
      Phenos = map(.x = Genotypes, .f = function(MatGenotypes) {
        Phenotypes <- read_excel(
          path = !!Phenotype_file,
          col_names = TRUE,
          na = "NA",
          guess_max = 9365
        ) %>%
          as.data.frame()
        rownames(Phenotypes) <- paste(Phenotypes$RUN, Phenotypes$ID, sep = "_")
        Phenotypes <- Phenotypes %>%
          filter(IID %in% MatGenotypes$IID)
        
        return(Phenotypes)
    
      })
    )
  message(paste("[MULTI-GENE]", "everything_old ok."))
  
  QC <- QC11

  everything_old$QC <- list(QC)

  message(paste("[MULTI-GENE]", "map DO_QC_clusterMultiGene"))
  DO_TheEndOfQCs <- pmap(
    .l = list(
      everything_old$group,
      everything_old$QC, 
      everything_old$Phenos, 
      everything_old$Genotypes, 
      everything_old$All_annos, 
      everything_old$Detail_Genos, 
      CustomFile = CustomFile,
      y_name = y_name,
      thresholdrare = 0.01
    ),
    .f = DO_QC_clusterMultiGene
  )

  message(paste("[MULTI-GENE]", "DO_TheEndOfQCs on the whole cluster."))
  everything_qcs <- DO_TheEndOfQCs %>%
    transpose() %>%
    as_tibble() %>%
    mutate(
      Analysis = analyse_LARGE_STRICT,
      genes_init = list(Group1, Group2, Group3, Group4, Group5, Group6, Group7),
      genes_final = list(Group1_final, Group2_final, Group3_final, Group4_final, Group5_final, Group6_final, Group7_final)
    )

  #### MiST analyses with cluster based on Detail_Genos$Gene
  message(paste("[MULTI-GENE]", "MIST analyses."))
  everything_qcs <- everything_qcs %>%
    mutate(
      stat_MiST = list(c(), c(), c(), c(), c(), c(), c()), 
      stat_rare = list(c(), c(), c(), c(), c(), c(), c())
    )
  if (VERSION == "ACMGCUSTOM") {
    everything_qcs$QC12 <- list(c(), c(), c(), c(), c(), c(), c())
  }
  
  

  for (igroup in 1:nrow(everything_qcs)) {
    message("GROUP ", igroup)
    Genotype <- everything_qcs$all_genotype_data[[igroup]]
    Detail_Geno <- everything_qcs$Detail_Geno[[igroup]]
    Phenotypes <- everything_qcs$Phenotypes[[igroup]]

    if (nrow(Genotype) != nrow(Phenotypes)) {
      stop(paste("[MULTI-GENE]", "//!\\ in MIST analyses: nrow(Genotype)!=nrow(Phenotypes)"))
    }
    
    message(VERSION)
    if (VERSION == "ACMGCUSTOM") {
      message("QC12 based on custom file")
      
      QC12_customfile <- read_excel(
        path = ".//Docs/CUSTOM_files/custom_22072020.xlsx",  
        sheet = 1,
        col_names = TRUE, 
        col_types = "text"
      )
      QC12_excluded <- setdiff(x = Detail_Geno$Variant, y = QC12_customfile$CUSTOM)
      Detail_Geno <- Detail_Geno %>% 
        filter(Variant %in% QC12_customfile$CUSTOM) 
      Genotype <- Genotype[, c("IID", Detail_Geno$Variant), drop = FALSE]
      
      ## save it 
      everything_qcs$all_genotype_data[[igroup]] <- Genotype 
      everything_qcs$Detail_Geno[[igroup]] <- Detail_Geno  
    } 
    
    variants <- names(Genotype)
    variants <- variants[-grep("IID", variants)]
    
    if (VERSION == "ACMGCUSTOM") {
      everything_qcs$QC12[[igroup]] <- QC12_excluded 
    }
    
    alldta <- inner_join(Genotype, Phenotypes, by = "IID")

    G_rare <- alldta[, variants, drop = FALSE]
    if (ncol(G_rare) > 1) {
      message(paste("[MULTI-GENE]", "After all QCs,", ncol(G_rare), "rare variants will be analysed."))
      info_filtred <- Detail_Geno[Detail_Geno$Variant %in% names(G_rare), ]
      Z <- model.matrix(~ Gene - 1, data = info_filtred)

      out <- try(
        mist(
          y = Phenotypes[, y_name], 
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
        rownames(stat_rare) <- gsub("^GZtype", "", rownames(stat_rare))
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
    } else {
      message(paste("[MULTI-GENE]", "rare analysis impossible"))
      stat_rare <- NA
      stat_MiST <- NA
    }
    
    stat_MiST$trait = y_name
    stat_MiST$covariates = paste(covar_names, collapse = ";")
    stat_MiST$sample_size = nrow(alldta)

    everything_qcs$stat_MiST[[igroup]] <- stat_MiST
    everything_qcs$stat_rare[[igroup]] <- stat_rare

    rownames(everything_qcs$all_genotype_data[[igroup]]) <- everything_qcs$all_genotype_data[[igroup]]$IID
    rownames(everything_qcs$Phenotypes[[igroup]]) <- everything_qcs$Phenotypes[[igroup]]$IID
  }

  everything_qcs$Y_name <- y_name
  everything_qcs$covar_names <- list(covar_names)
  
  message(paste("[MULTI-GENE]", "everything_finish bind_rows."))
  everything_finish <- bind_rows(everything_finish, everything_qcs)

  message(paste("[MULTI-GENE]", "Delete gene .RDS"))
  file.remove(list.files(path = outputPath, pattern = "*.RDS$", full.names = TRUE))
}

everything_finish$gene_final <- map(everything_finish$Detail_Geno, .f = function(df){
  unique(df$Gene)
})


#### clinicalData ####
message(paste("[MULTI-GENE]", "Compute Clinical table."))
everything_finish$clinicalData <- map(
  .x = everything_finish$Phenotypes,
  .f = Summary_clinique,
  x = y_name, 
  na_symbol = "NA"
)

#### Effectif_mutation ####
message(paste("[MULTI-GENE]", "Count mutations carriers."))
everything_finish$Effectif_mutation <- map2(
  .x = everything_finish$Phenotypes,
  .y = everything_finish$all_genotype_data,
  .f = ~select(
    .data = Effectif_mutation_YCC(Phenotypes = .x, G_mat = .y, y_name = y_name),
    -starts_with("genotype_Missing")
  ) %>% mutate(Variant = as.character(Variant))
)

#### MERGE DETAIL and EFFECTIF to write in xslx ####
message(paste("[MULTI-GENE]", "Format QC annonation."))
everything_finish$Annotation_qced <- map2(
  .x = everything_finish$Detail_Geno,
  .y = everything_finish$Effectif_mutation,
  .f = ~full_join(x = .x, y = .y, by = "Variant")
)

everything_finish$Annotation_qced_maj <- map2(
  .x = everything_finish$Annotation_qced,
  .y = everything_finish$all_genotype_data,
  .f = function(Annotation_qced, all_genotype_data) {
    Annotation_qced <- Annotation_qced %>%
      mutate(
        Samples_with_mutation = map_chr(
          .x = Variant,
          .f = function(x) {
            all_genotype_data[, c(x, "IID"), drop = FALSE] %>%
              filter(get(x) != 0) %>%
              .[["IID"]] %>%
              paste(., collapse = ";")
          }
        )
      ) %>%
      select(Variant, Variant_Old, everything())
    ## For frequent variant Samples_with_mutation is too long for excel (limit of 32,767 characters)
    Annotation_qced[Annotation_qced$MAF > 0.05, "Samples_with_mutation"] <- "..."

    return(Annotation_qced)
  }
)


#### SAVE everything_finish ####
message(paste("[MULTI-GENE]", "save everything_finish.RDS"))
write_rds(x = everything_finish, path = paste0(outputPath, "/", "everything_finish.RDS"))


message(paste("[MULTI-GENE]", "write_xlsx Annotation_qced (with mutates by default)."))


f3_adapt <- function(listf2) {
  nrowQCtable <- max(unlist(lapply(listf2, length)))
  out <- data.frame(matrix(nrow = nrowQCtable, ncol = length(listf2)))
  names(out) <- names(listf2)
  for (i_el in 1:length(listf2)) {
    el <- unlist(listf2[i_el], use.names = FALSE)
    if (length(el) < nrowQCtable) {
      out[, i_el] <- c(el, rep(NA, nrowQCtable - length(el)))
    } else {
      out[, i_el] <- el
    }
  }
  return(out)
}


for (iline in 1:nrow(everything_finish)) {
  
  listf2 <- everything_finish$QC[[iline]] %>% 
    unlist(recursive = FALSE) %>% 
    lapply(X = ., FUN = function(.x) {.x["Excluded"]})

  names(listf2) <- paste0("QC_",1:length(everything_finish$QC[[iline]]))

  listf2 <- listf2[!is.na(listf2)]
  listf2 <- listf2[unlist(lapply(listf2, length))!=0] 
  listf2 <- listf2[!unlist(
    lapply(X = listf2, FUN = `%in%` , "NA") %>%
    lapply(X = ., FUN = any)
  )]
  
  if (VERSION %in% "ACMGCUSTOM" ) {
    listf2$QC_12 <- everything_finish$QC12[[iline]]
  }

  QC_table <- f3_adapt(listf2 = listf2)
  
  write_xlsx(
    x = list(
      "Annot" = everything_finish$Annotation_qced_maj[[iline]],
      "QCs" = QC_table
    ),
    path = paste0(outputPath, "/", everything_finish$Analysis[[iline]], "_", everything_finish$group[[iline]], ".xlsx"),
    col_names = TRUE
  )
}

