## Do a meta analyses for gene present in réplication... ##

options(stringsAsFactors = FALSE)

script_name <- "07_Metaanalysis_inMODY"
working_directory <-  ".//Data"
output_directory <- paste(working_directory, script_name, sep = "/")
dir.create(path = output_directory, recursive = TRUE, showWarnings = FALSE, mode = "0777")

## load meta & grid
library(grid)
library(meta)
library(readxl) 
library(writexl) 
library(officer) ## for read_pptx => write pptx 
library(rvg) ## for ph_with_vg_at
library(tidyverse) ## to use `%>%`, mutate, select ...

#### DATASET ---------------------------------------------------------------------------------------
MODY_res <- readxl::read_xlsx(
  path = ".//Data/01_ACMG_CUSTOM/000-RESULTS/CUSTOM_T2Dmono_CCDiabete_CC_T2D61_rare.xlsx",
  sheet = 1
) %>% 
  rename(
    "Gene_Name" = "Gene_Symbol", 
    "CI_2.5" = "CI_2_5", 
    "CI_97.5" = "CI_97_5"
  ) %>% 
  mutate(POP = "MODY")

HNP_res <- readxl::read_xlsx(
  path = ".//Data/03_Replic_HealthyNevada/MiST_results.xlsx",
  sheet = 1
) %>% 
  mutate(POP = "HNP")

UK_res <- readxl::read_xlsx(
  path = ".//Data/02_Replic_UKBioBank/MiST_results.xlsx",
  sheet = 1
) %>% 
  mutate(POP = "UK")

AMP_res <- readxl::read_xlsx(
  path = ".//Data/07_Metaanalysis_inMODY/AMPpourmeta.xlsx",
  sheet = 1
) %>% 
  mutate(nb_rare_var = NA)

Common_gene <- intersect(MODY_res$Gene_Name, HNP_res$Gene_Name) %>%
  intersect(.,  UK_res$Gene_Name)  %>%
  intersect(.,  AMP_res$Gene_Name) 
message("want to check those genes: ", paste(Common_gene, collapse = ', '))

## gather all dataset  ... 
dta <- rbind.data.frame(
  subset(
    x = MODY_res, 
    subset = Gene_Name %in% Common_gene, 
    select = c("POP", "Gene_Name", "Pi_hat", "CI_2.5", "CI_97.5", "SE", "P_val", "OR", "nb_rare_var")
  ),
  subset(
    x = HNP_res, 
    subset = Gene_Name %in% Common_gene, 
    select = c("POP", "Gene_Name", "Pi_hat", "CI_2.5", "CI_97.5", "SE", "P_val", "OR", "nb_rare_var")
  ), 
  subset(
    x = UK_res, 
    subset = Gene_Name %in% Common_gene, 
    select = c("POP", "Gene_Name", "Pi_hat", "CI_2.5", "CI_97.5", "SE", "P_val", "OR", "nb_rare_var")
  ), 
  subset(
    x = AMP_res, 
    subset = Gene_Name %in% Common_gene, 
    select = c("POP", "Gene_Name", "Pi_hat", "CI_2.5", "CI_97.5", "SE", "P_val", "OR", "nb_rare_var")
  )
)

dta
# dput(dta)
## pour info voila la structure de mon object "dta" : 
# dta <- structure(list(
#   POP = c("MODY", "MODY", "MODY", "MODY", "MODY",
# "MODY", "MODY", "HNP", "HNP", "HNP", "HNP", "HNP", "HNP", "HNP",
# "UK", "UK", "UK", "UK", "UK", "UK", "UK"),
# Gene_Name = c("ABCC8",
# "GATA4", "GATA6", "GCK", "HNF1A", "HNF4A", "KCNJ11", "GCK", "ABCC8",
# "HNF4A", "GATA6", "HNF1A", "KCNJ11", "GATA4", "ABCC8", "GCK",
# "HNF4A", "KCNJ11", "HNF1A", "GATA4", "GATA6"),
#   Pi_hat = c(0.207925600198783,
# 1.88567818556954, -4.75992047369255, 1.95771721535179, 2.00145743780637,
# 0.124528656711712, 0.848809580481303, 1.59789998110591, 1.64855083051555,
# -12.8890923907403, 0.795758203217314, -0.181625934119027, 0.341094587458617,
# -10.473897226316, 1.33532111383067, 1.60732774438519, 0.952408769615033,
# 1.81088440774369, -1.25500902905987, 0.222207122098226, 0.143517675541615
# ),
#   CI_2.5 = c(-0.534048454761081, -0.601036782760681, -9.48657936302997,
# 0.802468153792509, 0.954530055954353, -1.26225416123246, -0.292510493937119,
# 0.477848745000833, 0.355359502910546, NA, -0.599772170653116,
# -1.25016750631846, -1.45734454639282, NA, 0.758935673070754,
# 0.777851529290706, 0.156776271206201, 0.287223470589621, -3.07086716199774,
# -2.71831815943812, -0.91354322082316),
#   CI_97.5 = c(0.909675613406493,
# 5.03863412188904, -0.990170289493387, 3.30670113474363, 3.27976643072443,
# 1.40632055899676, 2.00519731543251, 2.82409933081721, 3.05865416604063,
# 11.4381014246579, 2.08243933847273, 0.76076406306055, 2.142059012374,
# 29.6700378153818, 1.86550165793939, 2.34236439934898, 1.63511122911422,
# 3.05014771461879, -0.0806635423557972, 1.97194826680101, 0.964391305436382
# ), SE = c(0.365645379566775, 1.29917083078154, 2.81727715232674,
# 0.624226995071529, 0.578077881973835, 0.659777899226375, 0.575625391530438,
# 0.585692882919847, 0.67184114163164, 234.18057580129, 0.663807316790707,
# 0.504827709911626, 0.882461108612321, 196.967684763947, 0.280702039244763,
# 0.394530623372692, 0.372889953402636, 0.679451979926418, 0.72442343829915,
# 1.08044361712239, 0.469066126903179),
#   P_val = c(0.569591170272858,
# 0.146655305926441, 0.0911147675402399, 0.00171137003357648, 0.000535654215356047,
# 0.850294014467458, 0.140323676892844, 0.00636768127788481, 0.014136315991808,
# 0.956107299195403, 0.230613930200858, 0.719013111613838, 0.699106794118713,
# 0.957591905837756, 1.96416613280035e-06, 4.62074569510503e-05,
# 0.0106454033346911, 0.00769395203405579, 0.0831979804065563,
# 0.837054273395551, 0.759631525034932),
#   OR = c(1.2311215709413,
# 6.59082272754561, 0.00856629061587152, 7.08313931392432, 7.39983304009577,
# 1.13261447700174, 2.33686334774416, 4.94264187671269, 5.19943949583682,
# 2.52544722845588e-06, 2.2161206297635, 0.833913220563743, 1.40648627053251,
# 2.826469026029e-05, 3.80121637361168, 4.98946028395243, 2.59194554556803,
# 6.11585394968765, 0.285073274269172, 1.24883001137128, 1.1543272139855
# ),
#   nb_rare_var = c(33, 3, 2, 14, 19, 6, 13, 10, 7, 2, 5, 19,
# 4, 1, 17, 15, 4, 4, 36, 7, 5)),
#   row.names = c(NA, -21L), class = c("tbl_df", "tbl", "data.frame")
# )

#### META ------------------------------------------------------------------------------------------

## do the metaanalysis for each gene ... Go 

all_res <- tibble(
  Gene_Name = Common_gene
) %>% 
  mutate(
    meta_object = map(.x = Gene_Name, .f = function(igene) {
      meta::metagen(
        TE = dta$Pi_hat[dta$Gene_Name %in% igene], 
        seTE = dta$SE[dta$Gene_Name %in% igene], 
        studlab = dta$POP[dta$Gene_Name %in% igene], 
        data = dta[dta$Gene_Name %in% igene, ], 
        sm = "MD"
      )
    }), 
    base = map(.x = meta_object, .f = ~summary(.x)), 
    tmp = map2(.x = base, .y = meta_object, .f = function(base, meta_object){
      tmp <- unlist(base$fixed)
      # NB : utiliser la fonction random quand la pval d'hétérogénéité est inférieure à 0.05. ##
      if (meta_object$pval.Q < 0.05) {
        tmp		<- unlist(base$random)
      }
      return(tmp)
    }), 
    my_meta_res = map2(.x = tmp, .y = base, .f = function(tmp, base){
      tmp <- unlist(tmp)
      base <- unlist(base)
      
      beta <- tmp['TE']
      cil <- tmp['lower']
      ciu <- tmp['upper']
      p.value <- tmp['p']
      pv.het <- pchisq(q = base$Q, df = (base$k-1), lower.tail = FALSE)
      res <- c(beta, cil, ciu, p.value, pv.het)
	    names(res)	<- c("BETA", "L95", "U95", "p-value", "Heterogeneity")
      # res <- c("Beta" = beta, "cil" = cil, "ciu" = ciu, "p.value" = p.value, "pv.het" = pv.het)
      return(res)
    })
    
  )

## save it 
saveRDS(object = all_res, file = paste0(output_directory, "/all_res.RDS"))

## extract meta res ...
all_meta_res <- do.call(what = "rbind", all_res$my_meta_res %>% `names<-`(all_res$Gene_Name)) %>% 
  as.data.frame() %>% 
  rownames_to_column("Gene_Name")
## to write in xlsx ...
writexl::write_xlsx(x = all_meta_res, path = paste0(output_directory, "/all_meta_res.xlsx"), col_names = TRUE)


#### FOREST PLOT -----------------------------------------------------------------------------------

for (i in 1:length(Common_gene)) {

  ppt <- read_pptx()
  
  ppt <- ppt %>%
  add_slide(layout = "Blank", master = "Office Theme") %>%
  ph_with_vg_at(
    code = print(meta::forest(x = all_res$meta_object[[i]])), 
    width = 9.5, 
    height = 8, 
    left = 0.25, 
    top = 0.25
  )
  
  print(x = ppt, target = paste0(output_directory, "/forestplot_", all_res$Gene_Name[[i]], ".pptx"))
  
}

## write dta 
writexl::write_xlsx(x = dta %>% arrange(POP, Gene_Name), path = paste0(output_directory, "/detail_data_used_inMetaAnalyse.xlsx"), col_names = TRUE)
