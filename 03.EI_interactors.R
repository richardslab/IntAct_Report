library("tidyverse")
library("sparklyr")
library("sparklyr.nested")
library("ggsci")
library("ggplot2")
library("reshape2")
library("dplyr")
library("biomaRt")

# Read the positive control set
merged_positive_gene_list<-read.csv("/Users/dandantan/Desktop/IntAct_spark_dataset/merged_positive_control_gene_list.csv")
merged_positive_gene_list$gene_trait_pairs <- paste(
  merged_positive_gene_list$gene,
  ifelse(is.na(merged_positive_gene_list$trait_exwas_ei), 
         merged_positive_gene_list$trait_gps, 
         merged_positive_gene_list$trait_exwas_ei),
  sep = "_"
)
###########Data Preprocesses
# Read the EI algorithm prediction results
ei_path <-"/Users/dandantan/Desktop/IntAct_spark_dataset/ei_results_new.csv"
ei_results_all <- read.csv(ei_path)
ei_results_all$gene_trait_pairs <- paste(ei_results_all$names.genes, ei_results_all$all.trait, sep = "_")

ei_results_all$is_in_positive <- ifelse(ei_results_all$gene_trait_pairs %in% merged_positive_gene_list$gene_trait_pairs, 1, 0)

# Find all loci containing positive control genes for each trait
positive_loci <- ei_results_all %>%
  group_by(all.trait) %>%
  filter(is_in_positive == 1) %>%
  dplyr::select(all.trait, locus.name)

# Only keep the loci with at least one positive control gene for each trait
ei_loci <- ei_results_all %>%
  semi_join(positive_loci, by = c("all.trait","locus.name"))

# Find the highest EI score for each trait at each locus
highest_prob_per_locus <- ei_loci %>%
  group_by(all.trait, locus.name) %>%
  filter(all.locus.prob == max(all.locus.prob, na.rm = TRUE))

ei_unselect <- ei_results_all %>%
  anti_join(highest_prob_per_locus, by = "gene_trait_pairs")

# Convert Gene Symbol to ensembl ID
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
ensembl_gene_id <- data.frame(gene_name = highest_prob_per_locus$names.genes,
                              ensmble_genes = NA)

for (i in 1:nrow(ensembl_gene_id)) {
  gene_name <- ensembl_gene_id$gene_name[i]
  ensembl_id <- getBM(attributes = 'ensembl_gene_id', 
                      filters = 'hgnc_symbol', 
                      values = gene_name, 
                      mart = ensembl)
  if (length(ensembl_id) == 0) {
    ensembl_gene_id$ensmble_genes[i] <- NA}
  else{
    ensembl_gene_id$ensmble_genes[i] <- ensembl_id
  }
}
# manually correct ensbmle id
ensembl_gene_id[[24,2]] = "ENSG00000225663"
ensembl_gene_id[[44,2]] = "ENSG00000175164"
ensembl_gene_id[[74,2]] = "ENSG00000171862"
ensembl_gene_id[[84,2]] = "ENSG00000103489"
ensembl_gene_id[[100,2]] = "ENSG00000214575"
ensembl_gene_id[[101,2]] = "ENSG00000007968"
ensembl_gene_id[[137,2]] = "ENSG00000285346"
ensembl_gene_id[[145,2]] = "ENSG00000275410"
ensembl_gene_id[[150,2]] = "ENSG00000053918"
ensembl_gene_id[[168,2]] = "ENSG00000084674"

highest_prob_per_locus_new <- highest_prob_per_locus %>% 
  left_join(ensembl_gene_id, by = c("names.genes" = "gene_name"))

highest_prob_per_locus_new$ensmble_genes <- unlist(highest_prob_per_locus_new$ensmble_genes)
#######################--------------------------Add IntAct Information -------------------------#############################

### Find IntAct Partners
config <- spark_config()
sc <- spark_connect(master = "local", config = config)

local_path <- "/Users/dandantan/Desktop/Open_target_2021_FDA/Dataset"
interaction_path <- paste(
  local_path,
  "/interaction/",
  sep = ""
)
database<-"string" #functional/physical

# Data about molecular interactions
interactions <- spark_read_parquet(sc, interaction_path, memory = FALSE) %>%
  filter(sourceDatabase == database) %>%
  filter(!is.na(targetA)) %>%
  filter(!is.na(targetB)) %>%
  filter(scoring > 0.42) %>%
  filter(speciesA == speciesB)%>%  # only keep human species
  dplyr::select(targetA, targetB) %>%
  sdf_distinct()
#local_interaction <- as.data.frame(interactions %>% collect()) #view as dataframe

highest_prob_per_locus_spark <- copy_to(sc, highest_prob_per_locus_new, "highest_prob_per_locus_spark",overwrite = TRUE)

#Filter by interactions only
interactors_ass <- highest_prob_per_locus_spark %>%
  inner_join(interactions, by = c("ensmble_genes" = "targetA")) %>%
  dplyr::select(X_1,X,all_locus_prob, all_locus_y,names_genes,all_trait,locus_chrm,locus_start,locus_end, locus_name, gene_trait_pairs,targetB) %>%
  sdf_distinct() %>%
  collect() #968 rows

interactors_ass$gene_trait_pairs2 <- paste(interactors_ass$targetB, interactors_ass$all_trait,sep = "_")

#check if the interacting genes in the positive control set
interactors_ass$is_in_positive <- ifelse(interactors_ass$gene_trait_pairs2 %in% merged_positive_gene_list$gene_tait_pair_merged, 1, 0)
interactors_ass <- interactors_ass %>% rename(
  X.1 = X_1,
  all.locus.prob = all_locus_prob,
  all.locus.y = all_locus_y,
  names.genes = names_genes,
  all.trait = all_trait,
  locus.chrm = locus_chrm,
  locus.start = locus_start,
  locus.end = locus_end,
  locus.name = locus_name,
)
sum(interactors_ass$is_in_positive) #49/1613

highest_prob_per_locus_new$gene_trait_pair_ensg<-paste(highest_prob_per_locus_new$ensmble_genes,highest_prob_per_locus_new$all.trait,sep="_")
interactors_ass$original_EI_gene<-ifelse(interactors_ass$gene_trait_pairs2 %in% highest_prob_per_locus_new$gene_trait_pair_ensg, 1, 0) 

interactors_only <- interactors_ass %>%
  filter(original_EI_gene != 1)

TP_ei_interactor<-sum(interactors_only$is_in_positive)
FP_ei_interactor<-nrow(interactors_only)-sum(interactors_only$is_in_positive)
########################------------------Evaluate  with sensitivity and specificity and precision-------------------------##################################

### EI prediction
TP_ei<- sum(highest_prob_per_locus_new$is_in_positive) 
FP_ei <- nrow(highest_prob_per_locus_new)-sum(highest_prob_per_locus_new$is_in_positive)  
FN_ei<-sum(ei_unselect$is_in_positive)
TN_ei<-nrow(ei_unselect)-sum(ei_unselect$is_in_positive)

### EI prediction + IntAct
TP_ei_with_interactor<- TP_ei+TP_ei_interactor
FP_ei_with_interactor <- FP_ei+FP_ei_interactor
FN_ei_with_interactor<- FN_ei 
TN_ei_with_interactor <- TN_ei 


sensitivity_ei<- TP_ei/(TP_ei+FN_ei) 
specificity_ei<- TN_ei/(TN_ei+FP_ei) 
precision_ei<-TP_ei/(FP_ei+TP_ei) 

sensitivity_ei_intact<- TP_ei_with_interactor/(TP_ei_with_interactor+FN_ei_with_interactor) 
specificity_ei_intact<- TN_ei_with_interactor/(TN_ei_with_interactor+FP_ei_with_interactor) 
precision_ei_intact<-TP_ei_with_interactor/(FP_ei_with_interactor+TP_ei_with_interactor) 




