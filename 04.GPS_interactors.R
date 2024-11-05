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
gps_path <-"/Users/dandantan/Desktop/IntAct_spark_dataset/GPS_allgenes.csv"
gps_results_all <- read.csv(gps_path)
gps_results_all$gene_trait_pairs <- paste(gps_results_all$Gene, gps_results_all$Phecode.description, sep = "_")


threshold <- quantile(gps_results_all$Genetic.priority.score..GPS., 0.9981, na.rm = TRUE)
# Select the top 0.19% of rows based on the Genetic priority score
gps_results_top_0.19 <- gps_results_all[gps_results_all$Genetic.priority.score..GPS. > threshold, ]
#195 proteins selected
gps_results_unselect<-gps_results_all[gps_results_all$Genetic.priority.score..GPS. <= threshold, ]

gps_results_top_0.19$is_true<-ifelse(gps_results_top_0.19$gene_trait_pairs %in% merged_positive_gene_list$gene_trait_pairs,1,0)
gps_results_unselect$is_true<-ifelse(gps_results_unselect$gene_trait_pairs %in% merged_positive_gene_list$gene_trait_pairs,1,0)

TP_GPS<-sum(gps_results_top_0.19$is_true)
FP_GPS<-nrow(gps_results_top_0.19)-sum(gps_results_top_0.19$is_true)

FN_GPS<-sum(gps_results_unselect$is_true)
TN_GPS<-nrow(gps_results_unselect)-sum(gps_results_unselect$is_true)

#####Add interactor
ensembl_gene_id<-read.csv( "/Users/dandantan/Desktop/IntAct_spark_dataset/ensembl_gene_GPS.csv")

ensembl_gene_id[18,2]="ENSG00000084674"
ensembl_gene_id[27,2]="ENSG00000084674"
ensembl_gene_id[43,2]="ENSG00000197249"
ensembl_gene_id[44,2]="ENSG00000243649"
ensembl_gene_id[53,2]="ENSG00000197249"
ensembl_gene_id[56,2]="ENSG00000104044"
ensembl_gene_id[67,2]="ENSG00000177628" #GBA1
ensembl_gene_id[69,2]="ENSG00000084674"
ensembl_gene_id[73,2]="ENSG00000175164"
ensembl_gene_id[76,2]="ENSG00000206449"
ensembl_gene_id[103,2]="ENSG00000275410"
ensembl_gene_id[106,2]="ENSG00000177628" #GBA1
ensembl_gene_id[139,2]="ENSG00000204983"
ensembl_gene_id[149,2]="ENSG00000204267"
ensembl_gene_id[153,2]="ENSG00000084674"
ensembl_gene_id[193,2]="ENSG00000161011"
gps_results_top_0.19_with_geneid <- cbind(gps_results_top_0.19,ensembl_gene_id)


config <- spark_config()
sc <- spark_connect(master = "local", config = config)

local_path <- "/Users/dandantan/Desktop/Open_target_2021_FDA/Dataset"
interaction_path <- paste(
  local_path,
  "/interaction/",
  sep = ""
)
database<-"intact" #functional/physical
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

gps_results_top_0.19_with_geneid_spark <- copy_to(sc, gps_results_top_0.19_with_geneid, "gps_results_with_geneid_spark",overwrite = TRUE)

#Filter by interactions only
interactors_ass <- gps_results_top_0.19_with_geneid_spark %>%
  inner_join(interactions, by = c("ensmble_genes" = "targetA")) %>%
  #dplyr::select(X_1,X,all_locus_prob, all_locus_y,names_genes,all_trait,locus_chrm,locus_start,locus_end, locus_name, gene_trait_pairs,targetB) %>%
  sdf_distinct() %>%
  collect() #

interactors_ass$gene_trait_pairs2 <- paste(interactors_ass$targetB, interactors_ass$Phecode_description,sep = "_")
interactors_ass$is_in_positive <- ifelse(interactors_ass$gene_trait_pairs2 %in% merged_positive_gene_list$gene_tait_pair_merged, 1, 0)
sum(interactors_ass$is_in_positive)

gps_results_top_0.19_with_geneid$gene_trait_pair_hsgn<-paste(gps_results_top_0.19_with_geneid$ensmble_genes,gps_results_top_0.19_with_geneid$Phecode.description,sep="_")
interactors_ass$original_GPS_gene<-ifelse(interactors_ass$gene_trait_pairs2 %in% gps_results_top_0.19_with_geneid$gene_trait_pair_hsgn, 1, 0) 

interactors_only <- interactors_ass %>%
  filter(original_GPS_gene != 1)

TP_gps_interactor<-sum(interactors_only$is_in_positive)
FP_gps_interactor<-nrow(interactors_only)-sum(interactors_only$is_in_positive)

TP_gps_with_interactor<-TP_GPS+TP_gps_interactor
FP_gps_with_interactor<-FP_GPS+FP_gps_interactor
TN_gps_with_interactor<-TN_GPS
FN_gps_with_interactor<-FN_GPS

#########Calculate sensitivity, specificity and precision
sensitivity_gps<-(TP_GPS/(TP_GPS+FN_GPS))
precision_gps<-(TP_GPS / (TP_GPS + FP_GPS))
specificity_gps<- (TN_GPS / (TN_GPS + FP_GPS))


sensitivity_gps_with_interactor<-(TP_gps_with_interactor/(TP_gps_with_interactor+FN_gps_with_interactor))
precision_gps_with_interactor<-(TP_gps_with_interactor / (TP_gps_with_interactor + FP_gps_with_interactor))
specificity_gps_with_interactor<- (TN_gps_with_interactor / (TN_gps_with_interactor + FP_gps_with_interactor))


precision_gps_with_interactor
sensitivity_gps_with_interactor
specificity_gps_with_interactor


