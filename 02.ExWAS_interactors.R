library("tidyverse")
library("sparklyr")
library("sparklyr.nested")
library("ggsci")
library("ggplot2")


merged_positive_gene_list<-read.csv("/Users/dandantan/Desktop/IntAct_spark_dataset/merged_positive_control_gene_list.csv")


# Define the list of all_sig datasets and corresponding pathway datasets
all_sig_datasets <- list(
  all_sig_pLoF = "sig_gene_ids_pLoF(196)_new.csv",
  all_sig_pLoF_alpha = "sig_gene_ids_pLoF_alpha(341)_new.csv",
  all_sig_pLoF_5in5 = "sig_gene_ids_pLoF_5in5(317)_new.csv",
  all_sig_pLoF_5in5and1 = "sig_gene_ids_pLoF_5in5and1(407)_new.csv"
)

exwas_datasets<-list(
  exwas_pLoF="ExWAS_results_all_plof_new.csv",
  exwas_pLoF_alpha="ExWAS_results_all_alpha_new.csv",
  exwas_pLoF_5in5="ExWAS_results_all_5in5_new.csv",
  exwas_pLoF_5in5and1="ExWAS_results_all_5in5and1_new.csv"
)

TP_exwas_all<-c()
FP_exwas_all<-c()
FN_exwas_all<-c()
TN_exwas_all<-c()

TP_exwas_interactors_all<-c()
FP_exwas_interactors_all<-c()
FN_exwas_interactors_all<-c()
TN_exwas_interactors_all<-c()



ExWAS_methods <- c("pLoF","pLoF with \n AlphaMissense","pLoF with \n Missense (5/5)","pLoF with \n Missense (1/5)")

# Loop over the datasets
for (i in seq_along(all_sig_datasets)) {
  sig_genes_path<- paste0("/Users/dandantan/Desktop/IntAct_spark_dataset/", all_sig_datasets[[i]])
  sig_genes<-read.csv(sig_genes_path)
  sig_genes$gene_trait_pair_merged<-paste(sig_genes$gene_id,sig_genes$Traits,sep="_")
  sig_genes$is_true <- ifelse(sig_genes$gene_trait_pair_merged %in% merged_positive_gene_list$gene_tait_pair_merged, 1, 0) 
  
  TP_exwas<-sum(sig_genes$is_true)
  FP_exwas<-nrow(sig_genes)-sum(sig_genes$is_true)

  exwas_path<-paste0("/Users/dandantan/Desktop/IntAct_spark_dataset/", exwas_datasets[[i]])
  exwas_results<-read.csv(exwas_path)
  exwas_unselect <- exwas_results %>%
    anti_join(sig_genes, by = "trait_gene_pairs")
  
  exwas_unselect$gene_trait_pair_merged<-paste(exwas_unselect$gene,exwas_unselect$Traits,sep="_")
  exwas_unselect$is_true<-ifelse(exwas_unselect$gene_trait_pair_merged %in% merged_positive_gene_list$gene_tait_pair_merged, 1, 0) 
  FN_exwas<-sum(exwas_unselect$is_true)
  TN_exwas<-nrow(exwas_unselect)-sum(exwas_unselect$is_true)
  
  TP_exwas_all<-c(TP_exwas_all,TP_exwas)
  FP_exwas_all<-c(FP_exwas_all,FP_exwas)
  FN_exwas_all<-c(FN_exwas_all,FN_exwas)
  TN_exwas_all<-c(TN_exwas_all,TN_exwas)
  
  
  ###############ADD Interactors##########
  #Spark config
  config <- spark_config()
  # spark connect
  sc <- spark_connect(master = "local", config = config)
  
  #set significant gene ids in spark format
  sig_gene_ids <- spark_read_csv(
    sc,
    path = sig_genes_path,
    memory = FALSE
  )
  local_sig_gene_ids <- as.data.frame(sig_gene_ids %>% collect())
  
  ###### Read Platform data
  local_path <- "/Users/dandantan/Desktop/Open_target_2021_FDA/Dataset" # data released 21.11
  
  interaction_path <- paste(
    local_path,
    "/interaction/",
    sep = ""
  )
  
  database<-"intact" #functional/physical
  # Data about molecular interactions (Functional)
  interactions <- spark_read_parquet(sc, interaction_path, memory = FALSE) %>%
    filter(sourceDatabase == database) %>%
    filter(scoring > 0.42) %>%
    filter(!is.na(targetA)) %>%
    filter(!is.na(targetB)) %>%
    filter(speciesA == speciesB)%>%  # only keep human species
    dplyr::select(targetA, targetB) %>%
    sdf_distinct()
  
  #Filter by interactions
  interactors_ass <- sig_gene_ids %>%
    inner_join(interactions, by = c("gene_id" = "targetA")) %>%
    dplyr::select(trait_gene_pairs,Traits,gene_id,targetB) %>%
    sdf_distinct() %>%
    collect()
  #write.csv(interactors_ass, "~/Desktop/interactions.csv")
  
  # New gene_trait pairs with interacting genes
  interactors_ass$interactor_trait_pair <- paste(interactors_ass$targetB,interactors_ass$Traits,sep="_")
  interactors_ass$is_true<-ifelse(interactors_ass$interactor_trait_pair %in% merged_positive_gene_list$gene_tait_pair_merged, 1, 0) 
  interactors_ass$original_sig_gene<-ifelse(interactors_ass$interactor_trait_pair %in% sig_genes$gene_trait_pair_merged, 1, 0) 
  
  interactors_only <- interactors_ass %>%
    filter(original_sig_gene != 1)
  
  TP_exwas_interactor<-sum(interactors_only$is_true)
  FP_exwas_interactor<-nrow(interactors_only)-sum(interactors_only$is_true)
  
  TP_exwas_with_interactor<-TP_exwas_interactor+TP_exwas
  FP_exwas_with_interactor<-FP_exwas_interactor+TP_exwas
  
  TP_exwas_interactors_all<-c(TP_exwas_interactors_all,TP_exwas_with_interactor)
  FP_exwas_interactors_all<-c(FP_exwas_interactors_all,FP_exwas_with_interactor)
  TN_exwas_interactors_all<-c(TN_exwas_interactors_all,TN_exwas)
  FN_exwas_interactors_all<-c(FN_exwas_interactors_all,FN_exwas)
}

########Calculate sensitivity, specificity and precision
sensitivity_exwas<-(TP_exwas_all/(TP_exwas_all+FN_exwas_all))
precision_exwas<-(TP_exwas_all / (TP_exwas_all + FP_exwas_all))
specificity_exwas<- (TN_exwas_all / (TN_exwas_all + FP_exwas_all))


sensitivity_exwas_with_interactor<-(TP_exwas_interactors_all/(TP_exwas_interactors_all+FN_exwas_interactors_all))
precision_exwas_with_interactor<-(TP_exwas_interactors_all / (TP_exwas_interactors_all + FP_exwas_interactors_all))
specificity_exwas_with_interactor<- (TN_exwas_interactors_all / (TN_exwas_interactors_all + FP_exwas_interactors_all))

####Calculate averaged TP,FP,TN,FN
averaged_TP_exwas<-mean(TP_exwas_all)
averaged_FP_exwas<-mean(FP_exwas_all)
averaged_TN_exwas<-mean(TN_exwas_all)
averaged_FN_exwas<-mean(FN_exwas_all)

averaged_TP_exwas_with_interactor<-mean(TP_exwas_interactors_all)
averaged_FP_exwas_with_interactor<-mean(FP_exwas_interactors_all)
averaged_TN_exwas_with_interactor<-mean(TN_exwas_interactors_all)
averaged_FN_exwas_with_interactor<-mean(FN_exwas_interactors_all)

averaged_sensitivity_exwas<-(averaged_TP_exwas/(averaged_TP_exwas+averaged_FN_exwas))
averaged_precision_exwas<-(averaged_TP_exwas / (averaged_TP_exwas + averaged_FP_exwas))
averaged_specificity_exwas<- (averaged_TN_exwas / (averaged_TN_exwas + averaged_FP_exwas))

averaged_sensitivity_exwas_with_interactor<-(averaged_TP_exwas_with_interactor/(averaged_TP_exwas_with_interactor+averaged_FN_exwas_with_interactor))
averaged_precision_exwas_with_interactor<-(averaged_TP_exwas_with_interactor / (averaged_TP_exwas_with_interactor + averaged_FP_exwas_with_interactor))
averaged_specificity_exwas_with_interactor<- (averaged_TN_exwas_with_interactor / (averaged_TN_exwas_with_interactor + averaged_FP_exwas_with_interactor))

  
  

