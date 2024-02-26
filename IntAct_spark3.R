library("tidyverse")
library("sparklyr")
library("sparklyr.nested")
library("cowplot")
library("ggsci")
library("ggplot2")
library("likert")
library("reshape2")


# Define the list of all_sig datasets and corresponding pathway datasets
all_sig_datasets <- list(
  all_sig_pLoF = "sig_gene_ids_pLoF(196).csv",
  all_sig_pLoF_alpha = "sig_gene_ids_pLoF_alpha(341).csv",
  all_sig_pLoF_5in5 = "sig_gene_ids_pLoF_5in5(317).csv",
  all_sig_pLoF_5in5and1 = "sig_gene_ids_pLoF_5in5and1(407).csv"
)

pathway_datasets <- list(
  gs_sig_gene_ids_pLoF = "TP_gene_ids_pLoF(42).csv",
  gs_sig_gene_ids_pLoF_alpha = "TP_gene_ids_pLoF_alphamissense(60).csv",
  gs_sig_gene_ids_pLoF_5in5 = "TP_gene_ids_pLoF_5in5missense(57).csv",
  gs_sig_gene_ids_pLoF_5in5and1 = "TP_gene_ids_pLoF_5and1in5missense(60).csv"
)

# Initialize variables to store the results
True_positives_all <- vector("list", length = 0)
overlap_with_original_all <- vector("list", length = 0)
new_identified_all <- vector("list", length = 0)
overlap_False_positive_all <- vector("list", length = 0)
interactions_all <- vector("list", length = 0)
interactions_all2 <- vector("list", length = 0)
# Loop over the datasets
for (i in seq_along(all_sig_datasets)) {
  # all_sig and pathway datasets
  all_sig_pathway <- paste0("/Users/dandantan/Desktop/IntAct_spark_dataset/", all_sig_datasets[[i]])
  TP <- paste0("/Users/dandantan/Desktop/IntAct_spark_dataset/", pathway_datasets[[i]])
  
  #read csv
  all_sig <- read.csv(all_sig_pathway)
  
  #Spark config
  config <- spark_config()
  
  # spark connect
  sc <- spark_connect(master = "local", config = config)
  
  # read csv
  sig_gene_ids <- spark_read_csv(
    sc,
    path = TP,
    memory = FALSE
  )
  local_sig_gene_ids <- as.data.frame(sig_gene_ids %>% collect())
  ###### Read Platform data
  local_path <- "/Users/dandantan/Desktop/Open_target_2021_FDA/Dataset" # data released 21.11
  
  ass_indirectby_ds_path <- paste(
    local_path,
    "/associationByDatasourceIndirect/",
    sep = ""
  )
  
  interaction_path <- paste(
    local_path,
    "/interaction/",
    sep = ""
  )
  
  # Data about indirect associations
  ass_indirectby_ds <- spark_read_parquet(sc, ass_indirectby_ds_path)
  
  
  # Data about molecular interactions
  interactions <- spark_read_parquet(sc, interaction_path, memory = FALSE) %>%
    filter(sourceDatabase == "string") %>%
    filter(!is.na(targetA)) %>%
    filter(!is.na(targetB)) %>%
    filter(scoring > 0.42) %>%
    select(targetA, targetB) %>%
    sdf_distinct()

  #Filter by interactions and indirect ass.
  interactors_ass <- sig_gene_ids %>%
    inner_join(interactions, by = c("gene_id" = "targetA")) %>%
    inner_join(
      ass_indirectby_ds,
      by = c("targetB" = "targetId")
    ) %>%
    select(trait_gene_pairs,Traits,gene_id,targetB) %>%
    sdf_distinct() %>%
    collect()
  #write.csv(interactors_ass, "~/Desktop/interactions.csv")
  
  interactions_all <- c(interactions_all,nrow(interactors_ass))
  
  # New gene_trait pairs with interacting genes
  gene_interactions <- paste(interactors_ass$targetB,interactors_ass$Traits,sep="_")
  
  #compare to the true positive
  positive_control_set <- read_csv("/Users/dandantan/Desktop/IntAct_spark_dataset/positive_control_gene_list_with_opentarget.csv")
  positive_control_set$gene_trait_pairs <- paste(positive_control_set$ensg_gene_name, positive_control_set$Trait, sep = "_")
  
  # True Positives
  True_Positives <- length(intersect(gene_interactions,positive_control_set$gene_trait_pairs)) #23 20 23 20
  True_Positives1<- intersect(gene_interactions,positive_control_set$gene_trait_pairs)
  True_positives_all <- c(True_positives_all, True_Positives)
  
  # existing identified TP gene_trait pairs in the original results
  overlap_with_original <- length(intersect(True_Positives1, local_sig_gene_ids$trait_gene_pairs)) #9 7 9 8
  overlap_with_original1 <- intersect(True_Positives1, local_sig_gene_ids$trait_gene_pairs) 
  overlap_with_original_all<- c(overlap_with_original_all,overlap_with_original)
  
  # New identified gene_trait pairs with IntAct info
  new_identified <- length(setdiff(True_Positives1, overlap_with_original1)) #14 13 14 12
  new_identified_all <- c(new_identified_all, new_identified)
  
  # existing gene_trait pairs in original ExWAS results
  overlap_False_positive <- length(intersect(gene_interactions,all_sig$trait_gene_pairs)) # 24 34 27 28 
  overlap_False_positive_all <- c(overlap_False_positive_all,overlap_False_positive)
  
  # newly found interactions
  interactions_all2 <- c(interactions_all2,(nrow(interactors_ass)-overlap_False_positive))
}

#visualize the data
original_result <-c(196,341,317,407)
original_TP <- c(42,60,57,60)
false_positive <- original_result -original_TP

# interactions_all - overlap_False_positive
intact_result <- interactions_all2
intact_TP <- new_identified_all

withIntact_TP<- original_TP + intact_TP
withIntact_false_positive<- intact_result+original_result-withIntact_TP
  
x_axis<- c("pLoF","pLoF+alpha","pLoF+5in5miss","pLoF+5in5+1miss")

value<- c(original_TP,false_positive,withIntact_TP,withIntact_false_positive)
data <- data.frame(
  facet = rep(x_axis, each = 1), 
  group = rep(c("Original", "IntAct"), each = 8),
  stack = rep(c("True Positives","False Positives"), each = 4),
  value = value
)

ggplot(data, aes(x = group, y = value, fill = stack)) + 
  geom_bar(stat = "identity", position = "stack") +
  geom_text(aes(label = value), position = position_stack(vjust = 0.5), size = 3) +  # Add text labels
  facet_grid(~ facet) +
  labs(
    title = "Causal Gene Identification Across Four ExWAS Datasets With False Positives",
    x = "ExWAS Datasets",
    y = "Num of Significant Genes",
    fill = "Legend"
  ) +   scale_fill_manual(values = c("True Positives" = "lightblue", "False Positives" = "#FF9999")) +
  theme_minimal() + theme(plot.title = element_text(hjust = 0.25), text = element_text(size = 10))


## Only True Positives
value2<- c(original_TP,withIntact_TP)
data2 <- data.frame(
  facet = rep(x_axis, each = 1), 
  group = rep(c("Original", "IntAct"), each =4),
  stack = rep(c("True Positives"), times = 1),
  value = value2
)

ggplot(data2, aes(x = group, y = value, fill = stack)) + 
  geom_bar(stat = "identity", position = "stack") +
  geom_text(aes(label = value), position = position_stack(vjust = 0.5), size = 3) +  # Add text labels
  facet_grid(~ facet) +
  labs(
    title = "Causal Gene Identification Across Four ExWAS Datasets Without False Positives",
    x = "ExWAS Datasets",
    y = "Num of Causal Genes",
    fill = "Legend"
  ) + scale_fill_manual(values = "lightblue") +
  theme_minimal()+ theme(plot.title = element_text(hjust = 0.25), text = element_text(size = 10))
