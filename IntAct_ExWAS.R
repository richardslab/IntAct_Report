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


# Initialize variables to store the results
True_positives_all <- c()
overlap_with_original_all <- c()
new_identified_all <- c()
overlap_False_positive_all <-c()
interactions_all <- c()
interactions_all2 <- c()

# Loop over the datasets
for (i in seq_along(all_sig_datasets)) {
  # all_sig and pathway datasets
  all_sig_pathway <- paste0("your file path", all_sig_datasets[[i]])
  TP <- paste0("your file path", all_sig_datasets[[i]])
  
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
    filter(sourceDatabase == "intact") %>%
    filter(!is.na(targetA)) %>%
    filter(!is.na(targetB)) %>%
    filter(scoring > 0.42) %>%
    filter(speciesA == speciesB)%>%  # only keep human species
    select(targetA, targetB) %>%
    sdf_distinct()
  
  
  #Filter by interactions and indirect ass.
  interactors_ass <- sig_gene_ids %>%
    inner_join(interactions, by = c("gene_id" = "targetA")) %>%
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

x_axis<- c("pLoF","pLoF with AlphaMissense","pLoF with Missense (5/5)","pLoF with Missense (1/5)")

value<- c(original_TP,false_positive,withIntact_TP,withIntact_false_positive)
data <- data.frame(
  facet = rep(x_axis, each = 1), 
  group = rep(c("Original", "IntAct"), each = 8),
  stack = rep(c("True Positives","False Positives"), each = 4),
  value = value
)

ggplot(data, aes(x = group, y = value, fill = stack)) + 
  geom_bar(stat = "identity", position = "stack") +
  geom_text(aes(label = value), position = position_stack(vjust = 0.6), size = 4.5) + 
  facet_grid(~ facet) +
  labs(
    title = "Causal Gene Identification Across Four ExWAS Datasets With False Positives",
    x = "ExWAS Datasets",
    y = "Num of Significant Genes",
    fill = "Legend"
  ) +   scale_fill_manual(values = c("True Positives" = "lightblue", "False Positives" = "#FF9999")) +
  theme_minimal() + 
  theme(strip.text = element_text(size = 14, color = "black")) +
  theme(plot.title = element_text(hjust = 0.5), text = element_text(size = 18))


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
  geom_text(aes(label = value), position = position_stack(vjust = 0.6), size = 5) + 
  facet_grid(~ facet) +
  labs(
    title = "Causal Gene Identification Across Four ExWAS Datasets Without False Positives",
    x = "ExWAS Datasets",
    y = "Num of Causal Genes",
    fill = "Legend"
  ) + scale_fill_manual(values = "lightblue") +
  theme_minimal()+   
  theme(strip.text = element_text(size = 14, color = "black")) +
  theme(plot.title = element_text(hjust = 0.5), text = element_text(size = 18))


############################--------------Evaluate not with roc/prc--------------------################################
total_positive <- c(682,698,694,698)
total_num_rows <- c(388442,393081,392014,393288)

TP_original<- c(42,60,57,60)
FP_original<- false_positive
FN_original<- total_positive-TP_original
TN_original<- total_num_rows-total_positive-FP_original

sensitivity_original<- TP_original/(TP_original+FN_original) # 0.06158358 0.08595989 0.08213256 0.08595989
specificity_original<- TN_original/(TN_original+FP_original) #0.9996028 0.9992839 0.9993356 0.9991161
precision_original<-TP_original/(FP_original+TP_original) # 0.2142857 0.1759531 0.1798107 0.1474201

## with IntAct
TP_IntAct<- withIntact_TP
FP_IntAct<- withIntact_false_positive
FN_IntAct<- FN_original
TN_IntAct<- TN_original

sensitivity_IntAct<- TP_IntAct/(TP_IntAct+FN_IntAct) # 0.08045977 0.10267229 0.10028249 0.10140845
specificity_IntAct<- TN_IntAct/(TN_IntAct+FP_IntAct) # 0.9983978 0.9972151 0.9972637 0.9969069
precision_IntAct<-TP_IntAct/(FP_IntAct+TP_IntAct) # 0.08259587 0.06250000 0.06206294 0.05585725

# Draw the plot
x_axis1 <- c("sensitivity", "specificity", "precision")
groups <- rep(c("Original ExWAS", "ExWAS+IntAct"), each = length(x_axis1) * 4)
categories <- rep(c("pLoF","pLoF with \n AlphaMissense","pLoF with \n Missense(5/5)","pLoF with \n Missense (1/5)"), times = 2)
data <- c(sensitivity_original[1],specificity_original[1],precision_original[1],
          sensitivity_original[2],specificity_original[2],precision_original[2],
          sensitivity_original[3],specificity_original[3],precision_original[3],
          sensitivity_original[4],specificity_original[4],precision_original[4],
          sensitivity_IntAct[1], specificity_IntAct[1], precision_IntAct[1],
          sensitivity_IntAct[2], specificity_IntAct[2], precision_IntAct[2],
          sensitivity_IntAct[3], specificity_IntAct[3], precision_IntAct[3],
          sensitivity_IntAct[4], specificity_IntAct[4], precision_IntAct[4]
)

# Create a data frame
df <- data.frame(
  x_axis = rep(categories, each = length(x_axis1)),
  groups = groups,
  categories = rep(x_axis1, times = 8),
  data = data
)

# Plotting
ggplot(df, aes(x = x_axis, y = data, fill = groups)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.85), width = 0.8) + 
  geom_text(aes(label = round(data, 3)), position = position_dodge(width = 0.85), vjust = -0.5, size = 4.5) +
  
  facet_wrap(~ categories, scales = "free_x") +
  labs(
    title = "Comparison of Metrics between Original ExWAS and ExWAS+IntAct",
    x = "Datasets",
    y = "Value",
    fill = "Groups"
  ) +
  theme_minimal()+
  theme(strip.text = element_text(size = 14, color = "black")) +
  theme(strip.text.x = element_text(margin = margin(b = 20)), 
        legend.position = "bottom",  
        legend.spacing.x = unit(0.2, "cm"),
        plot.title = element_text(hjust = 0.5,size =18),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 15))

