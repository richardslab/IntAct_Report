library("tidyverse")
library("sparklyr")
library("sparklyr.nested")
library("ggsci")
library("ggplot2")


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
  all_sig_pathway <- paste0("ExWAS_Data/", all_sig_datasets[[i]])
  TP <- paste0("ExWAS_Data/", all_sig_datasets[[i]])
  
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
  local_path <- "YourPathway/Open_target_2021_FDA/Dataset" # data released 21.11
  
  interaction_path <- paste(
    local_path,
    "/interaction/",
    sep = ""
  )
  
  # Data about molecular interactions
  interactions <- spark_read_parquet(sc, interaction_path, memory = FALSE) %>%
    filter(sourceDatabase == "intact") %>%
    filter(!is.na(targetA)) %>%
    filter(!is.na(targetB)) %>%
    filter(scoring > 0.42) %>%
    filter(speciesA == speciesB)%>%  # only keep human species
    select(targetA, targetB) %>%
    sdf_distinct()
  
  #Filter by interactions
  interactors_ass <- sig_gene_ids %>%
    inner_join(interactions, by = c("gene_id" = "targetA")) %>%
    select(trait_gene_pairs,Traits,gene_id,targetB) %>%
    sdf_distinct() %>%
    collect()
  
  interactions_all <- c(interactions_all,nrow(interactors_ass))
  
  # New gene_trait pairs with interacting genes
  gene_interactions <- paste(interactors_ass$targetB,interactors_ass$Traits,sep="_")
  
  #compare to the true positive
  positive_control_set <- read_csv("ExWAS_Data/positive_control_gene_list.csv")
  #change trait names to match the data
  positive_control_set <- positive_control_set %>%
    mutate(
      Trait = ifelse(Trait == "ebmd", "ZBMD", Trait),
      Trait = ifelse(Trait == "tg", "IRNT_TG", Trait),
      Trait = ifelse(Trait == "t2d", "T2D", Trait),
      Trait = ifelse(Trait == "ldl", "IRNT_LDL", Trait),
      Trait = ifelse(Trait == "height", "IRNT_height", Trait),
      Trait = ifelse(Trait == "BMI", "IRNT_BMI", Trait),
      Trait = ifelse(Trait == "lowtsh", "lowtsh", Trait),
      Trait = ifelse(Trait == "rbc", "IRNT_RBC", Trait),
      Trait = ifelse(Trait == "dbp", "IRNT_DBP", Trait),
      Trait = ifelse(Trait == "calcium", "IRNT_Ca", Trait),
      Trait = ifelse(Trait == "WHR", "IRNT_WHR", Trait),
      Trait = ifelse(Trait == "sbp", "IRNT_SBP", Trait),
      Trait = ifelse(Trait == "glucose", "IRNT_glu", Trait),
      Trait = ifelse(Trait == "dbilirubin", "IRNT_biliru", Trait)
    )
  positive_control_set$gene_trait_pairs <- paste(positive_control_set$ensg_gene_name, positive_control_set$Trait, sep = "_")
  
  # True Positives
  True_Positives <- length(intersect(gene_interactions,positive_control_set$gene_trait_pairs)) 
  True_Positives1<- intersect(gene_interactions,positive_control_set$gene_trait_pairs)
  True_positives_all <- c(True_positives_all, True_Positives)
  
  # existing identified TP gene_trait pairs in the original results
  overlap_with_original <- length(intersect(True_Positives1, local_sig_gene_ids$trait_gene_pairs)) 
  overlap_with_original1 <- intersect(True_Positives1, local_sig_gene_ids$trait_gene_pairs) 
  overlap_with_original_all<- c(overlap_with_original_all,overlap_with_original)
  
  # New identified gene_trait pairs with IntAct info
  new_identified <- length(setdiff(True_Positives1, overlap_with_original1)) 
  new_identified_all <- c(new_identified_all, new_identified)
  
  # existing gene_trait pairs in original ExWAS results
  overlap_False_positive <- length(intersect(gene_interactions,all_sig$trait_gene_pairs)) 
  overlap_False_positive_all <- c(overlap_False_positive_all,overlap_False_positive)
  
  # newly found interactions
  interactions_all2 <- c(interactions_all2,(nrow(interactors_ass)-overlap_False_positive))
}

#####################---------------------------------Try with Randomly selected genes--------------------------------#########################
differences<- interactions_all2
ExWAS_results <- read.csv("/Users/dandantan/Desktop/IntAct_spark_dataset/ExWAS_results_all_5in5.csv")

total_average<-c()

sums <- numeric(10000)
for (i in 1:10000) {
    difference <- differences[x]
    random_indices <- sample(1:nrow(ExWAS_results), difference, replace = FALSE)
    random_genes <- ExWAS_results4[random_indices, ]$trait_gene_pairs
    random_genes <- as.data.frame(random_genes)
    random_genes$is_in_positive <- ifelse(random_genes$random_genes %in% positive_control_set$gene_trait_pairs, 1, 0)
    sums[i] <- sum(random_genes$is_in_positive)
  }
  average_sum <- mean(sums)

cat("Averaged increased of True Positive with randomly selected pLoF with Missense (5/5) dataset:",average_sum)
######################-------------------------------visualize the data-----------------------------------------------################################
original_result <-c(196,341,317,407)
TP_original <- c(42,60,57,60)
false_positive <- original_result -original_TP

intact_result <- interactions_all2
intact_TP <- new_identified_all

withIntact_TP<- TP_original + intact_TP
withIntact_false_positive<- intact_result+original_result-withIntact_TP

x_axis<- c("pLoF","pLoF with \nAlphaMissense","pLoF with \nMissense (5/5)","pLoF with\n Missense (1/5)")

value<- c(TP_original,false_positive,withIntact_TP,withIntact_false_positive)
data <- data.frame(
  facet = rep(x_axis, each = 1), 
  group = rep(c("Original \nExWAS", "ExWAS\nwith\nIntAct"), each = 8),
  stack = rep(c("True Positives","False Positives"), each = 4),
  value = value
)

ggplot(data, aes(x = group, y = value, fill = stack)) + 
  geom_bar(stat = "identity", position = "stack") +
  geom_text(aes(label = value), position = position_stack(vjust = 0.7), size = 6.5) + 
  facet_grid(~ facet) +
  labs(
    x = "ExWAS Datasets",
    y = "Num of Significant Genes",
    fill = "Legend"
  ) +
  theme_minimal() + 
  theme(strip.text = element_text(size = 20, color = "black"),
        axis.text.x = element_text(size = 16),
        axis.title.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        legend.text=element_text(size=20))

## Only True Positives
value2<- c(TP_original,withIntact_TP)
data2 <- data.frame(
  facet = rep(x_axis, each = 1), 
  group = rep(c("Original \nExWAS", "ExWAS\nwith\nIntAct"), each =4),
  stack = rep(c("True Positives"), times = 1),
  value = value2
)

ggplot(data2, aes(x = group, y = value, fill = stack)) + 
  geom_bar(stat = "identity", position = "stack") +
  geom_text(aes(label = value), position = position_stack(vjust = 0.6), size = 9) + 
  facet_grid(~ facet) +
  labs(
    x = "ExWAS Datasets",
    y = "Num of Causal Genes",
    fill = "Legend"
  ) + scale_fill_manual(values = "lightblue") +
  theme_minimal()+   
  theme(strip.text = element_text(size = 20, color = "black"),
        axis.text.x = element_text(size = 16),
        axis.title.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        legend.text=element_text(size=20))


############################--------------Evaluate with sensitivity, precison and specificity--------------------################################
total_positive <- c(682,698,694,698)
total_num_rows <- c(388442,393081,392014,393288)

TP_original<- c(42,60,57,60)
FP_original<- false_positive
FN_original<- total_positive-TP_original
TN_original<- total_num_rows-total_positive-FP_original

sensitivity_original<- TP_original/(TP_original+FN_original) 
specificity_original<- TN_original/(TN_original+FP_original) 
precision_original<-TP_original/(FP_original+TP_original)

## with IntAct
TP_IntAct<- withIntact_TP
FP_IntAct<- withIntact_false_positive
FN_IntAct<- FN_original
TN_IntAct<- TN_original

sensitivity_IntAct<- TP_IntAct/(TP_IntAct+FN_IntAct) 
specificity_IntAct<- TN_IntAct/(TN_IntAct+FP_IntAct) 
precision_IntAct<-TP_IntAct/(FP_IntAct+TP_IntAct) 

# Draw the plot
precision_all <- c(precision_original[1],precision_original[2],precision_original[3],precision_original[4],
                   precision_IntAct[1],precision_IntAct[2],precision_IntAct[3],precision_IntAct[4])
sensitivity_all<- c(sensitivity_original[1],sensitivity_original[2],sensitivity_original[3],sensitivity_original[4],
                    sensitivity_IntAct[1],sensitivity_IntAct[2],sensitivity_IntAct[3],sensitivity_IntAct[4])
specificity_all <- c(specificity_original[1],specificity_original[2],specificity_original[3],specificity_original[4],
                    specificity_IntAct[1],specificity_IntAct[2],specificity_IntAct[3],specificity_IntAct[4])

categories <- rep(c("pLoF","pLoF with \n AlphaMissense","pLoF with \n Missense(5/5)","pLoF with \n Missense (1/5)"))

plot_precision <- data.frame(
  facet = rep(categories, each = 1), 
  group = rep(c("Original \nExWAS", "ExWAS\nwith\nIntAct"), each =4),
  value = precision_all
)

plot_sensitivity <- data.frame(
  facet = rep(categories, each = 1), 
  group = rep(c("Original \nExWAS", "ExWAS\nwith\nIntAct"), each =4),
  value = sensitivity_all
)

plot_specificity <- data.frame(
  facet = rep(categories, each = 1), 
  group = rep(c("Original \nExWAS", "ExWAS\nwith\nIntAct"), each =4),
  value = specificity_all
)


##Precision
ggplot(plot_precision, aes(x = group, y = value, fill = group)) + 
  geom_bar(stat = "identity", position = "stack") +
  geom_text(aes(label = round(value, 2)), position = position_stack(vjust = 0.6), size = 9) + 
  facet_grid(~ facet) +
  labs(
    x = "ExWAS Datasets",
    y = "Precision",
    fill = "Legend"
  ) +
  theme_minimal()+   
  theme(strip.text = element_text(size = 20, color = "black"),
        axis.text.x = element_text(size = 16),
        axis.title.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        legend.text=element_text(size=16))

##Sensitivity
ggplot(plot_sensitivity, aes(x = group, y = value, fill = group)) + 
  geom_bar(stat = "identity", position = "stack") +
  geom_text(aes(label = round(value, 2)), position = position_stack(vjust = 0.6), size = 9) + 
  facet_grid(~ facet) +
  labs(
    x = "ExWAS Datasets",
    y = "Sensitivity",
    fill = "Legend"
  ) +
  theme_minimal()+   
  theme(strip.text = element_text(size = 20, color = "black"),
        axis.text.x = element_text(size = 16),
        axis.title.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        legend.text=element_text(size=16))


##Specificity
ggplot(plot_specificity, aes(x = group, y = value, fill = group)) + 
  geom_bar(stat = "identity", position = "stack") +
  geom_text(aes(label = round(value, 2)), position = position_stack(vjust = 0.6), size = 9) + 
  facet_grid(~ facet) +
  labs(
    x = "ExWAS Datasets",
    y = "Specificity",
    fill = "Legend"
  ) +
  theme_minimal()+   
  theme(strip.text = element_text(size = 20, color = "black"),
        axis.text.x = element_text(size = 16),
        axis.title.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        legend.text=element_text(size=16))



