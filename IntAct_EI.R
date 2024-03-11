library("tidyverse")
library("sparklyr")
library("sparklyr.nested")
library("cowplot")
library("ggsci")
library("ggplot2")
library("reshape2")
library("dplyr")
library("PRROC")
library("biomaRt")

###################################----------------------------Data Preprocesses--------------------------------#############################
# Read the EI algorithm prediction results
ei_path <-"/Users/dandantan/Desktop/IntAct_spark_dataset/ei_results.csv"
ei_results_all <- read.csv(ei_path)

# Read the positive control set (validate Ei result)
result_path<- "~/Desktop/IntAct_spark_dataset/positive_control_set_ei.csv"
positive_control_set_ei <- read.csv(result_path)
positive_control_set_ei$gene_trait_pairs <- paste(positive_control_set_ei$hgnc_gene_name, positive_control_set_ei$Trait, sep = "_")

# Cutoff with 0.46 to optimize performance (indicated in the paper)
ei_cutoff <- ei_results_all %>%
  filter(all.locus.prob >= 0.46)
ei_cutoff$gene_trait_pairs <- paste(ei_cutoff$names.genes, ei_cutoff$all.trait, sep = "_")
# Check if the EI-identified is in the positive control set
ei_cutoff$is_in_positive <- ifelse(ei_cutoff$gene_trait_pairs %in% positive_control_set_ei$gene_trait_pairs, 1, 0)

# Find all loci containing positive control genes for each trait
positive_loci <- ei_cutoff %>%
  group_by(all.trait) %>%
  filter(is_in_positive == 1) %>%
  dplyr::select(all.trait, locus.name)

# Only keep the loci with at least one positive control gene for each trait
ei_cutoff_loci <- ei_cutoff %>%
  semi_join(positive_loci, by = c("all.trait","locus.name"))

#######################--------------------------Add IntAct Information -------------------------#############################
# Find the highest EI score for each trait at each locus
highest_prob_per_locus <- ei_cutoff_loci %>%
  group_by(all.trait, locus.name) %>%
  filter(all.locus.prob == max(all.locus.prob, na.rm = TRUE)) #113 genes
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
ensembl_gene_id[[44,2]] = "ENSG00000171862"
ensembl_gene_id[[51,2]] = "ENSG00000103489"
ensembl_gene_id[[61,2]] = "ENSG00000214575"
ensembl_gene_id[[92,2]] = "ENSG00000275410"
ensembl_gene_id[[109,2]] = "ENSG00000084674"


highest_prob_per_locus_new <- highest_prob_per_locus %>% 
  left_join(ensembl_gene_id, by = c("names.genes" = "gene_name"))

highest_prob_per_locus_new$ensmble_genes <- unlist(highest_prob_per_locus_new$ensmble_genes)

### Find IntAct Partners
config <- spark_config()
sc <- spark_connect(master = "local", config = config)

local_path <- "/Users/dandantan/Desktop/Open_target_2021_FDA/Dataset"
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
positive_control_set_ei$gene_trait_pairs2<- paste(positive_control_set_ei$ensg_gene_name, positive_control_set_ei$Trait, sep = "_")

#check if the interacting genes in the positive control set
interactors_ass$is_in_positive <- ifelse(interactors_ass$gene_trait_pairs2 %in% positive_control_set_ei$gene_trait_pairs2, 1, 0)
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
sum(interactors_ass$is_in_positive) # 48/968

ei_cutoff_loci_IntAct <- bind_rows(interactors_ass, ei_cutoff_loci)%>%
  distinct()

####################-------------------------------- With Random selected genes---------------------------------################
ExWAS_results <- read.csv("/Users/dandantan/Desktop/IntAct_spark_dataset/ExWAS_results_all_5in5.csv")
difference <- nrow(ei_cutoff_loci_IntAct) - nrow(ei_cutoff_loci) # number of interacting genes
positive_control_set_ei2 <- positive_control_set_ei %>%
  mutate(
    Trait = ifelse(Trait == "ebmd", "ZBMD", Trait),
    Trait = ifelse(Trait == "tg", "IRNT_TG", Trait),
    Trait = ifelse(Trait == "t2d", "T2D", Trait),
    Trait = ifelse(Trait == "ldl", "IRNT_LDL", Trait),
    Trait = ifelse(Trait == "height", "IRNT_height", Trait),
    Trait = ifelse(Trait == "lowtsh", "lowtsh", Trait),
    Trait = ifelse(Trait == "rbc", "IRNT_RBC", Trait),
    Trait = ifelse(Trait == "dbp", "IRNT_DBP", Trait),
    Trait = ifelse(Trait == "calcium", "IRNT_Ca", Trait),
    Trait = ifelse(Trait == "sbp", "IRNT_SBP", Trait),
    Trait = ifelse(Trait == "glucose", "IRNT_glu", Trait),
    Trait = ifelse(Trait == "dbilirubin", "IRNT_biliru", Trait)
  )
positive_control_set_ei2$gene_trait_pairs2 <- paste(positive_control_set_ei2$ensg_gene_name, positive_control_set_ei2$Trait, sep = "_")
sums <- numeric(10000)
for (i in 1:10000) {
  random_indices <- sample(1:nrow(ExWAS_results), difference, replace = FALSE)
  random_genes <- ExWAS_results[random_indices, ]$trait_gene_pairs
  random_genes <- as.data.frame(random_genes)
  random_genes$is_in_positive <- ifelse(random_genes$random_genes %in% positive_control_set_ei2$gene_trait_pairs2, 1, 0)
  sums[i] <- sum(random_genes$is_in_positive)
}
average_sum <- mean(sums) # about 1

########################------------------Evaluate  with sensitivity and specificity and precision-------------------------##################################

### EI prediction
TP_ei<- sum(ei_cutoff_loci$is_in_positive) #117
FP_ei <- nrow(ei_cutoff_loci)-TP_ei #334

ei_results_all$gene_trait_pairs <- paste(ei_results_all$names.genes, ei_results_all$all.trait, sep = "_")
ei_results_all$is_in_positive <- ifelse(ei_results_all$gene_trait_pairs %in% positive_control_set_ei$gene_trait_pairs, 1, 0)

FN_ei<- sum(ei_results_all$is_in_positive) - TP_ei #63
TN_ei <- nrow(ei_results_all)-sum(ei_results_all$is_in_positive)-FP_ei #28411

sensitivity_ei<- TP_ei/(TP_ei+FN_ei) #0.65
specificity_ei<- TN_ei/(TN_ei+FP_ei) #0.9883
precision_ei<-TP_ei/(FP_ei+TP_ei) #0.2594

### EI prediction + IntAct
TP_ei_intact<- sum(ei_cutoff_loci_IntAct$is_in_positive) #165
FP_ei_intact <- nrow(ei_cutoff_loci_IntAct)-TP_ei_intact #1254

FN_ei_intact<- FN_ei #63
TN_ei_intact <- TN_ei #28411

sensitivity_ei_intact<- TP_ei_intact/(TP_ei_intact+FN_ei_intact) #0.72368
specificity_ei_intact<- TN_ei_intact/(TN_ei_intact+FP_ei_intact) #0.95772
precision_ei_intact<-TP_ei_intact/(FP_ei_intact+TP_ei_intact) #0.11627


## Visualize the data 
x_axis <- c("sensitivity","specificity","precision")
data_ei <- c(sensitivity_ei,specificity_ei,precision_ei)
data_ei_intact <- c(sensitivity_ei_intact,specificity_ei_intact,precision_ei_intact)

data <- data.frame(Method = rep(c("EI", "EI with IntAct"), each = 3),
                   Metric = rep(x_axis, 2),
                   Value = c(data_ei, data_ei_intact))

ggplot(data, aes(x = Metric, y = Value, fill = Method)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_text(aes(label = round(Value, 3)), vjust = -0.3, position = position_dodge(width = 0.9), size = 4) +
  labs(x = "Metric", y = "Value", title = "Performance Comparison of EI Methods") +
  theme_minimal() +
  theme(
    text = element_text(size = 12),  
    axis.text.x = element_text(size = 16),  
    axis.text.y = element_text(size = 12),
    axis.title.x = element_text(size = 14), 
    axis.title.y = element_text(size = 14), 
    plot.title = element_text(size = 18, hjust = 0.5),
    legend.text = element_text(size = 12),  
    legend.title = element_text(size = 14),
    legend.position = "bottom"
  )


###################################----------------------------Plot TP vs FP -----------------------#############################
total_TP <- c(TP_ei, TP_ei_intact)
total_FP <- c(FP_ei, FP_ei_intact)
groups =rep(c("True Positive","False Positive"), each =2)
methods <- rep(c("EI", "EI with IntAct"), each=1)

df2 <- data.frame(Method = methods, group =groups,value=c(total_TP,total_FP))


ggplot(df2, aes(x = Method, y = value, fill = group)) +
  geom_bar(stat = "identity", position = "stack") +
  geom_text(aes(label = round(value, 3)), vjust = 1.5, size = 4, color = "black",position = position_stack(vjust = 0.6)) +
  labs(x = "Method", y = "Num of Significant Genes", fill = "Group", title = "Performance of Causal Gene Identification with Different Methods") +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 12),  
        axis.text.y = element_text(size = 14),  
        axis.title.x = element_text(size = 14), 
        axis.title.y = element_text(size = 14),  
        plot.title = element_text(size = 18, hjust = 0.5),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        legend.position = "right")


