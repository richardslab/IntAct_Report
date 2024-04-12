library("tidyverse")
library("sparklyr")
library("sparklyr.nested")
library("ggsci")
library("ggplot2")
library("reshape2")
library("dplyr")
library("biomaRt")

###################################----------------------------Data Preprocesses--------------------------------#############################
# Read the EI algorithm prediction results
ei_path <-"YourPathways/ei_results.csv"
ei_results_all <- read.csv(ei_path)

# Read the positive control set (validate Ei result)
result_path<- "YourPathways/ei_positive_control_gene_list.csv"
positive_control_set_ei <- read.csv(result_path)
positive_control_set_ei$gene_trait_pairs <- paste(positive_control_set_ei$Gene.Symbol, positive_control_set_ei$Trait, sep = "_")


ei_results_all$gene_trait_pairs <- paste(ei_results_all$names.genes, ei_results_all$all.trait, sep = "_")

# Check if the EI-identified is in the positive control set
ei_results_all$is_in_positive <- ifelse(ei_results_all$gene_trait_pairs %in% positive_control_set_ei$gene_trait_pairs, 1, 0)

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
  filter(all.locus.prob == max(all.locus.prob, na.rm = TRUE)) #150 genes

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
ensembl_gene_id[[22,2]] = "ENSG00000225663"
ensembl_gene_id[[41,2]] = "ENSG00000175164"
ensembl_gene_id[[71,2]] = "ENSG00000171862"
ensembl_gene_id[[81,2]] = "ENSG00000103489"
ensembl_gene_id[[95,2]] = " ENSG00000186073"
ensembl_gene_id[[98,2]] = "ENSG00000214575"
ensembl_gene_id[[99,2]] = "ENSG00000007968"
ensembl_gene_id[[135,2]] = "ENSG00000275410"
ensembl_gene_id[[140,2]] = "ENSG00000053918"
ensembl_gene_id[[146,2]] = "ENSG00000084674"

highest_prob_per_locus_new <- highest_prob_per_locus %>% 
  left_join(ensembl_gene_id, by = c("names.genes" = "gene_name"))

highest_prob_per_locus_new$ensmble_genes <- unlist(highest_prob_per_locus_new$ensmble_genes)

#######################--------------------------Add IntAct Information -------------------------#############################

### Find IntAct Partners
config <- spark_config()
sc <- spark_connect(master = "local", config = config)

local_path <- "YourPathways/Open_target_2021_FDA/Dataset"
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
positive_control_set_ei$gene_trait_pairs2<- paste(positive_control_set_ei$ensemble_gene, positive_control_set_ei$Trait, sep = "_")

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
sum(interactors_ass$is_in_positive) # 47/1417

ei_loci_IntAct <- bind_rows(interactors_ass, highest_prob_per_locus_new)%>%
  distinct()

###########------------------ With Random selected genes--------------##########
ExWAS_results <- read.csv("YourPathways/ExWAS_results_all_5in5.csv")
difference <- nrow(ei_loci_IntAct) - nrow(highest_prob_per_locus_new) # number of interacting genes
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
TP_random <- 59 + average_sum
FP_random <- 91 + difference -average_sum
FN_random <- 100
TN_random <- 28675
precision_random <- TP_random/(TP_random+FP_random) 
sensitivity_random <- TP_random/(TP_random+FN_random)
specificity_random<- TN_random/(TN_random+FP_random)

########################------------------Evaluate  with sensitivity and specificity and precision-------------------------##################################

### EI prediction
TP_ei<- sum(highest_prob_per_locus_new$is_in_positive) #59
FP_ei <- nrow(highest_prob_per_locus_new)-TP_ei #91

ei_results_all$gene_trait_pairs <- paste(ei_results_all$names.genes, ei_results_all$all.trait, sep = "_")
ei_results_all$is_in_positive <- ifelse(ei_results_all$gene_trait_pairs %in% positive_control_set_ei$gene_trait_pairs, 1, 0)

FN_ei<- sum(ei_results_all$is_in_positive) - TP_ei #100
TN_ei <- nrow(ei_results_all)-sum(ei_results_all$is_in_positive)-FP_ei #28675

sensitivity_ei<- TP_ei/(TP_ei+FN_ei) 
specificity_ei<- TN_ei/(TN_ei+FP_ei) 
precision_ei<-TP_ei/(FP_ei+TP_ei) 

### EI prediction + IntAct
TP_ei_intact<- sum(ei_loci_IntAct$is_in_positive) 
FP_ei_intact <- nrow(ei_loci_IntAct)-TP_ei_intact 

FN_ei_intact<- FN_ei 
TN_ei_intact <- TN_ei 

sensitivity_ei_intact<- TP_ei_intact/(TP_ei_intact+FN_ei_intact) 
specificity_ei_intact<- TN_ei_intact/(TN_ei_intact+FP_ei_intact) 
precision_ei_intact<-TP_ei_intact/(FP_ei_intact+TP_ei_intact) 


## Visualize the data 
x_axis <- c("sensitivity","specificity","precision")
data_ei <- c(sensitivity_ei,specificity_ei,precision_ei)
data_ei_intact <- c(sensitivity_ei_intact,specificity_ei_intact,precision_ei_intact)

data <- data.frame(Method = rep(c("EI", "EI with IntAct"), each = 3),
                   Metric = rep(x_axis, 2),
                   Value = c(data_ei, data_ei_intact))


ggplot(data, aes(x = Metric, y = Value, fill = Method)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_text(aes(label = round(Value, 3)), vjust = -0.2, position = position_dodge(width = 0.9), size = 9) +
  labs(x = "Metric", y = "Value") +
  theme_minimal() +
  theme(
    text = element_text(size = 12),  
    axis.text.x = element_text(size = 22),  
    axis.text.y = element_text(size = 20),
    axis.title.x = element_text(size = 20), 
    axis.title.y = element_text(size = 20), 
    legend.text = element_text(size = 20),  
    legend.title = element_text(size = 20),
    legend.position = "bottom"
  ) +
  scale_fill_manual(values = c("#d6604d","#4393c3"))

###################################----------------------------Plot TP vs FP -----------------------#############################
total_TP <- c(TP_ei, TP_ei_intact)
total_FP <- c(FP_ei, FP_ei_intact)
groups =rep(c("True Positive","False Positive"), each =2)
methods <- rep(c("EI", "EI with IntAct"), each=1)

df2 <- data.frame(Method = methods, group =groups,value=c(total_TP,total_FP))


ggplot(df2, aes(x = Method, y = value, fill = group)) +
  geom_bar(stat = "identity", position = "stack") +
  #geom_text(aes(label = round(value, 3)), vjust = 1, size = 9, color = "black",position = position_stack(vjust = 0.6)) +
  labs(x = "Method", y = "Num of Significant Genes", fill = "Legend",) +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 22),  
        axis.text.y = element_text(size = 20),  
        axis.title.x = element_text(size = 20), 
        axis.title.y = element_text(size = 20),  
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20),
        legend.position = "bottom")+
  scale_fill_manual(values = c("#9970ab","#5aae61"))


