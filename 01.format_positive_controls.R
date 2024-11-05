library("dplyr")

###In the study we have three methods: ExWAS, Effector Index, and Genetic Priority Score
#Load dataset used to evaluate ExWAS
ref_genes_exwas<-read.csv("/Users/dandantan/Desktop/IntAct_spark_dataset/positive_control_gene_list_with_opentarget.csv")
#Load dataset used to evaluate EI
ref_genes_EI<-read.csv("/Users/dandantan/Desktop/IntAct_spark_dataset/positive_control_set_ei.csv")
#Load dataset used to evaluate GPS
ref_genes_GPS<-read.csv("/Users/dandantan/Desktop/IntAct_spark_dataset/GPS_positive_results.csv")

#create gene_trait pair
ref_genes_exwas$gene_trait_pair <- paste(ref_genes_exwas$ensg_gene_name, ref_genes_exwas$Trait, sep = "_")
ref_genes_EI$gene_trait_pair <- paste(ref_genes_EI$ensg_gene_name, ref_genes_EI$Trait, sep = "_")
ref_genes_GPS$gene_trait_pair <- paste(ref_genes_GPS$ensmble_genes, ref_genes_GPS$Phecode.description, sep = "_")



###merge exwas and ei
merge_exwas_ei <- ref_genes_exwas %>%
  left_join(ref_genes_EI, by = "gene_trait_pair") %>%
  distinct(gene_trait_pair, .keep_all = TRUE) %>%  # Keep unique rows based on 'gene_trait_pair'
  select(hgnc_gene_name.x, ensg_gene_name.x,Trait.x, gene_trait_pair, source.x)  # Select the desired columns
#write.csv(merge_exwas_ei,file = "/Users/dandantan/Desktop/IntAct_spark_dataset/merged_exwas_ei.csv", row.names = FALSE)


###Change traits name to match with GPS
merge_exwas_ei <- merge_exwas_ei %>%
  mutate(
    Trait.x = ifelse(Trait.x == "hypertension", "Hypertension", Trait.x),
    Trait.x = ifelse(Trait.x == "Hypercholesterolemia", "Disorders of lipoid metabolism", Trait.x),
    Trait.x = ifelse(Trait.x == "Cataract", "Cataract", Trait.x),
    Trait.x = ifelse(Trait.x == "t2d", "Diabetes mellitus", Trait.x),
    Trait.x = ifelse(Trait.x == "Major_depressive_disorder", "Mood disorders", Trait.x),
    Trait.x = ifelse(Trait.x == "Atrial_fibrillation", "Cardiac dysrhythmias", Trait.x),
    Trait.x = ifelse(Trait.x == "Cancer_of_prostate", "Cancer of prostate", Trait.x),
    Trait.x = ifelse(Trait.x == "Breast_cancer", "Breast cancer", Trait.x),
    Trait.x = ifelse(Trait.x == "Osteoarthritis_localized", "Osteoarthrosis", Trait.x)
  )
merge_exwas_ei$gene_trait_pair<-paste(merge_exwas_ei$ensg_gene_name.x, merge_exwas_ei$Trait.x, sep = "_")

###merge exwas,ei with gps (just directly add the gene_traits_pair)
# Extract columns from merge_exwas_ei
exwas_data <- data.frame(
  gene = merge_exwas_ei$hgnc_gene_name.x,
  ensmble_gene_gps = merge_exwas_ei$ensg_gene_name.x,
  trait_exwas_ei=merge_exwas_ei$Trait.x,
  trait_gps=NA,
  gene_trait_pair_exwas_ei = merge_exwas_ei$gene_trait_pair,
  gene_trait_pair_gps = NA,  # Placeholder for missing GPS data
  gene_tait_pair_merged=merge_exwas_ei$gene_trait_pair,
  source="ExWAS, EI",
  stringsAsFactors = FALSE
)

# Extract columns from ref_genes_GPS
gps_data <- data.frame(
  gene = ref_genes_GPS$Gene,
  ensmble_gene_gps = ref_genes_GPS$ensmble_genes,
  trait_exwas_ei=NA,
  trait_gps=ref_genes_GPS$Phecode.description,
  gene_trait_pair_exwas_ei = NA,  # Placeholder for missing EXWAS/ei data
  gene_trait_pair_gps = ref_genes_GPS$gene_trait_pair,
  gene_tait_pair_merged=ref_genes_GPS$gene_trait_pair,
  source="GPS",
  stringsAsFactors = FALSE
)

# Append the two datasets vertically
merged_positive_gene_list <- rbind(exwas_data, gps_data)
