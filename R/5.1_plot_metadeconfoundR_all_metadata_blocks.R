library(devtools)
library(tidyr)
library(tibble)
library(purrr)
require(phyloseq)
require(tidyverse)
require(magrittr)

ps_0 <- readRDS("/fast/AG_Forslund/Lola/CIPF_2018/Analysis/Rats_HE_CCL4model/output/ps_0.Rds")
library(readxl)
metadata <- read_excel("/fast/AG_Forslund/Lola/CIPF_2018/Analysis/Rats_HE_CCL4model/metadata.xlsx")
View(metadata)


get_latest_annotation <- function(phyloseq_obj) {
  tax_table <- phyloseq_obj@tax_table %>%
    as.data.frame() %>%
    rownames_to_column('ASV') %>%
    as_tibble() %>%
    tidyr::gather('Rank', 'Value', -ASV) %>%
    nest(data = -ASV) %>%
    mutate(TaxaID = map_chr(data, function(x) {
      x %>%
        dplyr::filter(!is.na(Value)) %>%
        magrittr::use_series(Value) %>%
        tail(1)
    })) %>%
    mutate(TaxaUp = map_chr(data, function(x) {
      x %>%
        dplyr::filter(!is.na(Value)) %>%
        magrittr::use_series(Value) %>%
        tail(2) %>% head(1)
    })) %>%
    mutate(TaxaID = paste(ASV, TaxaID)) %>%
    dplyr::select(-c(data, TaxaUp))
  
  return(tax_table)
}

#Heatmap including all variables
##MetadeconfoundR
library(metadeconfoundR)
library(gtools)
library(dplyr)
ps.rel<-microbiome::transform(ps_0, "compositional")
latest_annotations <- get_latest_annotation(ps.rel)
latest_annotations<-as.data.frame(latest_annotations)
short_names<-as.data.frame(otu_table(ps.rel))
rownames(short_names)<-latest_annotations$TaxaID
Meta_Species <- t(short_names) 

metadata_R <- subset(metadata, select = -c(SampleID, Week, FisabioID, GABAa1, GABAa5, Ymaze_trialstolearn, OLM, IL1_b, Group))
metadata_R <- metadata_R %>% mutate(Batch = ifelse(Batch=="T5",1,0))
metadata_R <- metadata_R %>% mutate(Antibiotic = ifelse(Antibiotic=="Yes",1,0))
metadata_R <- metadata_R %>% mutate(Treated_CCL = ifelse(Treated_CCL=="Yes",1,0))
metadata_R <- as.data.frame(metadata_R)
rownames(metadata_R) <- (metadata$SampleID)

# check correct ordering
all(rownames(metadata_R) == rownames(Meta_Species))
## [1] TRUE
all(order(rownames(metadata_R)) == order(rownames(Meta_Species)))
## [1] TRUE
Meta_Species<- as.data.frame(Meta_Species)

Output_batch <- MetaDeconfound(featureMat = as.data.frame(Meta_Species),
                               metaMat = as.data.frame(metadata_R), nnodes = 14, randomVar = list("+ (1|Batch)",
                                                                                                  c("Batch")))
View(Output_batch)
saveRDS(Output_batch, "~/Rats_HE_CCL4model/output/MetadeconfoundR/Output_batch_figure")

right <- BuildHeatmap(Output_batch_figure)
Metadec_table<-data.frame(right[["data"]])
Metadec_table<-Metadec_table%>%
  mutate(Blocks = case_when(
    metaVariable %in% c("Ambulatory.Counts", "Vertcal.Counts", "Stereotipic.Counts", "Learning_index_Rmaze", "Average.velocity") ~ "Cognition",
    metaVariable %in% c("Antibiotic", "Treated_CCL") ~ "Treatments",
    metaVariable %in% c("Ammonia") ~ "Ammonia levels",
    metaVariable %in% c("CCL2", "CCL5","CCR5","CCL20","CCR2", "CX3CL1", "CX3CR1", "IFN.gamma", "IL_10", "IL.15", "IL.17", "IL_4", "IL.6", "Occludine",  "TNF_a", "TGF.b") ~ "Cytokines",
    metaVariable %in% c("GLAST", "GLT1", "GAT1", "GAT3", "TNFR1", "NR2B", "NR2A", "NR1", "GLUA1", "GLUA2") ~ "Memb. receptors",
    metaVariable %in% c("SCFA_AA", "SCFA_BA", "SCFA_CA", "SCFA_PA", "SCFA_VA") ~ "SCFAs"))%>%
  mutate(ordernames1 = case_when(
    metaVariable %in% c("Antibiotic") ~ Ds))
# Define the desired order of levels for facet_wrap
facet_order <- c("Treatments", "SCFAs", "Memb. receptors", "Cytokines", "Cognition")
# Convert the Blocks column to a factor with the desired order
Metadec_table$Blocks <- factor(Metadec_table$Blocks, levels = facet_order)
Metadec_table$status_ok <- Metadec_table$status
# Replace "_" with "."
Metadec_table$metaVariable <- gsub("_", ".", Metadec_table$metaVariable)
Metadec_table <- Metadec_table %>%
  select(stars, status_ok, feature, metaVariable)


View(Output_batch_figure)
raw_p <- Output_batch_figure[1]
corr_p <- Output_batch_figure[2]
effect_size<- Output_batch_figure[3]
status<- Output_batch_figure[4]

raw_p_df <- data.frame(raw_p$Ps)
raw_p_df  <- raw_p_df  %>% 
  select(-c(Batch))%>%
  rename_with(~ gsub("_", ".", .))%>%
  rename_all(~paste0("p_",.))%>%
  rownames_to_column()


corr_p_df <- data.frame(corr_p$Qs)
corr_p_df  <- corr_p_df %>%  
  select(-c(Batch))%>%
  rename_with(~ gsub("_", ".", .))%>%
  rename_all(~paste0("q_",.))%>%
  rownames_to_column()


effect_size_df <- data.frame(effect_size$Ds)
effect_size_df  <- effect_size_df %>%  
  select(-c(Batch))%>%
  rename_with(~ gsub("_", ".", .))%>%
  rename_all(~paste0("effect_size_",.))%>%
  rownames_to_column()


status_df <- data.frame(status$status)
status_df  <- status_df %>%  
  select(-c(Batch))%>%
  rename_with(~ gsub("_", ".", .))%>%
  rename_all(~paste0("status_",.))%>%
  rownames_to_column()

#create two-column-dataframe containing corresponding "human-readable" names to the "machine-readable" feature names used as row.names in metaDeconfOutput.  
taxtable <- latest_annotations
taxtable$rowname <- latest_annotations$ASV
taxtable<- cbind(rowname=taxtable$rowname,subset(taxtable,select = -c(rowname)))

effect_table <- raw_p_df%>%
  full_join(corr_p_df, by="rowname")%>%
  full_join(effect_size_df, by="rowname")%>%
  full_join(status_df, by="rowname")%>%
  full_join(taxtable, by="rowname")

# remove the entries which have NS and AD in status
effect_table_sig <- effect_table%>%
  filter(if_any(starts_with("status_"), ~. != "NS"))


#pivot long format

effect_table_sig_long <- effect_table_sig%>%
  pivot_longer(cols = starts_with("status"), names_to = "comparison_status", values_to = "status")%>%
  separate(comparison_status, c("name" ,"variable"), "_")%>%
  mutate(comparison_status=paste(variable, sep="_"))%>%
  select(-c(name, variable))%>%
  pivot_longer(cols = starts_with("p_"), names_to = "comparison_p", values_to = "raw_p")%>%
  separate(comparison_p, c("name","variable"), "_")%>%
  mutate(comparison_p=paste(variable, sep="_"))%>%
  select(-c(name, variable))%>%
  filter(comparison_p==comparison_status)%>%
  pivot_longer(cols = starts_with("effect_size"), names_to = "comparison_effectSize", values_to = "effectSize")%>%
  separate(comparison_effectSize, c("name","size", "variable"), "_")%>%
  mutate(comparison_effectSize=paste(variable, sep="_"))%>%
  select(-c(name, size, variable))%>%
  filter(comparison_p==comparison_effectSize)%>%
  pivot_longer(cols = starts_with("q"), names_to = "comparison_q", values_to = "corr_p")%>%
  separate(comparison_q, c("name","variable"), "_")%>%
  mutate(comparison_q=paste(variable, sep="_"))%>%
  select(-c(name))%>%
  filter(comparison_p==comparison_q)%>%
  select(-c(comparison_q, comparison_effectSize,comparison_status, ASV, TaxaID))%>%
  mutate(Blocks = case_when(
    comparison_p %in% c("Ambulatory.Counts", "Vertcal.Counts", "Stereotipic.Counts", "Learning.index.Rmaze", "Average.velocity") ~ "Cognitive tests",
    comparison_p %in% c("Antibiotic", "Treated.CCL") ~ "Treatments",
    comparison_p %in% c("Ammonia") ~ "Ammonia levels",
    comparison_p %in% c("CCL2", "CCL5","CCR5","CCL20","CCR2", "CX3CL1", "CX3CR1", "IFN.gamma", "IL.10", "IL.15", "IL.17", "IL.4", "IL.6", "Occludine","IL_4", "TNF.a", "TGF.b") ~ "Cytokines brain",
    comparison_p %in% c("GLAST", "GLT1", "GAT1", "GAT3", "TNFR1", "NR2B", "NR2A", "NR1", "GLUA1", "GLUA2") ~ "Membrane expression receptors",
    comparison_p %in% c("SCFA.AA", "SCFA.BA", "SCFA.CA", "SCFA.PA", "SCFA.VA") ~ "SCFAs")) %>%
  left_join(Metadec_table, by = c("rowname" = "feature", "variable" = "metaVariable")) %>%
  mutate(stars = ifelse(is.na(stars), "", stars))%>%
  mutate(ordernames1 = case_when(
    comparison_p %in% c("Antibiotic") ~ effectSize))


#Heatmap all variables
heatmap_all<-ggplot(effect_table_sig_long, aes(x = variable, y = fct_reorder(rowname, ordernames1, .fun = function(x) mean(x, na.rm = TRUE)))) + geom_tile(aes(fill = effectSize)) +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0) +
  geom_text(aes(label = stars), size = 2, key_glyph = "point") +
  theme_classic()  +
  theme(
    axis.text.x = element_text(size = 7, angle = 90, hjust = 1, vjust = 0.3),
    axis.text.y = element_text(size = 7, angle = 0, hjust = 1, vjust = 0.35),
    plot.title.position = "plot",
    plot.title = element_text(hjust = 0),
    plot.subtitle = element_text(size = 8) ,
  ) +
  labs(
    title = "Significant metadata correlation heatmap",
    subtitle = "FDR-values: < 0.001 = ***, < 0.01 = **, < 0.1 = *",
    x = "Variables",
    y = "ASVs agglomerated by Genus",
    fill = "Effect size", 
    shape = "Confounding status"
  )+
  facet_grid(. ~ (effect_table_sig_long$Blocks), scales = "free", space='free', switch = "x") 

heatmap_all
ggsave("/fast/AG_Forslund/Lola/CIPF_2018/Analysis/Rats_HE_CCL4model/output/Figures/All_variables_blockplots.svg", heatmap_all, width = 12.38, height = 7.38, device = "svg")  # Adjust width and height as needed
