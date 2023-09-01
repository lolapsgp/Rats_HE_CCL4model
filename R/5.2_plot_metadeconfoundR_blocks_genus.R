#Heatmap genus
ps_0 <- readRDS("/fast/AG_Forslund/Lola/CIPF_2018/Analysis/Rats_HE_CCL4model/output/ps_0.Rds")
ps_0_genus <- tax_glom(ps_0, taxrank = 'Genus', NArm = FALSE)
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

##Usage latest_annotations <- get_latest_annotation(phyloseq_obj)
latest_annotations <- get_latest_annotation(ps_0_genus)
latest_annotations<-as.data.frame(latest_annotations)
short_names<-as.data.frame(otu_table(ps_0_genus))
rownames(short_names)<-latest_annotations$TaxaID

# generate a vector containing the full taxonomy path for all OTUs
wholetax <- do.call(paste, c(as.data.frame(tax_table(ps_0_genus))
                             [c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")], 
                             sep = "__"))  # to distinguish from "_" within tax ranks

# turn the otu_table into a data.frame
otu_export <- as.data.frame(t(otu_table(ps_0_genus)))
tmp <- names(otu_export)

# paste wholetax and OTU_ids together
for(i in 1:length(tmp)){
  names(tmp)[i] = paste(wholetax[i], tmp[i], sep = "__")
}

# overwrite old names
names(otu_export) <- names(tmp)
otu_export<-t(otu_export)

head(otu_export)[5]


##MetadeconfoundR
library(metadeconfoundR)
library(gtools)
library(dplyr)
ps.rel<-microbiome::transform(ps_0_genus, "compositional")
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

right <- BuildHeatmap(Output_batch)
plot(right)

#Heatmap blocks
Metadec_table<-data.frame(right[["data"]])
row_names <-Metadec_table$feature
# Remove "ASV-number" prefix
cleaned_row_names <- sub("^ASV\\d+\\s", "", row_names)
Metadec_table$cleaned_row_names <- as.factor(cleaned_row_names)
Metadec_table<-Metadec_table%>%
  mutate(Blocks = case_when(
    metaVariable %in% c("Ambulatory.Counts", "Vertcal.Counts", "Stereotipic.Counts", "Learning_index_Rmaze", "Average.velocity") ~ "Cognitive tests",
    metaVariable %in% c("Antibiotic", "Treated_CCL") ~ "Treatments",
    metaVariable %in% c("Ammonia") ~ "Ammonia levels",
    metaVariable %in% c("CCL2", "CCL5","CCR5","CCL20","CCR2", "CX3CL1", "CX3CR1", "IFN.gamma", "IL_10", "IL.15", "IL.17", "IL_4", "IL.6", "Occludine",  "TNF_a", "TGF.b") ~ "Cytokines brain",
    metaVariable %in% c("GLAST", "GLT1", "GAT1", "GAT3", "TNFR1", "NR2B", "NR2A", "NR1", "GLUA1", "GLUA2") ~ "Memb. receptors",
    metaVariable %in% c("SCFA_AA", "SCFA_BA", "SCFA_CA", "SCFA_PA", "SCFA_VA") ~ "SCFAs"))%>%
  mutate(ordernames1 = case_when(
    metaVariable %in% c("Antibiotic") ~ Ds))


# Define the desired order of levels for facet_wrap
facet_order <- c("Treatments", "SCFAs", "Memb. receptors", "Cytokines brain", "Cognitive tests")

# Convert the Blocks column to a factor with the desired order
Metadec_table$Blocks <- factor(Metadec_table$Blocks, levels = facet_order)

ggplot(Metadec_table, aes(x = metaVariable, y = reorder(featureNames, ordernames1))) +
  geom_tile(aes(fill = Ds)) +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0) +
  geom_text(aes(label = stars, colour = status), size = 2, key_glyph = "point") +
  scale_color_manual(values = c("black", "gray22"), labels = c("Confounded", "Deconfounded")) +
  guides(color = guide_legend(override.aes = list(shape = c(1, 8)))) +
  theme_classic() +
  facet_wrap(vars(Metadec_table$Blocks), nrow = 1, scales = "free_x") +
  theme(
    axis.text.x = element_text(size = 7, angle = 90, hjust = 1, vjust = 0.3),
    axis.text.y = element_text(size = 7, angle = 0, hjust = 1, vjust = 0.35),
    plot.title.position = "plot",
    plot.title = element_text(hjust = 0),
    plot.subtitle = element_text(size = 8)
  ) +
  labs(
    title = "Significant metadata correlation heatmap",
    subtitle = "FDR-values: < 0.001 = ***, < 0.01 = **, < 0.1 = *",
    x = "Variables",
    y = "ASVs agglomerated by Genus",
    fill = "Effect size", 
    color = "Confounding status"
  )



#Same heatmap changing aesthetics
names<-get_latest_annotation(ps_0_genus)
families<-tax_table(ps_0_genus) %>%
      as.data.frame() %>%
     rownames_to_column('ASV') %>%select(c(ASV, Family))
# Perform a left join based on the ASV column
Families_data <- left_join(names, families, by = "ASV")
Families_data$feature<-Families_data$TaxaID
Families_data<-Families_data[,c("feature", "Family")]
Metadec_table<-left_join(Metadec_table, Families_data, by = "feature")


heatmap<-ggplot(Metadec_table, aes(x = metaVariable, y = reorder(cleaned_row_names, ordernames1))) +
  geom_tile(aes(fill = Ds)) +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0) +
  geom_text(aes(label = stars, colour = status), size = 2, key_glyph = "point") +
  scale_color_manual(values = c("black", "gray22"), labels = c("Confounded", "Deconfounded")) +
  guides(color = guide_legend(override.aes = list(shape = c(1, 8)))) +
  theme_classic()  +
  theme(
    axis.text.x = element_text(size = 7, angle = 90, hjust = 1, vjust = 0.3),
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
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
  facet_grid(. ~ (Metadec_table$Blocks), scales = "free", space='free', switch = "x") 

scatter_plot <- ggplot(data = Metadec_table, aes(x = "Family", y = reorder(cleaned_row_names, ordernames1), color = Family)) +
  geom_point() +
  xlab("Family") +
  ylab("Feature Value")+ theme(legend.position = "left", panel.background = element_blank(),
                               legend.text = element_text(size = 10),  # Adjust the legend text size if needed
                               legend.title = element_text(size = 12),  # Adjust the legend title size if needed
                               legend.box.spacing = unit(0.2, "cm"),  # Adjust the spacing between legend items if needed
                               legend.spacing = unit(0.2, "cm"),  # Adjust the spacing between legend columns if needed
                               legend.direction = "vertical",  # Set the direction to vertical
                               legend.box = "vertical",  # Set the box layout to vertical
                               )+
  guides(color = guide_legend(ncol = 1)  )
print(scatter_plot)
#Change relwith to 1.5 if we keep the legend of scater_plot
plot_fam<-cowplot::plot_grid(scatter_plot, heatmap, align = "h", axis = "ltb", rel_widths = c(1,3))
plot_fam
ggsave("/fast/AG_Forslund/Lola/CIPF_2018/Analysis/Rats_HE_CCL4model/output/Figures/plot_fam.svg", plot_fam, width = 12.38, height = 7.38, device = "svg")  # Adjust width and height as needed


#Heatmap including all variables
View(Output_batch)
raw_p <- Output_batch[1]
corr_p <- Output_batch[2]
effect_size<- Output_batch[3]
status<- Output_batch[4]

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
    comparison_p %in% c("SCFA.AA", "SCFA.BA", "SCFA.CA", "SCFA.PA", "SCFA.VA") ~ "SCFAs"))%>%
  mutate(ordernames1 = case_when(
    comparison_p %in% c("Antibiotic") ~ effectSize))%>%
  mutate(fdr= as_factor(case_when(corr_p <= 0.05 ~ "*", corr_p <= 0.01 ~ "**", corr_p <= 0.001 ~ "***", corr_p <= 0.1 ~ ".")))



#Heatmap all variables
ggplot(effect_table_sig_long, aes(x = variable, y = fct_reorder(rowname, ordernames1, .fun = function(x) mean(x, na.rm = TRUE)))) +
  # do the heatmap tile coloring based on effect sizes
  geom_tile(aes(fill = effectSize)) +
  scale_fill_gradient2 (low = "blue", high = "red", mid = "white", midpoint = 0) +
  # add significance stars/circles for deconfounded/confounded associations
  geom_text (aes (label = stars.pval(corr_p)))+
  guides(color = guide_legend(override.aes = list(shape = c(1,8)) ) ) +
  
  # make it pretty
  theme_classic() +
  facet_wrap(vars(effect_table_sig_long$Blocks), nrow =1, scales = "free_x")+
  theme(axis.text.x = element_text(size = 7,
                                   angle = 90,
                                   hjust = 1,
                                   vjust = 0.3),
        axis.text.y = element_text(size = 7,
                                   angle = 0,
                                   hjust = 1,
                                   vjust = 0.35),
        plot.title.position = "plot",
        plot.title = element_text(hjust = 0),
        plot.subtitle=element_text(size=8)) +
  labs(title="All metadata and significant genus correlation heatmap",
       subtitle="FDR-values: < 0.001 = ***, < 0.01 = **, < 0.1 = * ",
       x = "Variables",
       y = "ASVs")