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
# Convert 'featureNames' to a factor and reorder it based on 'Family' to order the rows by Family
Metadec_table$featureNames <- factor(Metadec_table$featureNames, levels = unique(Metadec_table$featureNames[order(Metadec_table$Family)]))

heatmap<-ggplot(Metadec_table, aes(x = metaVariable, y = reorder(featureNames, ordernames1))) + geom_tile(aes(fill = Ds)) +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0) +
  geom_text(aes(label = stars, colour = status), size = 5, key_glyph = "point") +
  scale_color_manual(values = c("black", "gray22"), labels = c("Confounded", "Deconfounded")) +
  guides(color = guide_legend(override.aes = list(shape = c(1, 8)))) +
  theme_classic()  +
  theme(
    axis.text.x = element_text(size = 12, angle = 90, hjust = 1, vjust = 0.3),
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    plot.title.position = "plot",
    plot.title = element_text(hjust = 0),
    plot.subtitle = element_text(size = 8) ,
    strip.text = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12)
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

scatter_plot <- ggplot(data = Metadec_table, aes(x = "Family", y = reorder(featureNames, ordernames1), color = Family)) +
  geom_point() +
  xlab("Family") +
  ylab("Feature Value")+ theme(legend.position = "bottom", panel.background = element_blank(),
                               legend.text = element_text(size = 10),  # Adjust the legend text size if needed
                               legend.title = element_text(size = 12),  # Adjust the legend title size if needed
                               legend.box.spacing = unit(0.2, "cm"),  # Adjust the spacing between legend items if needed
                               legend.spacing = unit(0.2, "cm"),  # Adjust the spacing between legend columns if needed
                               legend.direction = "vertical",  # Set the direction to vertical
                               legend.box = "vertical",  # Set the box layout to vertical
                               )+
  guides(color = guide_legend(ncol = 5)  )
print(scatter_plot)
#Change relwith to 1.5 if we keep the legend of scater_plot
plot_fam<-cowplot::plot_grid(scatter_plot, heatmap, align = "h", axis = "ltb", rel_widths = c(1,3))
plot_fam
ggsave("/fast/AG_Forslund/Lola/CIPF_2018/Analysis/Rats_HE_CCL4model/output/Figures/plot_fam.svg", plot_fam, width = 12.38, height = 7.38, device = "svg")  # Adjust width and height as needed

#New aesthetics
# Replace specific substrings in featureNames based on feature column
Metadec_table$featureNames <- ifelse(Metadec_table$feature == "ASV36 NK4A214 group", 
                                     "ASV36 Oscillospiraceae NK4A214 group", 
                                     ifelse(Metadec_table$feature == "ASV196 UBA1819", 
                                            "ASV196 Ruminococcaceae UBA1819",
                                            ifelse(Metadec_table$feature == "ASV538 A2", 
                                                   "ASV538 Lachnospiraceae UBA1819",
                                                   Metadec_table$feature)))
# Replace NAs in Family based on feature column
Metadec_table$Family <- ifelse(Metadec_table$feature == "ASV31 Clostridia UCG-014", 
                                     "Clostridia UCG-014", 
                                     ifelse(Metadec_table$feature == "ASV328 Clostridia", 
                                            "Order Clostridia; Family NA",
                                            ifelse(Metadec_table$feature == "ASV335 Clostridia vadinBB60 group", 
                                                   "Clostridia vadinBB60 group",
                                                   ifelse(Metadec_table$feature == "ASV326 Gastranaerophilales", 
                                                                                       "Order Gastranaerophilales; Family NA",
                                                   Metadec_table$Family))))
# Replace substrings in the metaVariable column
Metadec_table$metaVariable <- gsub("SCFA_AA", "Acetic Acid", Metadec_table$metaVariable)
Metadec_table$metaVariable <- gsub("SCFA_BA", "Butyric Acid", Metadec_table$metaVariable)
Metadec_table$metaVariable <- gsub("SCFA_CA", "Caproic Acid", Metadec_table$metaVariable)
Metadec_table$metaVariable <- gsub("Treated_CCL", "Treated CCl4", Metadec_table$metaVariable)
Metadec_table$Blocks <- gsub("Treatments", "Treat.", Metadec_table$Blocks)
Metadec_table$Blocks <- gsub("Memb. receptors", "Receptors", Metadec_table$Blocks)
Metadec_table$Blocks <- gsub("Cytokines", "Cyt", Metadec_table$Blocks)
Metadec_table$Blocks <- gsub("Cognition", "Cog", Metadec_table$Blocks)

# Remove ASV number from featureNames
Metadec_table$featureNames <- gsub("^ASV[0-9]+\\s", "", Metadec_table$featureNames)
# Convert featureNames back to factor if necessary
Metadec_table$featureNames <- as.factor(Metadec_table$featureNames)
# Convert 'featureNames' to a factor and reorder it based on 'Family' to order the rows by Family
Metadec_table$featureNames <- factor(Metadec_table$featureNames, levels = unique(Metadec_table$featureNames[order(Metadec_table$Family)]))

heatmap<-ggplot(Metadec_table, aes(x = reorder(featureNames, ordernames1), y = metaVariable)) + geom_tile(aes(fill = Ds)) +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0) +
  geom_text(aes(label = stars, colour = status), size = 6.5, key_glyph = "point") +
  scale_color_manual(values = c("black", "gray22"), labels = c("Confounded", "Deconfounded")) +
  guides(color = guide_legend(override.aes = list(shape = c(1, 8)))) +
  theme_classic()  +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 12),
    axis.title.y = element_blank(),
    plot.title.position = "plot",
    plot.title = element_text(hjust = 0),
    plot.subtitle = element_text(size = 8) ,
    strip.text = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12)
  ) +
  labs(
    title = "Significant metadata correlation heatmap",
    subtitle = "FDR-values: < 0.001 = ***, < 0.01 = **, < 0.1 = *",
    x = "Variables",
    y = "ASVs agglomerated by Genus",
    fill = "Effect size", 
    shape = "Confounding status"
  )+
  facet_grid(Blocks~., scales = "free",space='free', switch = "y")
scatter_plot <-ggplot(data = Metadec_table, aes(x = reorder(featureNames, ordernames1), y = "Family", color = Family)) +
  geom_point() +
  xlab("Family") +
  ylab("Feature Value")+ theme(legend.position = "bottom", panel.background = element_blank(),
                               legend.text = element_text(size = 10),  # Adjust the legend text size if needed
                               legend.title = element_text(size = 12),  # Adjust the legend title size if needed
                               legend.box.spacing = unit(0.2, "cm"),  # Adjust the spacing between legend items if needed
                               legend.spacing = unit(0.2, "cm"),  # Adjust the spacing between legend columns if needed
                               legend.direction = "vertical",  # Set the direction to vertical
                               legend.box = "vertical",  # Set the box layout to vertical
                               axis.text.x = element_text(size = 12, angle = 45, hjust = 1, vjust = 0.999),
                               axis.text.y = element_text(size = 12))+
  guides(color = guide_legend(ncol = 5) )
plot_fam<-cowplot::plot_grid(heatmap, scatter_plot, align = "v", axis = "tblr", rel_widths = c(1,3), nrow = 2)
plot_fam
ggsave("~/plot_fam.svg", plot_fam, width = 13, height = 9, device = "svg")  # Adjust width and height as needed


