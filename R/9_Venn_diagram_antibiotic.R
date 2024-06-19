library(eulerr)
library(microbiome)
library(tidyverse)
library(microbiomeutilities)
library(VennDiagram)

Output2_batch_asRandVar_ctrvsctrrif <- readRDS("~/path/Rats_HE_CCL4model/output/MetadeconfoundR/Output2_batch_asRandVar_ctrvsctrrif.Rds")
Output2_batch_ccl4vsccl4rif <- readRDS("~/path/Rats_HE_CCL4model/output/MetadeconfoundR/Output2_batch_ccl4vsccl4rif.Rds")

ccl4_status<-data.frame(Output2_batch_ccl4vsccl4rif[["status"]])
ctr_status<-data.frame(Output2_batch_asRandVar_ctrvsctrrif[["status"]])

ccl4_status<-ccl4_status[,c("Antibiotic", "Group")]
ctr_status<-ctr_status[,c("Antibiotic", "Ammonia")]
# Remove rows with "NS" in the "Antibiotic" column
ctr_status <- ctr_status[ctr_status$Antibiotic != "NS", ]
ccl4_status <- ccl4_status[ccl4_status$Antibiotic != "NS", ]

# Create a list of row names
ctr_status$list <- rownames(ctr_status)
ccl4_status$list <- rownames(ccl4_status)
ctr_asvs <- (ctr_status$list)
ccl4_asvs <- (ccl4_status$list)

# Make a venn diagram from the lists
library(ggvenn)
x<-(list(Control = ctr_asvs, CCl4 = ccl4_asvs))
Venn_diagram<-ggvenn(x, show_elements = T, digits = 5, label_sep = "\n", fill_color =  c(("red"), ('turquoise')), fill_alpha = 0.45, set_name_color = c("darkred", 'blue'), set_name_size = 8, text_size = 2, stroke_linetype = 1, stroke_size = 0.009) + theme(plot.tag = element_text(face = "bold"))

# Save the plot as an SVG file
ggsave("VennDiagram_antibioticeffect.svg", plot = Venn_diagram, device = "svg")



#Volcano plot
library(dplyr)
library(ggplot2)
tax_table<-data.frame(tax_table(ps_0))
tax_table$ASV_numbers<-rownames(tax_table)
ccl4_status<-data.frame(Output2_batch_ccl4vsccl4rif[["status"]])
ctr_status<-data.frame(Output2_batch_asRandVar_ctrvsctrrif[["status"]])
ccl4_Ds<-data.frame(Output2_batch_ccl4vsccl4rif[["Ds"]])
ctr_Ds<-data.frame(Output2_batch_asRandVar_ctrvsctrrif[["Ds"]])

ccl4_status<-ccl4_status[,c("Antibiotic", "Group")]
ctr_status<-ctr_status[,c("Antibiotic", "Group")]
ccl4_Ds<-ccl4_Ds[,c("Antibiotic", "Group")]
ctr_Ds<-ctr_Ds[,c("Antibiotic", "Group")]

ccl4_status$CCl4_status<-ccl4_status$Antibiotic
ccl4_status$Antibiotic<-NULL
ccl4_status$Group<-NULL
ctr_status$Ctr_status<-ctr_status$Antibiotic
ctr_status$Antibiotic<-NULL
ctr_status$Group<-NULL
ccl4_Ds$CCl4<-ccl4_Ds$Antibiotic
ccl4_Ds$Antibiotic<-NULL
ccl4_Ds$Group<-NULL
ctr_Ds$Ctr<-ctr_Ds$Antibiotic
ctr_Ds$Antibiotic<-NULL
ctr_Ds$Group<-NULL

Ds_data <- data.frame(ccl4_status, ctr_status, ccl4_Ds, ctr_Ds)
#New column with ASV numbers
row_names<-rownames(Ds_data)
ASV_numbers <- sub("^([A-Za-z0-9]+) .*$", "\\1", row_names)
Ds_data$ASV_numbers <- ASV_numbers

# Filter the data frame to include only points that are NOT "NS" in either CCl4_status or Ctr_status
labeled_data <- Ds_data %>%
  filter(CCl4_status != "NS" | Ctr_status != "NS")
#New column with ASV numbers
row_names<-rownames(labeled_data)
ASV_numbers <- sub("^([A-Za-z0-9]+) .*$", "\\1", row_names)
labeled_data$ASV_numbers <- ASV_numbers
# Create a new data frame with the desired columns
labeled_data<- data.frame(
  labeled_data,
  Phylum = tax_table$Phylum[match(labeled_data$ASV_numbers, tax_table$ASV_numbers)],
  Order = tax_table$Order[match(labeled_data$ASV_numbers, tax_table$ASV_numbers)]
)
Ds_data<- data.frame(
  Ds_data,
  Phylum = tax_table$Phylum[match(Ds_data$ASV_numbers, tax_table$ASV_numbers)],
  Order = tax_table$Order[match(Ds_data$ASV_numbers, tax_table$ASV_numbers)]
)

Ds_data$Significance <- ifelse(Ds_data$CCl4_status != "NS" & Ds_data$Ctr_status == "NS", "Significant effec of Abx in CCl4",
                               ifelse(Ds_data$CCl4_status == "NS" & Ds_data$Ctr_status != "NS", "Significant effec of Abx in Ctr",
                                      ifelse(Ds_data$CCl4_status != "NS" & Ds_data$Ctr_status != "NS", "Significant effec of Abx in both", "Not significant")))
Ds_data$Alpha_values <- ifelse(Ds_data$CCl4_status != "NS" & Ds_data$Ctr_status == "NS", "Significant",
                               ifelse(Ds_data$CCl4_status == "NS" & Ds_data$Ctr_status != "NS", "Significant",
                                      ifelse(Ds_data$CCl4_status != "NS" & Ds_data$Ctr_status != "NS", "Significant", "Not significant")))
# Define the mapping of shapes to categories
shape_mapping <- c("Not significant" = 1, "Significant effec of Abx in both" = 9, "Significant effec of Abx in Ctr" = 8, "Significant effec of Abx in CCl4" = 24) # You can adjust the shapes as needed
#Phylum
plot<-ggplot(data = Ds_data, aes(x = Ctr, y = CCl4)) +
  geom_point(aes(color = Phylum, alpha = Alpha_values, shape = Significance, fill = Phylum), size = 4) +
  geom_vline(xintercept = 0, col = "black", linetype = 'solid') +
  geom_hline(yintercept = 0, col = "black", linetype = 'solid') +
  geom_abline(intercept = 0, slope = 1, col = "black", linetype = 'dashed') +
  scale_shape_manual(values = shape_mapping) + # Apply the manual shape mapping
  theme(legend.text = element_text(size = 12),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16))+ 
  theme_light()
# Save the plot as an SVG file
ggsave("Volcanoplot.svg", plot = plot, device = "svg", height = 7, width = 9, units = "in")

#Order only significant
ggplot(data = labeled_data, aes(x = Ctr, y = CCl4)) +
  geom_point(aes(color = Order, shape = Significance, fill = Order), size = 4) +
  geom_vline(xintercept = 0, col = "black", linetype = 'solid') +
  geom_hline(yintercept = 0, col = "black", linetype = 'solid') +
  geom_abline(intercept = 0, slope = 1, col = "black", linetype = 'dashed') +
  scale_shape_manual(values = shape_mapping) + # Apply the manual shape mapping
  # Apply the new color palette
  theme(legend.text = element_text(size = 12),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16)) + 
  theme_light()

#Order only significant colour blind palette
labeled_data$Significance <- ifelse(labeled_data$CCl4_status != "NS" & labeled_data$Ctr_status == "NS", "Significant effec of Abx in CCl4",
                               ifelse(labeled_data$CCl4_status == "NS" & labeled_data$Ctr_status != "NS", "Significant effec of Abx in Ctr",
                                      ifelse(labeled_data$CCl4_status != "NS" & labeled_data$Ctr_status != "NS", "Significant effec of Abx in both", "Not significant")))

# Define the new color palette
new_palette <- c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499",  "#F0E442", "#999933", "#882255", "#44AA99", "#6699CC", "#888888", "#661100", "#E69F00")
# Define the mapping of shapes to categories
shape_mapping <- c("Significant effec of Abx in both" = 9, "Significant effec of Abx in Ctr" = 19, "Significant effec of Abx in CCl4" = 17) # You can adjust the shapes as needed
# Apply the new color palette to your plot
plot <- ggplot(data = labeled_data, aes(x = Ctr, y = CCl4)) +
  geom_point(aes(color = Order, shape = Significance, fill = Order), size = 4) +
  geom_vline(xintercept = 0, col = "black", linetype = 'solid') +
  geom_hline(yintercept = 0, col = "black", linetype = 'solid') +
  geom_abline(intercept = 0, slope = 1, col = "black", linetype = 'dashed') +
  scale_shape_manual(values = shape_mapping) + # Apply the manual shape mapping
  scale_color_manual(values = new_palette) +   # Apply the new color palette
  theme(legend.text = element_text(size = 12),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16)) + 
  theme_light()

plot

#Not used
# Create a new column 'Colour' in Ds_data using nested ifelse statements
# Calculate dynamic random nudge values for each point
num_labels <- nrow(labeled_data)
nudge_x_values <- runif(num_labels, min = -0.3, max = 0.3)
nudge_y_values <- runif(num_labels, min = -0.3, max = 0.3)

Ds_data$Colour <- ifelse(Ds_data$CCl4_status != "NS" & Ds_data$Ctr_status == "NS", "blue",
                         ifelse(Ds_data$CCl4_status == "NS" & Ds_data$Ctr_status != "NS", "purple",
                                ifelse(Ds_data$CCl4_status != "NS" & Ds_data$Ctr_status != "NS", "red", "lightgray")))


library(dplyr)
library(ggplot2)
# Create ggplot with all points colored according to the 'Colour' column
ggplot(data = Ds_data, aes(x = Ctr, y = CCl4)) +
  geom_point(aes(color = Ds_data$Colour)) +
  geom_vline(xintercept = 0, col = "darkgray", linetype = 'dashed') +
  geom_hline(yintercept = 0, col = "darkgray", linetype = 'dashed') +
  geom_abline(intercept = 0, slope = 1, col = "black", linetype = 'dashed')+
  geom_text(data = labeled_data, aes(label = rownames(labeled_data),
                                     x = Ctr + nudge_x, y = CCl4 + nudge_y, 
                                     color = "black"),
            vjust = 0.1) +
  geom_segment(data = labeled_data, aes(x = Ctr, y = CCl4, xend = Ctr+ nudge_x, yend = CCl4+ nudge_y,
                                        color = "gray50")) +
  scale_color_manual(values = c("black","blue","gray50", "lightgray", "purple", "red"),
                     labels = c("Labels","Significant effect of Abx in CCl4","Segments","Other",
                                
                                "Significant effect of Abx in Ctr",
                                "Significant effect of Abx in both",
                                "")) +
  theme(
    legend.position = "bottom",
    legend.margin = margin(6, 6, 6, 6)
  )+
  labs(color = "Colour")

# Create ggplot with all points colored according to the 'Colour' column
Ds_data$Colour <- ifelse(Ds_data$CCl4_status != "NS" & Ds_data$Ctr_status == "NS", "lightblue",
                         ifelse(Ds_data$CCl4_status == "NS" & Ds_data$Ctr_status != "NS", "purple",
                                ifelse(Ds_data$CCl4_status != "NS" & Ds_data$Ctr_status != "NS", "red", "lightgray")))


# Labels withASV numbers instead of full
row_names <- rownames(labeled_data)
ASV_numbers <- sub("^([A-Za-z0-9]+) .*$", "\\1", row_names)
# Update row names of labeled_data with ASV numbers
rownames(labeled_data) <- ASV_numbers
ggplot(data = Ds_data, aes(x = Ctr, y = CCl4)) +
  geom_point(aes(color = Ds_data$Colour)) +
  geom_vline(xintercept = 0, col = "darkgray", linetype = 'dashed') +
  geom_hline(yintercept = 0, col = "darkgray", linetype = 'dashed') +
  geom_abline(intercept = 0, slope = 1, col = "black", linetype = 'dashed')+
  geom_text(data = labeled_data, aes(label = rownames(labeled_data),
                                     x = Ctr + nudge_x, y = CCl4 + nudge_y, 
                                     color = "black"),
            vjust = 0.1) +
  geom_segment(data = labeled_data, aes(x = Ctr, y = CCl4, xend = Ctr+ nudge_x, yend = CCl4+ nudge_y,
                                        color = "gray50")) +
  scale_color_manual(values = c("black","gray50","lightblue", "lightgray", "purple", "red"),
                     labels = c("Segments","Labels","Significant effect of Abx in CCl4",
                                "Other",
                                "Significant effect of Abx in Ctr",
                                "Significant effect of Abx in both",
                                "")) +
  theme(
    legend.position = "bottom",
    legend.margin = margin(6, 6, 6, 6)
  )+
  labs(color = "Colour")

ggplot(data = Ds_data, aes(x = Ctr, y = CCl4)) +
  geom_point(aes(color = Ds_data$Colour)) +
  geom_vline(xintercept = 0, col = "darkgray", linetype = 'dashed') +
  geom_hline(yintercept = 0, col = "darkgray", linetype = 'dashed') +
  geom_abline(intercept = 0, slope = 1, col = "black", linetype = 'dashed') +
  geom_text(data = labeled_data, aes(label = rownames(labeled_data)), nudge_x = -0.05, nudge_y = -0.05, check_overlap = TRUE) +
  scale_color_manual(values = c("lightblue", "lightgray", "purple", "red"),
                     labels = c("Significant effect of Abx in CCl4",
                                "Other",
                                "Significant effect of Abx in Ctr",
                                "Significant effect of Abx in both",
                                "")) +
  theme(
    legend.position = "bottom",
    legend.margin = margin(6, 6, 6, 6)
  )+
  labs(color = "Colour")

