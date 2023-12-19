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

# Make a venn diagram from the lists
library(ggvenn)
Venn_diagram<-ggvenn(x, show_elements = T, digits = 5, label_sep = "\n", fill_color =  c(("red"), ('turquoise')), fill_alpha = 0.45, set_name_color = c("darkred", 'blue'), set_name_size = 8, text_size = 2, stroke_linetype = 1, stroke_size = 0.009) + theme(plot.tag = element_text(face = "bold"))

# Save the plot as an SVG file
ggsave("VennDiagram_antibioticeffect.svg", plot = Venn_diagram, device = "svg")

