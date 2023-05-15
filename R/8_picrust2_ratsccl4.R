##Aim: Generate input for picrust2
###Otu table in biomformat 

require("biomformat")

library(phyloseq); packageVersion("phyloseq")
library(ggplot2); packageVersion("ggplot2")
library(readr)
library(readxl)
library(readr)
library(Biostrings)


otu_table <- read_delim("~/Documents/Lola/CIPF_2018/Rats_HE_CCL4model/output/Tables/otu_table_with_seq.csv", 
                        delim = ";", escape_double = FALSE, trim_ws = TRUE)
dna<- readDNAStringSet("~/Documents/Lola/CIPF_2018/Rats_HE_CCL4model/output/ASV.fasta")


otu_table <- as.data.frame(otu_table)
otu_table <-na.omit(otu_table)
rownames(otu_table) <- otu_table$...1
otu_ok <- subset(t(otu_table), select = -c(CV1_S1, CV16_S2, CR23_S3, CR24_S4, CC34_36_S5, CC37_1_S6, CCR55_S7, CCR56_S8))
otu_ok <- subset(t(otu_ok), select = -c( ...1))


##Biom file asv table
asvmat<- as.matrix(otu_ok)
asvmat <- t(asvmat) ##Rows should be the ASVs and columns the samples

tmp<- as.data.frame(rownames(asvmat))
tmp[,2]<- paste0("ASV", 1:nrow(asvmat))
colnames(tmp)<- c("Sequence", "ASV")
rownames(asvmat) <- paste0("ASV", 1:nrow(asvmat))
biom.tmp<- make_biom(asvmat, matrix_element_type = "int")
write_biom(biom.tmp, biom_file = "~/Documents/Lola/CIPF_2018/Rats_HE_CCL4model/input/asv_ratsccl4.biom")


##Sequences with corrected names
library(DECIPHER)

dna <- readDNAStringSet("~/Documents/Lola/CIPF_2018/Rats_HE_CCL4model/output/ASV.fasta")
names(dna)<- paste0("ASV", 1:length(dna)) 
writeXStringSet(dna, filepath = "~/Documents/Lola/CIPF_2018/Rats_HE_CCL4model/input/dna_ratsccl4.fasta") 