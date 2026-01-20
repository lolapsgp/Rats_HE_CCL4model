# Rats_HE_CCL4model

### Analysis Code for:
Giner‑Pérez, L. et al. (2025). Rifaximin-induced changes in the gut microbiome associated to improvement of neurotransmission alterations and learning in rats with chronic liver disease. Scientific Reports, 15(1), 34382.
https://doi.org/10.1038/s41598-025-17229-1

### 📘 Overview
This repository contains all scripts and input files used to perform the microbiome and functional analyses described in the above publication.
The study investigated how rifaximin, a non‑absorbable antibiotic, modulates the gut microbiome and gut–brain axis in a CCl₄‑induced chronic liver disease rat model.
The repository includes: R Markdown scripts for the full analysis workflow.
Scripts for DADA2, diversity analysis, batch‑effect correction, PERMANOVA, SCFA correlations, PICRUSt2, and MetadeconfoundR.


### 📄 Study Abstract
Rifaximin, a gut-targeted antibiotic, improves cognitive function and reduces the risk of hepatic encephalopathy (HE), yet its effects on the gut–brain axis remain unknown. This study explores how rifaximin influences gut microbiota functions and its association with cognitive function and molecular alterations in rats with liver injury. Liver injury was induced by chronic administration of carbon tetrachloride (CCl₄), and rifaximin was administered daily. Fecal samples were collected after eight weeks of CCl₄ administration, and taxonomic and functional changes in the gut microbiome were analyzed.
Rifaximin altered microbiota diversity and composition, increasing α‑diversity in liver‑injured rats but reducing diversity in healthy rats. It influenced microbiota interactions with neurotransmission alterations, where Dorea, Lachnospiraceae A2, and possibly Erysipelotrichaceae might be key contributors. Functionally, butyric acid levels negatively correlated with gene orthologues associated with GABA, tryptophan, and glutamate degradation pathways.
In healthy rats, short‑chain fatty acids (SCFAs) were strongly inter‑correlated—an association absent in injured or rifaximin‑treated groups. Overall, rifaximin promoted bacterial groups linked to improved cognition and neurotransmission in liver disease, highlighting the relationship between a healthy microbiome and balanced SCFA levels.

### 📁 Repository Structure
```
Rats_HE_CCL4model/
│
├── HEAD
├── README.md
│
├── input/
│   ├── asv_ratsccl4.biom
│   ├── dna_ratsccl4.fasta
│   ├── KO_metagenome_out_ccl4/
│   └── omixerRpm/
│
└── R/
    ├── 1_Dada2_pipeline_ratsccl4.Rmd
    ├── 2_phyloseq_object.Rmd
    ├── 3_AlphaDiv.Rmd
    ├── 4_BatchEffectRemoval_and_Betadiv.Rmd
    ├── 5_MetadeconfoundR_taxa.Rmd
    ├── 5.1_plot_metadeconfoundR_all_metadata_blocks.R
    ├── 5.2_plot_metadeconfoundR_blocks_genus.R
    ├── 6_PERMANOVA_and_CDR.Rmd
    ├── 7_1_CorrelationsSCFAs_plot.R
    ├── 7_DescriptiveAnalysis_and_CorrelationsMetadata.Rmd
    ├── 8_Picrust2_analysis_modules.Rmd
    └── 9_Venn_diagram_antibiotic.R
```

### 🔍 Analysis Workflow
The analysis is organized into modular R Markdown scripts located in the /R directory:


#### 1. DADA2 Pipeline
ASV inference, filtering, chimera removal, FASTA export.


#### 2. Phyloseq Object Construction
Integration of taxonomy, ASVs, metadata, and phylogenetic tree.


#### 3. Alpha Diversity
Shannon, Simpson, Observed richness.

#### 4. Beta Diversity & Batch Correction
Bray–Curtis and UniFrac distances
Batch removal (MMUPHin / limma)
Ordination plots (PCoA, NMDS)

#### 5. MetadeconfoundR Analysis
Identification of taxa robust to confounding variables; block‑level plots.

#### 6. PERMANOVA and CDR
Variance partitioning for microbiome differences across groups.

#### 7. Metadata & SCFA Correlations
Correlation matrices, SCFA–microbe associations, statistical tests.

#### 8. PICRUSt2 Functional Analysis
KO pathways, modules, gene‑level predictions.

#### 9. Venn Diagrams
Antibiotic‑related comparisons across groups.


