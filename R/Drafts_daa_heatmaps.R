#HeatMaps of significant taxa
###Coin
Dif abundance coin genus
```{r Heatmap 1}
library("RColorBrewer")
ps.t<- ps_0 %>% subset_samples(Group %in% c("ccl4+rif", "ccl4"))
matrix_1 <- as.matrix(data.frame(otu_table(ps.t)))
matrix_2 <- subset(t(matrix_1), select = c("ASV483", "ASV114", "ASV208", "ASV115", "ASV134", "ASV583", "ASV127", "ASV226", "ASV19", "ASV252", "ASV274", "ASV866", "ASV367", "ASV719", "ASV6", "ASV462", "ASV109", "ASV359", "ASV298", "ASV690"))

metadata_sub <- data.frame(sample_data(ps_0))
matrix <- t(matrix_2)


my_group <- as.numeric(as.factor(meta(ps.t)$Group))
colSide <- t(brewer.pal(8, "Set1")[my_group])
colMain <- colorRampPalette(brewer.pal(8, "Blues"))(16)
heatmap(matrix, Colv = NA, Rowv = NA, scale="row" , ColSideColors=(colSide), col=colMain, xlab="SampleID", ylab="ASV", main="Dif abundance coin genus ccl4+rif vs ccl4")


heatmap(matrix, Colv = NA, Rowv = NA, scale="row", xlab="SampleID", ylab="ASV", main="Dif abundance coin genus ccl4+rif vs ccl4")
```

Dif abundance coin genus
```{r Heatmap 2}
library("RColorBrewer")
ps.t<- ps_0 %>% subset_samples(Group %in% c("ctr", "ccl4"))
matrix_1 <- as.matrix(data.frame(otu_table(ps.t)))
matrix_2 <- subset(t(matrix_1), select = c("ASV358", "ASV583", "ASV18", "ASV127", "ASV58", "ASV213"))

metadata_sub <- data.frame(sample_data(ps_0))
matrix <- t(matrix_2)


my_group <- as.numeric(as.factor(meta(ps.t)$Group))
colSide <- t(brewer.pal(8, "Blues")[my_group])
colMain <- colorRampPalette(brewer.pal(2, "Blues"))(16)
heatmap(matrix, Colv = NA, Rowv = NA, scale="row" , ColSideColors=(colSide), col=colMain, xlab="SampleID", ylab="ASV", main="Dif abundance coin genus ctr vs ccl4")


heatmap(matrix, Colv = NA, Rowv = NA, scale="row", xlab="SampleID", ylab="ASV", main="
        Dif abundance coin genus ctr vs ccl4")
```


###mvabund
Dif abundance mvabund ASVs
```{r Heatmap 3}
library("RColorBrewer")
ps.t<- ps_0 %>% subset_samples(Group %in% c("ctr", "ccl4"))
matrix_1 <- as.matrix(data.frame(otu_table(ps.t)))
matrix_2 <- subset(t(matrix_1), select = c("ASV36", "ASV88", "ASV103", "ASV394", "ASV412", "ASV512", "ASV205", "ASV17", "ASV45", "ASV74", "ASV118", "ASV123", "ASV228", "ASV256", "ASV400", "ASV286", "ASV379", "ASV358", "ASV99", "ASV172", "ASV290", "ASV11", "ASV530", "ASV143", "ASV547", "ASV570", "ASV461", "ASV349", "ASV306", "ASV124", "ASV296", "ASV330", "ASV341", "ASV441", "ASV683", "ASV717", "ASV538", "ASV70", "ASV625", "ASV526", "ASV254", "ASV346", "ASV736", "ASV607", "ASV201", "ASV712", "ASV698", "ASV721", "ASV64", "ASV414"))

metadata_sub <- data.frame(sample_data(ps_0))
matrix <- t(matrix_2)


my_group <- as.numeric(as.factor(meta(ps.t)$Group))
colSide <- t(brewer.pal(8, "Set1")[my_group])
colMain <- colorRampPalette(brewer.pal(8, "Blues"))(16)
heatmap(matrix, Colv = NA, Rowv = NA, scale="row" , ColSideColors=(colSide), col=colMain, xlab="SampleID", ylab="ASV", main="Dif abundance mvabund ASVs ctr vs ccl4" )


heatmap(matrix, Colv = NA, Rowv = NA, scale="row", ColSideColors=(colSide), xlab="SampleID", ylab="ASV", main="Dif abundance mvabund ASVs ctr vs ccl4")
```

#Dif abundance mvabund ASVs
```{r Heatmap 4}
library("RColorBrewer")
ps.t<- ps_0 %>% subset_samples(Group %in% c("ctr", "ccl4"))
matrix_1 <- as.matrix(data.frame(otu_table(ps.t)))%>% microbiome::transform("compositional")
matrix_2 <- subset(t(matrix_1), select = c("ASV36","ASV88","ASV103","ASV394","ASV412","ASV512","ASV205","ASV17","ASV45","ASV74","ASV118","ASV123","ASV228","ASV256","ASV400","ASV286","ASV379","ASV358","ASV99","ASV172","ASV290","ASV11","ASV530","ASV143","ASV547","ASV570","ASV461","ASV349","ASV306","ASV124","ASV296","ASV330","ASV341","ASV441","ASV683","ASV717","ASV538","ASV70","ASV625","ASV526","ASV254", "ASV346","ASV736","ASV607", "ASV201", "ASV712", "ASV698", "ASV721", "ASV64", "ASV414"))

metadata_sub <- data.frame(sample_data(ps_0))
matrix <- t(matrix_2)


my_group <- as.numeric(as.factor(meta(ps.t)$Group))
colSide <- t(brewer.pal(8, "Set1")[my_group])
colMain <- colorRampPalette(brewer.pal(8, "Blues"))(16)
heatmap(matrix, Colv = NA, Rowv = NA, scale="row" , ColSideColors=(colSide), col=colMain, xlab="SampleID", ylab="ASV", main="Dif abundance mvabund ASV ccl4+rif vs ccl4")
```

#Correlations
```{r Heatmap 1}
library(microbiome) # Load libraries
library(phyloseq)
library(dplyr)
library(reshape2)
library(knitr)

#If you wanted to select taxa and group
pseq2 <- ps_0 %>%
  subset_taxa(Phylum == "Bacteroidetes") %>%
  subset_samples(group == "LGG")


# Z transformed abundance data
pseqz <- microbiome::transform(ps_0, "Z")

# Load example data 
ps_3 <- tax_glom(ps_0, taxrank = 'Genus', NArm = FALSE) 
otu <- as.matrix(otu_table(pseqz))

# Define data sets to cross-correlate
metadata1 <- subset(metadata, select = -c(SampleID, Group, Antibiotic, Treated_CCL, FisabioID, Week, Batch))
metadata1 <- as.matrix(metadata1) 
rownames(metadata1) <- (metadata$SampleID)
x <- log10(t(otu))
y <- as.matrix(metadata1) 

# Cross correlate data sets
correlations <- associate(x, y, method = "spearman", mode = "matrix", p.adj.threshold = 0.05, n.signif = 1)

# Or, alternatively, the same output is also available in a handy table format
correlation.table <- associate(x, y, method = "spearman", mode = "table", p.adj.threshold = 0.05, n.signif = 1)

kable(head(correlation.table))
p <- heat(correlation.table, "X1", "X2", 
          fill = "Correlation", 
          star = "p.adj", 
          p.adj.threshold = 0.05) 

p + theme(text=element_text(size=10), 
          axis.text.x = element_text(angle = 90, hjust = 1), 
          legend.key.size = unit(1.3, "cm"))
```

```{r Heatmap 2}
# Z transformed abundance data
pseqz <- microbiome::transform(ps_0, "Z")
otu <- as.matrix(otu_table(pseqz))

# Define data sets to cross-correlate
metadata1 <- subset(metadata, select = -c(SampleID, Group, Antibiotic, Treated_CCL, FisabioID, Week, Batch))
metadata1 <- as.matrix(metadata1) 
rownames(metadata1) <- (metadata$SampleID)
x <- as.data.frame(t(otu))
y <- as.data.frame(metadata1) 
# Set black-and-white theme
library(ggplot2)
theme_set(theme_bw())

# Pick only the correlations with q<0.05
# Note: this will leave other cells empty
library(dplyr)
correlation.table <- associate(x, y, method = "spearman", mode = "table", p.adj.threshold = 0.95, n.signif = 1)
subtable <- filter(correlation.table, p.adj < 0.6)

# Arrange the figure
p <- ggplot(subtable, aes(x = X1, y = X2, fill = Correlation))
p <- p + geom_tile() 
p <- p + scale_fill_gradientn("Correlation", 
                              breaks = seq(from = -1, to = 1, by = 0.2), 
                              colours = c("darkblue", "blue", "white", "red", "darkred"), 
                              limits = c(-1,1)) 

# Polish texts
p <- p + theme(axis.text.x=element_text(angle = 90, hjust=1, face = "italic"),
               axis.text.y=element_text(size = 8))
p <- p + xlab("") + ylab("")

# Mark the most significant cells with stars
p <- p + geom_text(data = subset(correlation.table, p.adj < 0.1), 
                   aes(x = X1, y = X2, label = "+"), col = "white", size = 5)

# Plot
print(p)
```
