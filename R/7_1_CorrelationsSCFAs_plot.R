library(readxl)
metadata <- read_excel("OneDrive - Universitat de Valencia/Doctorado/Estudios/CIPF_Estudio1_2018/Rats_HE_CCL4model/metadata.xlsx")
View(metadata)
library(tidyr)
library(dplyr)
newmetadata <- metadata %>%
  subset(select=c(SCFA_AA, SCFA_BA, SCFA_PA, SCFA_CA, SCFA_VA, Group)) %>% data.frame()%>%
  rename(
    "Acetic acid" = SCFA_AA,
    "Butyric acid" = SCFA_BA,
    "Propionic acid" = SCFA_PA,
    "Caproic acid" = SCFA_CA,
    "Valeric acid" = SCFA_VA
  )
#Associations Controls
library(corrplot)
filtered_data_ctr <- subset(newmetadata, Group == "ctr")
filtered_data_ctr<-filtered_data_ctr[, c("Acetic acid", "Butyric acid", "Propionic acid", "Caproic acid", "Valeric acid")]
M = cor(filtered_data_ctr)
testRes = cor.mtest(filtered_data_ctr, conf.level = 0.95)
## leave blank on non-significant coefficient
## add significant correlation coefficients
tiff("~/OneDrive - Universitat de Valencia/Doctorado/Estudios/CIPF_Estudio1_2018/Rats_HE_CCL4model/output/Figures/SCFAs_final/control1.tiff", width = 6, height = 6, units = 'in', res = 600, compression = "lzw")
corrplot(M, p.mat = testRes$p, method = 'circle', type = 'lower', insig='blank',
         addCoef.col ='white', number.cex = 0.8, order = 'alphabet', diag=FALSE, tl.col = 'black', tl.cex = 1.2, main = "Control")

dev.off()
## add significant level stars
tiff("~/OneDrive - Universitat de Valencia/Doctorado/Estudios/CIPF_Estudio1_2018/Rats_HE_CCL4model/output/Figures/SCFAs_final/control2.tiff", width = 6, height = 6, units = 'in', res = 600, compression = "lzw")
corrplot(M, p.mat = testRes$p, diag = FALSE, type = 'upper',method = 'circle',
         sig.level = c(0.001, 0.01, 0.05), pch.cex = 1.9,
         insig = 'label_sig', pch.col = 'white', order = 'alphabet', tl.col = 'black', tl.cex = 1.2, main = "Control")

dev.off()


#Associations Controls+Rifaximin
library(corrplot)
filtered_data_ctrrif <- subset(newmetadata, Group == "ctr+rif")
filtered_data_ctrrif <- na.omit(filtered_data_ctrrif)
filtered_data_ctrrif<-filtered_data_ctrrif[, c("Acetic acid", "Butyric acid", "Propionic acid", "Caproic acid", "Valeric acid")]
M = cor(filtered_data_ctrrif)
testRes = cor.mtest(filtered_data_ctrrif, conf.level = 0.95)
## leave blank on non-significant coefficient
## add significant correlation coefficients
tiff("~/OneDrive - Universitat de Valencia/Doctorado/Estudios/CIPF_Estudio1_2018/Rats_HE_CCL4model/output/Figures/SCFAs_final/controlrif1.tiff", width = 6, height = 6, units = 'in', res = 600, compression = "lzw")
corrplot(M, p.mat = testRes$p, method = 'circle', type = 'lower', insig='blank',
         addCoef.col ='white', number.cex = 0.8, order = 'alphabet', diag=FALSE, tl.col = 'black', tl.cex = 1.2, main = "Control+Rifaximin")

dev.off()

## add significant level stars
tiff("~/OneDrive - Universitat de Valencia/Doctorado/Estudios/CIPF_Estudio1_2018/Rats_HE_CCL4model/output/Figures/SCFAs_final/controlrif2.tiff", width = 6, height = 6, units = 'in', res = 600, compression = "lzw")
corrplot(M, p.mat = testRes$p, diag = FALSE, type = 'upper',method = 'circle',
         sig.level = c(0.001, 0.01, 0.05), pch.cex = 1.9,
         insig = 'label_sig', pch.col = 'white', order = 'alphabet', tl.col = 'black', tl.cex = 1.2, main = "Control+Rifaximin")

dev.off()


#Associations CCl4
library(corrplot)
filtered_data_ccl4 <- subset(newmetadata, Group == "ccl4")
filtered_data_ccl4<-filtered_data_ccl4[, c("Acetic acid", "Butyric acid", "Propionic acid", "Caproic acid", "Valeric acid")]
M = cor(filtered_data_ccl4)
testRes = cor.mtest(filtered_data_ccl4, conf.level = 0.95)
## leave blank on non-significant coefficient
## add significant correlation coefficients
tiff("~/OneDrive - Universitat de Valencia/Doctorado/Estudios/CIPF_Estudio1_2018/Rats_HE_CCL4model/output/Figures/SCFAs_final/ccl4_1.tiff", width = 6, height = 6, units = 'in', res = 600, compression = "lzw")
corrplot(M, p.mat = testRes$p, method = 'circle', type = 'lower', insig='blank',
         addCoef.col ='white', number.cex = 0.8, order = 'alphabet', diag=FALSE, tl.col = 'black', tl.cex = 1.2, main = "CCl4")

dev.off()

## add significant level stars
tiff("~/OneDrive - Universitat de Valencia/Doctorado/Estudios/CIPF_Estudio1_2018/Rats_HE_CCL4model/output/Figures/SCFAs_final/cccl4_2.tiff", width = 6, height = 6, units = 'in', res = 600, compression = "lzw")
corrplot(M, p.mat = testRes$p, diag = FALSE, type = 'upper',method = 'circle',
         sig.level = c(0.001, 0.01, 0.05), pch.cex = 1.9,
         insig = 'label_sig', pch.col = 'white', order = 'alphabet', tl.col = 'black', tl.cex = 1.2, main = "CCl4")

dev.off()


#Associations CCl4+Rif
library(corrplot)
filtered_data_ccl4rif <- subset(newmetadata, Group == "ccl4+rif")
filtered_data_ccl4rif<-filtered_data_ccl4rif[, c("Acetic acid", "Butyric acid", "Propionic acid", "Caproic acid", "Valeric acid")]
M = cor(filtered_data_ccl4rif)
testRes = cor.mtest(filtered_data_ccl4rif, conf.level = 0.95)
## leave blank on non-significant coefficient
## add significant correlation coefficients
tiff("~/OneDrive - Universitat de Valencia/Doctorado/Estudios/CIPF_Estudio1_2018/Rats_HE_CCL4model/output/Figures/SCFAs_final/ccl4rif_1.tiff", width = 6, height = 6, units = 'in', res = 600, compression = "lzw")
corrplot(M, p.mat = testRes$p, method = 'circle', type = 'lower', insig='blank',
         addCoef.col ='white', number.cex = 0.8, order = 'alphabet', diag=FALSE, tl.col = 'black', tl.cex = 1.2, main = "CCl4+Rifaximin")
dev.off()
## add significant level stars
tiff("~/OneDrive - Universitat de Valencia/Doctorado/Estudios/CIPF_Estudio1_2018/Rats_HE_CCL4model/output/Figures/SCFAs_final/cccl4rif_2.tiff", width = 6, height = 6, units = 'in', res = 600, compression = "lzw")
corrplot(M, p.mat = testRes$p, diag = FALSE, type = 'upper',method = 'circle',
              sig.level = c(0.001, 0.01, 0.05), pch.cex = 1.9,
              insig = 'label_sig', pch.col = 'white', order = 'alphabet', tl.col = 'black', tl.cex = 1.2, main = "CCl4+Rifaximin")
dev.off()

#Boxplots SCFAs
library(ggpubr)
# Acetic acid
AA<-ggplot(newmetadata, aes(x = Group, y = `Acetic acid`)) +   geom_boxplot() + 
  geom_jitter(aes(color = Group), position = position_jitterdodge(dodge.width = 0.65, jitter.width = 0.1, jitter.height = 0.1)) +
  theme_classic()+ theme(strip.background = element_blank(), axis.text.x.bottom = element_text(angle = -90), axis.text = element_text(size = 16), axis.title = element_text(size = 16), legend.text = element_text(size = 16), legend.title = element_text(size = 16))  +
  geom_pwc(label = "{p.adj}{p.adj.signif}", hide.ns = TRUE, method = "wilcox_test", p.adjust.method = "BH") 

# Butyric acid
BA<-ggplot(newmetadata, aes(x = Group, y = `Butyric acid`)) +   geom_boxplot() + 
  geom_jitter(aes(color = Group), position = position_jitterdodge(dodge.width = 0.65, jitter.width = 0.1, jitter.height = 0.1)) +
  theme_classic()+ theme(strip.background = element_blank(), axis.text.x.bottom = element_text(angle = -90), axis.text = element_text(size = 16), axis.title = element_text(size = 16), legend.text = element_text(size = 16), legend.title = element_text(size = 16))  +
  geom_pwc(label = "{p.adj}{p.adj.signif}", hide.ns = TRUE, method = "wilcox_test", p.adjust.method = "BH") 

# Caproic acid
CA<-ggplot(newmetadata, aes(x = Group, y = `Caproic acid`)) +   geom_boxplot() + 
  geom_jitter(aes(color = Group), position = position_jitterdodge(dodge.width = 0.65, jitter.width = 0.1, jitter.height = 0.1)) +
  theme_classic()+ theme(strip.background = element_blank(), axis.text.x.bottom = element_text(angle = -90), axis.text = element_text(size = 16), axis.title = element_text(size = 16), legend.text = element_text(size = 16), legend.title = element_text(size = 16))  +
  geom_pwc(label = "{p.adj}{p.adj.signif}", hide.ns = TRUE, method = "wilcox_test", p.adjust.method = "BH") 

# Propionic acid
PA<-ggplot(newmetadata, aes(x = Group, y = `Propionic acid`)) +   geom_boxplot() + 
  geom_jitter(aes(color = Group), position = position_jitterdodge(dodge.width = 0.65, jitter.width = 0.1, jitter.height = 0.1)) +
  theme_classic()+ theme(strip.background = element_blank(), axis.text.x.bottom = element_text(angle = -90), axis.text = element_text(size = 16), axis.title = element_text(size = 16), legend.text = element_text(size = 16), legend.title = element_text(size = 16))  +
  geom_pwc(label = "{p.adj}{p.adj.signif}", hide.ns = TRUE, method = "wilcox_test", p.adjust.method = "BH") 

# Valeric acid
VA<-ggplot(newmetadata, aes(x = Group, y = `Valeric acid`)) +   geom_boxplot() + 
  geom_jitter(aes(color = Group), position = position_jitterdodge(dodge.width = 0.65, jitter.width = 0.1, jitter.height = 0.1)) +
  theme_classic()+ theme(strip.background = element_blank(), axis.text.x.bottom = element_text(angle = -90), axis.text = element_text(size = 16), axis.title = element_text(size = 16), legend.text = element_text(size = 16), legend.title = element_text(size = 16))  +
  geom_pwc(label = "{p.adj}{p.adj.signif}", hide.ns = TRUE, method = "wilcox_test", p.adjust.method = "BH") 

ggarrange(AA, BA, CA, PA, VA, common.legend = TRUE, legend="bottom",nrow = 2)

ggsave(filename = "~/OneDrive - Universitat de Valencia/Doctorado/Estudios/CIPF_Estudio1_2018/Rats_HE_CCL4model/output/Figures/SCFAs_final/Boxplots.tiff")

