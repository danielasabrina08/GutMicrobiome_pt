#Cargar librerias
library("phyloseq")
library("ggplot2")
library("metagenomeSeq")
library("vegan")
library("MASS")
library("RColorBrewer")
library("dplyr")

#Ubicar la carpeta
setwd("C:/Users/Daniela Rivera/Desktop/Datos pt/Physeq")

#Import data function
import_data <- function(file) {
  data <- as.matrix(read.csv(file, header = T, row.names = NULL, sep = ",", check.names = F))
  return(data)
}

#assign rownames function
assign_dlt <- function(data) {
  rownames(data) <- data[,1]
  return(data[,-1])
}

otu.out <- import_data("OTUtable.csv")
taxa <- assign_dlt(otu.out)
otu_table <- taxa

#convert to numeric
otu.table <- matrix(
  as.numeric(otu_table), ncol = ncol(otu_table))

#assign names to numeric table
rownames(otu.table) <- rownames(taxa)
colnames(otu.table) <- colnames(otu_table)

#import data from taxonomy
taxa.table <- import_data("taxa.csv")
taxa.full <- assign_dlt(taxa.table)

#add sample data
sample.out <- as.data.frame(import_data("Metadata.full.csv"))
sample.out <- assign_dlt(sample.out)


#prepare to make phyloseq object
OTU = otu_table(otu.table, taxa_are_rows = TRUE)
TAX = tax_table(taxa.full)

#prepare a phyloseq_object
samplingdata <- sample_data(sample.out)
physeq = phyloseq(OTU, TAX, samplingdata)

#Phyloseq object with merged data
physeq_mg <- merge_samples(physeq, "percentile")

#Plot relative abundance
p <- plot_bar(physeq_mg, fill = "Species")
png("div_sp.png", width = 18*200, height = 12*200, res = 400, pointsize = 6)
p + geom_bar(aes(fill = Species, color = Species), stat = "identity") +
  theme(legend.position = "none")
dev.off()
#Generar tabla con alpha-indices
diversity_index<-as.data.frame(round(estimate_richness(physeq, measures=c("Shannon", "Simpson", "InvSimpson")), 4))
diversity_index2<-data.frame(Sample=rownames(diversity_index))
diversity_index3<-cbind(diversity_index2, diversity_index)
write.csv(diversity_index3, "alpha_diversity_index_vf.csv", quote = FALSE, sep = ",", row.names = FALSE, col.names = TRUE)

#Cálculos manuales de los indices de diversidad
Shannon.index <- estimate_richness(physeq_mg, measures = "Shannon")

#Gráfico diversidad alfa, agrupado por percentiles y en formato de boxplot
png("alpha_indexesX.png", width = 18*200, height = 12*200, res = 400, pointsize = 6)
plot_richness(physeq, x = "percentile", color = "percentile", measures=c("Shannon", "Simpson", "InvSimpson")) + 
              geom_boxplot(aes(group = percentile, fill = percentile), alpha = 0.8, color = "black") + 
              scale_color_manual(values = brewer.pal(3, name = "Accent")) +
              scale_fill_manual(values = brewer.pal(3, name = "Accent"))
dev.off()


#2. Estimate Bray-curtis
Bray_distance_otus <- vegdist(transp_otu_mat, method="bray")
#Make_clusters with UPGMA
otus.hclust <- hclust(Bray_distance_otus, method="average")
#Plot UPGMA
pdf("UPGMA_plot_average.pdf", width=7, height=7)
plot(otus.hclust, hang=-1, main="UPGMA Bray Curtis OTUs")
dev.off()

#Beta-diversity
mds.bray <- ordinate(physeq, method = "MDS", distance = "bray")
evals <- mds.bray$values$Eigenvalues

#Graph: MDS w Bray distance
png("beta_divv5.png", width = 18*200, height = 12*200, res = 400, pointsize = 6)
plot_ordination(physeq, mds.bray, color = "percentile") +
  labs(col = "percentile") +
  geom_point(size = 4) +
  scale_color_manual(values = brewer.pal(3, name = "Accent")) +
  scale_fill_manual(values = brewer.pal(3, name = "Accent")) +
  coord_fixed(sqrt(evals[2] / evals[1]))
dev.off()

#Kruskal-Wallis test
matrix_shannon <- import_data("alphadv_shannon.csv")
df_shannon <- as.data.frame(matrix_shannon)
kruskal.test(Shannon ~Group, data = df_shannon)
