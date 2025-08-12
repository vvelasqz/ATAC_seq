#atac seq analysis



#plot bars with DEG accessible and non accessible

library(tidyverse)

table_promoters_v0<- read.csv("~/iowa_state/lab/Chromatin_accessibility_paper/epistatic_inputs/table_promoters_epistasis.csv")
table_regions<- read.csv("~/iowa_state/lab/Chromatin_accessibility_paper/epistatic_inputs/table_regions_epistasis.csv")
table_any_regions<- read.csv("~/iowa_state/lab/Chromatin_accessibility_paper/epistatic_inputs/table_any_regions_epistasis.csv")




#Now do number of DE genes per timepoint 
time=c(0,16, 20, 32)
genotypes<-c("mla6", "bln1", "dm")

#padj_table<- read.csv("~/iowa_state/lab/RNAseq/RNASeq_paper/vvz manuscript/analysis/hv_R2_genes_deseq2_results_pairwise_genotype_padj_tax_sp.csv", row.names = 1)
padj_table<- read.csv("~/iowa_state/lab/RNAseq/DESeq2 HvR3/hv_R3_genes_deseq2_results_pairwise_genotype_padj_tax_sp.csv", row.names = 1)
padj_table_horvu<- padj_table[grep("HORVU", rownames(padj_table)), grep("wt", colnames(padj_table))]
padj_table_horvu<- padj_table_horvu[, grep("mla6|bln1|dm", colnames(padj_table_horvu))]
padj_table_horvu<- padj_table_horvu[, grep("t0|t16|t20|t32", colnames(padj_table_horvu))]
padj_table_horvu<- apply(padj_table_horvu, 2, FUN=function(x)ifelse(x>0.001 | is.na(x), 0,1))
#padj_table_horvu<- padj_table_horvu[rowSums(padj_table_horvu)>0,]
ncol(padj_table_horvu)
nrow(padj_table_horvu)

table_promoters<- read.csv("~/iowa_state/lab/Chromatin_accessibility_paper/epistatic_inputs/table_promoters_epistasis.csv", row.names = 1)

table_any_regions<- read.csv("~/iowa_state/lab/Chromatin_accessibility_paper/epistatic_inputs/table_any_regions_epistasis.csv")
table_any_regions<- table_any_regions[!duplicated(table_any_regions$gene),]
rownames(table_any_regions)<- table_any_regions$gene
table_any_regions$gene<-NULL


table_promoters <- table_any_regions


table_promoters<- table_promoters[grep("HORVU", rownames(table_promoters)), grep("wt", colnames(table_promoters))]
table_promoters<- table_promoters[, grep("padj", colnames(table_promoters))]

table_promoters<- apply(table_promoters, 2, FUN=function(x)ifelse(x>0.005 | is.na(x), 0,1))

table_promoters<- apply(table_promoters, 2, FUN=function(x)ifelse(x>0.001 | is.na(x), 0,1))


ncol(table_promoters)
nrow(table_promoters)


sum(!rownames(table_promoters) %in% rownames(padj_table_horvu))

summary<- data.frame(time=rep(paste0("t", time), 3),
                     reference= "wt",
                     genotype=rep(genotypes, each=4),
                     DE_acc=0, DE_nonAcc=0, NonDE_Acc=0, NonDE_nonAcc=0)

for(sample in  colnames(padj_table_horvu)){
  info<- strsplit(sample, "_")[[1]]
  t<- info[1]
  gen<- info[3]
  DE_acc=sum(rownames(padj_table_horvu[padj_table_horvu[,sample]==1,]) %in% rownames(table_promoters[table_promoters[,sample]==1,]))
  #DE_nonAcc=sum(rownames(padj_table_horvu[padj_table_horvu[,sample]==1,]) %in% rownames(table_promoters[table_promoters[,sample]!=1,]))
  DE_nonAcc=sum(padj_table_horvu[,sample]==1) - DE_acc
  
  NonDE_Acc=sum(rownames(padj_table_horvu[padj_table_horvu[,sample]!=1,]) %in% rownames(table_promoters[table_promoters[,sample]==1,]))
  
  #NonDE_nonAcc=sum(rownames(padj_table_horvu[padj_table_horvu[,sample]!=1,]) %in% rownames(table_promoters[table_promoters[,sample]!=1,]))
  NonDE_nonAcc=sum(padj_table_horvu[,sample]!=1) - NonDE_Acc
  
  summary[summary$time==t &summary$genotype==gen, 4:7]<- c(DE_acc, DE_nonAcc, NonDE_Acc, NonDE_nonAcc)
}

summary_001<- summary
summary_005<- summary
summary_005_any_region<- summary

write.csv(summary_001, "~/iowa_state/lab/Chromatin_accessibility_paper/ATACSeq/correlation_RNASeq/summary_001.csv")
write.csv(summary_005, "~/iowa_state/lab/Chromatin_accessibility_paper/ATACSeq/correlation_RNASeq/summary_005.csv")
write.csv(summary_005_any_region, "~/iowa_state/lab/Chromatin_accessibility_paper/ATACSeq/correlation_RNASeq/summary_005_any_region.csv")

summary_005_any_region<- read.csv("~/iowa_state/lab/Chromatin_accessibility_paper/ATACSeq/correlation_RNASeq/summary_005_any_region.csv", row.names = 1)

library(reshape2)
summary_005_any_region_long<- gather(summary_005_any_region, key = "Type", value = "Num_genes", 4:7)

pdf(paste0("~/iowa_state/lab/Chromatin_accessibility_paper/overlap_ATAC_RNASeq.pdf"), width = 5, height = 5, fonts = "ArialMT", pointsize = 30)

ggplot(summary_005_any_region_long, aes(x=Type, y=Num_genes, fill = Type)) +
  geom_bar(stat="identity") +
  facet_grid(~time~genotype) +
  scale_y_log10() + 
  theme(axis.text.x = element_text(angle = 90)) 

ggplot(summary_005_any_region_long[summary_005_any_region_long$Type!="NonDE_nonAcc",], aes(x=Type, y=Num_genes, fill = Type)) +
  geom_bar(stat="identity") +
  facet_grid(~time~genotype) +
  #scale_y_log10() + 
  theme(axis.text.x = element_text(angle = 90)) 

dev.off()


#distribution of regions in the atacseq data pie chart, need to get the metrics
# Load necessary library
library(tidyverse)

table_all_regions_meta_padj005<- read.csv("~/iowa_state/lab/Chromatin_accessibility_paper/ATACSeq/stats/table_all_regions_meta_padj005.csv")
table_all_regions_meta_padj005[is.na(table_all_regions_meta_padj005)]<-0

data <- data.frame(
  Category = table_all_regions_meta_padj005$Annotation,
  Count = rowSums(table_all_regions_meta_padj005[,-1])
)

# Compute percentages
data$Percentage <- round(data$Count / sum(data$Count) * 100, 1)
data$ypos <- cumsum(data$Count) - data$Count / 2 # Adjust for better placement


# Create pie chart with percentage labels
pdf(paste0("~/iowa_state/lab/Chromatin_accessibility_paper/piechart_peak_distribution.pdf"), width = 15, height = 5, fonts = "ArialMT", pointsize = 30)
ggplot(data, aes(x = "", y = Count, fill = Category)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar("y", start = 0) +
  theme_void() +
  geom_text(aes(x = 1.7, y = ypos, label = paste0(Percentage, "%")), 
            size = 5) +  # Adjust x value to move text outward
  labs(fill = "Category", title = "Pie Chart with Percentages Outside") +
  theme(plot.title = element_text(hjust = 0.5))


ggplot(data, aes(x = "", y = Count, fill = Category)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar("y", start = 0) +
  theme_void() +
  #geom_text(aes(x = 1.7, y = ypos, label = paste0(Percentage, "%")), size = 5) +  # Adjust x value to move text outward
  labs(fill = "Category", title = "Pie Chart with Percentages Outside") +
  theme(plot.title = element_text(hjust = 0.5))

dev.off()


#distribution of distance to TSS
library(tidyverse)
library(reshape2)
table_distance<- read.csv("~/iowa_state/lab/Chromatin_accessibility_paper/epistatic_inputs/table_promoters_distanceTSS.csv")

table_distance1<- table_distance[, c(1,grep("TSS", colnames(table_distance)))]
table_distance_long1<- gather(table_distance1, key = "comparison", value = "TSS_distance", 2:25)

table_distance2<- table_distance[, c(1,grep("padj", colnames(table_distance)))]
table_distance_long2<- gather(table_distance2, key = "comparison", value = "padj", 2:25)

table_distance_long1$padj<- table_distance_long2$padj
table_distance_long1<- table_distance_long1[!is.na(table_distance_long1$padj),]

pdf(paste0("~/iowa_state/lab/Chromatin_accessibility_paper/TSS_distance_distribution.pdf"), width = 5, height = 5, fonts = "ArialMT", pointsize = 30)

ggplot(table_distance_long1[table_distance_long1$padj<0.05 &
                              abs(table_distance_long1$TSS_distance)<5000,], aes(x = TSS_distance)) +
  geom_density(fill = "lightblue", alpha = 0.6) +
  theme_minimal(base_size = 14) 

dev.off()

write.csv(table_distance_long1, "~/iowa_state/lab/Chromatin_accessibility_paper/TSS_distance_table_long.csv")

#repeat dE vs accessible but only with the promoters in the 5k range

#plot bars with DEG accessible and non accessible

library(tidyverse)
#Now do number of DE genes per timepoint 
time=c(0,16, 20, 32)
genotypes<-c("mla6", "bln1", "dm")

#padj_table<- read.csv("~/iowa_state/lab/RNAseq/RNASeq_paper/vvz manuscript/analysis/hv_R2_genes_deseq2_results_pairwise_genotype_padj_tax_sp.csv", row.names = 1)
padj_table<- read.csv("~/iowa_state/lab/RNAseq/DESeq2 HvR3/hv_R3_genes_deseq2_results_pairwise_genotype_padj_tax_sp.csv", row.names = 1)
padj_table_horvu<- padj_table[grep("HORVU", rownames(padj_table)), grep("wt", colnames(padj_table))]
padj_table_horvu<- padj_table_horvu[, grep("mla6|bln1|dm", colnames(padj_table_horvu))]
padj_table_horvu<- padj_table_horvu[, grep("t0|t16|t20|t32", colnames(padj_table_horvu))]
padj_table_horvu<- apply(padj_table_horvu, 2, FUN=function(x)ifelse(x>0.001 | is.na(x), 0,1))
#padj_table_horvu<- padj_table_horvu[rowSums(padj_table_horvu)>0,]
ncol(padj_table_horvu)
nrow(padj_table_horvu)

table_distance_long<- read.csv("~/iowa_state/lab/Chromatin_accessibility_paper/TSS_distance_table_long.csv", row.names = 1)
table_distance_long<- table_distance_long[grep("wt", table_distance_long$comparison),]
table_distance_long<- table_distance_long[abs(table_distance_long$TSS_distance)<5000 & table_distance_long$padj<0.005,]
table_distance_long$comparison<- gsub("_distanceTSS", "_padj", table_distance_long$comparison)


sum(colnames(padj_table_horvu) %in% unique(table_distance_long$comparison))

sum(!rownames(table_promoters) %in% rownames(padj_table_horvu))

summary<- data.frame(time=rep(paste0("t", time), 3),
                     reference= "wt",
                     genotype=rep(genotypes, each=4),
                     DE_acc=0, DE_nonAcc=0, NonDE_Acc=0, NonDE_nonAcc=0)

for(sample in  colnames(padj_table_horvu)){
  info<- strsplit(sample, "_")[[1]]
  t<- info[1]
  gen<- info[3]
  DE_acc=sum(rownames(padj_table_horvu[padj_table_horvu[,sample]==1,]) %in% table_distance_long[table_distance_long$comparison==sample,"gene"])
  
  #DE_acc=sum(rownames(padj_table_horvu[padj_table_horvu[,sample]==1,]) %in% rownames(table_promoters[table_promoters[,sample]==1,]))
  #DE_nonAcc=sum(rownames(padj_table_horvu[padj_table_horvu[,sample]==1,]) %in% rownames(table_promoters[table_promoters[,sample]!=1,]))
  DE_nonAcc=sum(padj_table_horvu[,sample]==1) - DE_acc
  
  NonDE_Acc=sum(rownames(padj_table_horvu[padj_table_horvu[,sample]!=1,]) %in% table_distance_long[table_distance_long$comparison==sample,"gene"])
  
  #NonDE_nonAcc=sum(rownames(padj_table_horvu[padj_table_horvu[,sample]!=1,]) %in% rownames(table_promoters[table_promoters[,sample]!=1,]))
  NonDE_nonAcc=sum(padj_table_horvu[,sample]!=1) - NonDE_Acc
  
  summary[summary$time==t &summary$genotype==gen, 4:7]<- c(DE_acc, DE_nonAcc, NonDE_Acc, NonDE_nonAcc)
}


write.csv(summary, "~/iowa_state/lab/Chromatin_accessibility_paper/ATACSeq/correlation_RNASeq/summary_005_propoters_5000bp_TSS.csv")

library(reshape2)
summary_005_any_region_long<- gather(summary, key = "Type", value = "Num_genes", 4:7)

pdf(paste0("~/iowa_state/lab/Chromatin_accessibility_paper/overlap_ATAC_RNASeq_TSS_5000bp.pdf"), width = 5, height = 5, fonts = "ArialMT", pointsize = 30)

ggplot(summary_005_any_region_long, aes(x=Type, y=Num_genes, fill = Type)) +
  geom_bar(stat="identity") +
  facet_grid(~time~genotype) +
  scale_y_log10() + 
  theme(axis.text.x = element_text(angle = 90)) 

ggplot(summary_005_any_region_long[summary_005_any_region_long$Type!="NonDE_nonAcc",], aes(x=Type, y=Num_genes, fill = Type)) +
  geom_bar(stat="identity") +
  facet_grid(~time~genotype) +
  #scale_y_log10() + 
  theme(axis.text.x = element_text(angle = 90)) 

dev.off()




#HG test in bind for de genes and accesible genes

#HG test DE genes

HG_test<- function(trts, bin_table, DE_table_trts, overlap_table, padj=0.001){
  #table trts by bins
  hypergeom_test <- data.frame(matrix(ncol = 1+2*length(trts), nrow = nrow(bin_table))) 
  hypergeom_test[,1]<- bin_table[,1]
  colnames(hypergeom_test)<- c("binID", unlist(sapply(trts, function(x) c(paste0(x,"_pvalue"), paste0(x,"_padj")))))
  
  for(m in 1:length(trts)){
    DE_table<- rownames(DE_table_trts[DE_table_trts[,m]< padj & !is.na(DE_table_trts[,m]),])
    if(length(DE_table)>0){
      overlap_table_filtered<- overlap_table[overlap_table$seqname %in% DE_table, c("seqname", colnames(bin_table)[1])]
      DE_genes_bins<-data.frame(table(overlap_table_filtered[,colnames(bin_table)[1]]))
      overlap_table_filtered<- merge(overlap_table_filtered, DE_genes_bins, by.x=colnames(bin_table)[1], by.y="Var1", all.x=T)
      DE_bins<- bin_table[unlist(bin_table[,1]) %in% DE_genes_bins$Var1,]
      
      #total DE in bin, total DE trt, sum not DE in bins listed in trt, genes in bin
      hypergeom_test[match(DE_genes_bins$Var1, unlist(hypergeom_test$binID)), 
                     paste0(trts[m],"_pvalue")]<- phyper(DE_genes_bins$Freq, length(DE_table), 
                                                         sum(DE_bins$Freq - DE_genes_bins$Freq), 
                                                         DE_bins$Freq, lower.tail = FALSE)
      
      hypergeom_test[unlist(hypergeom_test$binID) %in% DE_genes_bins$Var1,paste0(trts[m],"_padj")]<- p.adjust(hypergeom_test[unlist(hypergeom_test$binID) %in% DE_genes_bins$Var1,paste0(trts[m],"_pvalue")], method = "BH")
      print(m)
    }
    
    
  }
  hypergeom_test<- hypergeom_test[rowSums(is.na(hypergeom_test[,-1])) != ncol(hypergeom_test[,-1]), ]
  return(hypergeom_test)
}

#test for all bin tables
padj_table<- read.csv("~/iowa_state/lab/RNAseq/DESeq2 HvR3/hv_R3_genes_deseq2_results_pairwise_genotype_padj_tax_sp.csv", row.names = 1)
padj_table<- padj_table[grep("HORVU", rownames(padj_table)),]
DE_table_trts<- padj_table

trts<- colnames(padj_table)

folder<- "~/iowa_state/lab/RNAseq/RNASeq_paper/vvz manuscript/analysis/chr analysis/R3/"

complete_overlap_df<- read.csv(paste0(folder, "Binned_genes_HR3.csv"))

overlap_table<- complete_overlap_df

bins_1Mb<- read.csv( paste0(folder, "bins_1Mb.csv"))
bins_10Mb<-read.csv(paste0(folder, "bins_10Mb.csv"))
bins_100Mb<-read.csv( paste0(folder, "bins_100Mb.csv"))
bins_chr<-read.csv( paste0(folder, "bins_chrs.csv"))

folder_out<- "~/iowa_state/lab/Chromatin_accessibility_paper/ATACSeq/correlation_RNASeq/DE_HG_tests_RNASeq/"
#bin_table, DE_table_trts, overlap_table

HG_1Mb<- HG_test(trts, bin_table=bins_1Mb, DE_table_trts=padj_table, overlap_table, padj=0.001)
write.csv(HG_1Mb, paste0(folder_out, "hypergeometric_test_bins_1Mb.csv"), row.names = F)

HG_10Mb<- HG_test(trts, bins_10Mb, padj_table, overlap_table, padj=0.001)
write.csv(HG_10Mb, paste0(folder_out, "hypergeometric_test_bins_10Mb.csv"))

HG_100Mb<- HG_test(trts, bins_100Mb, padj_table, overlap_table, padj=0.001)
write.csv(HG_100Mb, paste0(folder_out, "hypergeometric_test_bins_100Mb.csv"))

HG_chrs<- HG_test(trts, bins_chr, padj_table, overlap_table, padj=0.001)
write.csv(HG_chrs, paste0(folder_out, "hypergeometric_test_bins_chrs.csv"))


#differentially accessible 

table_promoters<- read.csv("~/iowa_state/lab/Chromatin_accessibility_paper/epistatic_inputs/table_promoters_epistasis.csv", row.names = 1)

table_any_regions<- read.csv("~/iowa_state/lab/Chromatin_accessibility_paper/epistatic_inputs/table_any_regions_epistasis.csv")
table_any_regions<- table_any_regions[!duplicated(table_any_regions$gene),]
rownames(table_any_regions)<- table_any_regions$gene
table_any_regions$gene<-NULL


padj_table<- table_promoters
padj_table<- table_any_regions

padj_table<- padj_table[grep("HORVU", rownames(padj_table)),grep("padj", colnames(padj_table))]
DE_table_trts<- padj_table

trts<- colnames(padj_table)

folder<- "~/iowa_state/lab/RNAseq/RNASeq_paper/vvz manuscript/analysis/chr analysis/R3/"

overlap_table<- read.csv(paste0(folder, "Binned_genes_HR3.csv"))

bins_1Mb<- read.csv( paste0(folder, "bins_1Mb.csv"))
bins_10Mb<-read.csv(paste0(folder, "bins_10Mb.csv"))
bins_100Mb<-read.csv( paste0(folder, "bins_100Mb.csv"))
bins_chr<-read.csv( paste0(folder, "bins_chrs.csv"))


folder_out<- "~/iowa_state/lab/Chromatin_accessibility_paper/ATACSeq/correlation_RNASeq/DE_HG_tests_ATACSeq/promoters_005/"
folder_out<- "~/iowa_state/lab/Chromatin_accessibility_paper/ATACSeq/correlation_RNASeq/DE_HG_tests_ATACSeq/Any_region_005/"

HG_1Mb_atac<- HG_test(trts, bin_table=bins_1Mb, DE_table_trts=padj_table, overlap_table, padj=0.05)
write.csv(HG_1Mb_atac, paste0(folder_out, "hypergeometric_test_bins_1Mb.csv"), row.names = F)

HG_10Mb_atac<- HG_test(trts, bins_10Mb, padj_table, overlap_table, padj=0.05)
write.csv(HG_10Mb_atac, paste0(folder_out, "hypergeometric_test_bins_10Mb.csv"))

HG_100Mb_atac<- HG_test(trts, bins_100Mb, padj_table, overlap_table, padj=0.05)
write.csv(HG_100Mb_atac, paste0(folder_out, "hypergeometric_test_bins_100Mb.csv"))

HG_chrs_atac<- HG_test(trts, bins_chr, padj_table, overlap_table, padj=0.05)
write.csv(HG_chrs_atac, paste0(folder_out, "hypergeometric_test_bins_chrs.csv"))


#Take DE and DA results and make matrices to visualize the correlation, by trt
#do it just ofr wt vs mla6
library(reshape2)
library(tidyverse)
HG_100Mb_atac<- read.csv("~/iowa_state/lab/Chromatin_accessibility_paper/ATACSeq/correlation_RNASeq/DE_HG_tests_ATACSeq/Any_region_005/hypergeometric_test_bins_100Mb.csv", row.names = 1)
HG_100Mb<- read.csv("~/iowa_state/lab/Chromatin_accessibility_paper/ATACSeq/correlation_RNASeq/DE_HG_tests_RNASeq/hypergeometric_test_bins_100Mb.csv", row.names = 1)

colnames(HG_100Mb) %in% colnames(HG_100Mb_atac)

colnames(HG_100Mb_atac)[!colnames(HG_100Mb_atac) %in% colnames(HG_100Mb)]


HG_100Mb_filtered<- HG_100Mb[,grep("wt_mla6|binID", colnames(HG_100Mb))]
HG_100Mb_filtered<- HG_100Mb_filtered[,grep("t0|t16|t20|t32|binID", colnames(HG_100Mb_filtered))]
HG_100Mb_filtered<- HG_100Mb_filtered[,grep("padj_padj|binID", colnames(HG_100Mb_filtered))]
HG_100Mb_filtered[is.na(HG_100Mb_filtered)]<-1


HG_100Mb_atac_filtered<- HG_100Mb_atac[,grep("wt_mla6|binID", colnames(HG_100Mb_atac))]
HG_100Mb_atac_filtered<- HG_100Mb_atac_filtered[,grep("padj_padj|binID", colnames(HG_100Mb_atac_filtered))]
HG_100Mb_atac_filtered[is.na(HG_100Mb_atac_filtered)]<-1

vec1<- as.numeric(as.vector(as.matrix(HG_100Mb_filtered[,-1])))
vec2<- as.numeric(as.vector(as.matrix(HG_100Mb_atac_filtered[,-1])))

cor<- cor(vec1, vec2, method = "pearson")
cor_df<- data.frame(DE=vec1, DA=vec2)

ggplot(cor_df, aes(x=log(vec1), y=log(vec2))) +
  geom_point()

plot(vec1, vec2)

#make heatmap with bins and significant non significant? thr 0.05
#one per comp and an average 0,1

# Example binary matrices
mat1 <- HG_100Mb_filtered[,2:5]
rownames(mat1)<- HG_100Mb_filtered$binID
mat1<- ifelse(mat1<0.05, 1,0)

mat2 <- HG_100Mb_atac_filtered[,2:5]
rownames(mat2)<- HG_100Mb_atac_filtered$binID
mat2<- ifelse(mat2<0.05, 1,0)
# Generate agreement matrix (1 if equal, 0 if different)
agreement_matrix <- ifelse(mat1 == mat2, 1, 0)

# Visualize with a heatmap
library(ggplot2)
library(reshape2)

# Convert matrix to long format for ggplot
df <- melt(agreement_matrix)

# Heatmap plot
ggplot(df, aes(Var2, Var1, fill=value)) +
  geom_tile() +
  scale_fill_gradient(low="red", high="blue") + # Corrected scale
  theme_minimal() +
  labs(title="Binary Agreement Heatmap", fill="Agreement")


# Recode agreement matrix
recoded_matrix <- ifelse(mat1 == 1 & mat2 == 1, "DE-DA", 
                         ifelse(mat1 == 0 & mat2 == 0, "NotDE-NotDA", 
                                ifelse(mat1 == 0 & mat2 == 1, "NotDE-DA", 
                                       "DE-NotDA")))

# Visualize heatmap
library(ggplot2)
library(reshape2)
library(viridis)

df <- melt(recoded_matrix)

ggplot(df, aes(Var2, Var1, fill=value)) +
  geom_tile() +
  scale_fill_manual(values=c("DE-DA"="darkgreen", "NotDE-NotDA"="gray90", "NotDE-DA"="darkorange", "DE-NotDA"="turquoise")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(title="DE vs DA genes", fill="Category")


ggplot(df, aes(Var2, Var1, fill=value)) +
  geom_tile() +
  scale_fill_manual(values=c("DE-DA"="darkgreen", "NotDE-NotDA"="gray90", "NotDE-DA"="darkorange", "DE-NotDA"="turquoise")) +
  scale_x_discrete(labels=c("t0_wt_mla6_padj_padj"="0",  "t16_wt_mla6_padj_padj"="16", "t20_wt_mla6_padj_padj"="20", "t32_wt_mla6_padj_padj"="32")) +  # Custom labels
  labs(title="DE vs DA genes", fill="Category", y="Location", x="Hours after inoculation")

#"t0_wt_mla6_padj_padj"="0",  "t16_wt_mla6_padj_padj"="16", "t20_wt_mla6_padj_padj"="20", "t32_wt_mla6_padj_padj"="32"

df$chr<- substr(df$Var1, 1,2)
df$bin<- substr(df$Var1, 4,4)

pdf(paste0("~/iowa_state/lab/Chromatin_accessibility_paper/DE_vs_DA_bins.pdf"), width = 5, height = 8, fonts = "ArialMT", pointsize = 30)

ggplot(df, aes(Var2, Var1, fill=value)) +
  geom_tile() +
  scale_fill_manual(values=c("DE-DA"="darkgreen", "NotDE-NotDA"="gray90", "NotDE-DA"="darkorange", "DE-NotDA"="turquoise")) +
  scale_x_discrete(labels=c("t0_wt_mla6_padj_padj"="0",  "t16_wt_mla6_padj_padj"="16", "t20_wt_mla6_padj_padj"="20", "t32_wt_mla6_padj_padj"="32")) +  # Custom labels
  labs(title="DE vs DA genes", fill="Category", y="Location", x="Hours after inoculation")

ggplot(df, aes(Var2, Var1, fill=value)) +
  geom_tile() +
  scale_y_discrete(limits=rev) +
  scale_fill_manual(values=c("DE-DA"="darkgreen", "NotDE-NotDA"="gray90", "NotDE-DA"="darkorange", "DE-NotDA"="turquoise")) +
  scale_x_discrete(labels=c("t0_wt_mla6_padj_padj"="0",  "t16_wt_mla6_padj_padj"="16", "t20_wt_mla6_padj_padj"="20", "t32_wt_mla6_padj_padj"="32")) +  # Custom labels
  labs(title="DE vs DA genes", fill="Category", y="Location", x="Hours after inoculation")


ggplot(df, aes(Var2, bin, fill=value)) +
  geom_tile() +
  scale_y_discrete(limits=rev) +
  facet_grid(~chr) +
  scale_fill_manual(values=c("DE-DA"="darkgreen", "NotDE-NotDA"="gray90", "NotDE-DA"="darkorange", "DE-NotDA"="turquoise")) +
  scale_x_discrete(labels=c("t0_wt_mla6_padj_padj"="0",  "t16_wt_mla6_padj_padj"="16", "t20_wt_mla6_padj_padj"="20", "t32_wt_mla6_padj_padj"="32")) +  # Custom labels
  labs(title="DE vs DA genes", fill="Category", y="Bin Location", x="Hours after inoculation")


dev.off()


pdf(paste0("~/iowa_state/lab/Chromatin_accessibility_paper/DE_vs_DA_bins_facet.pdf"), width = 8, height = 5, fonts = "ArialMT", pointsize = 30)

ggplot(df, aes(Var2, bin, fill=value)) +
  geom_tile() +
  scale_y_discrete(limits=rev) +
  facet_grid(~chr) +
  scale_fill_manual(values=c("DE-DA"="darkgreen", "NotDE-NotDA"="gray90", "NotDE-DA"="darkorange", "DE-NotDA"="turquoise")) +
  scale_x_discrete(labels=c("t0_wt_mla6_padj_padj"="0",  "t16_wt_mla6_padj_padj"="16", "t20_wt_mla6_padj_padj"="20", "t32_wt_mla6_padj_padj"="32")) +  # Custom labels
  labs(title="DE vs DA genes", fill="Category", y="Bin Location (100Mb)", x="Hours after inoculation")

ggplot(df, aes(Var2, bin, fill=value)) +
  geom_tile() +
  scale_y_discrete(limits=rev) +
  facet_wrap(~chr,ncol = 1) +
  scale_fill_manual(values=c("DE-DA"="darkgreen", "NotDE-NotDA"="gray90", "NotDE-DA"="darkorange", "DE-NotDA"="turquoise")) +
  scale_x_discrete(labels=c("t0_wt_mla6_padj_padj"="0",  "t16_wt_mla6_padj_padj"="16", "t20_wt_mla6_padj_padj"="20", "t32_wt_mla6_padj_padj"="32")) +  # Custom labels
  labs(title="DE vs DA genes", fill="Category", y="Bin Location (100Mb)", x="Hours after inoculation")


dev.off()



ggplot(df, aes(Var2, Var1, fill=value)) +
  geom_tile() +
  scale_fill_viridis_d(option="turbo") +
  #scale_fill_manual(values=c("DE-DA"="green", "NotDE-NotDA"="white", "NotDE-DA"="yellow", "DE-NotDA"="turquoise")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(title="Recoded Binary Agreement Heatmap", fill="Category")






#try with 10 mb
library(tidyverse)
HG_10Mb_atac<- read.csv("~/iowa_state/lab/Chromatin_accessibility_paper/ATACSeq/correlation_RNASeq/DE_HG_tests_ATACSeq/Any_region_005/hypergeometric_test_bins_10Mb.csv", row.names = 1)
HG_10Mb<- read.csv("~/iowa_state/lab/Chromatin_accessibility_paper/ATACSeq/correlation_RNASeq/DE_HG_tests_RNASeq/hypergeometric_test_bins_10Mb.csv", row.names = 1)

colnames(HG_10Mb) %in% colnames(HG_10Mb_atac)

colnames(HG_10Mb_atac)[!colnames(HG_10Mb_atac) %in% colnames(HG_10Mb)]


HG_10Mb_filtered<- HG_10Mb[,grep("wt_mla6|binID", colnames(HG_10Mb))]
HG_10Mb_filtered<- HG_10Mb_filtered[,grep("t0|t16|t20|t32|binID", colnames(HG_10Mb_filtered))]
HG_10Mb_filtered<- HG_10Mb_filtered[,grep("padj_padj|binID", colnames(HG_10Mb_filtered))]
HG_10Mb_filtered[is.na(HG_10Mb_filtered)]<-1


HG_10Mb_atac_filtered<- HG_10Mb_atac[,grep("wt_mla6|binID", colnames(HG_10Mb_atac))]
HG_10Mb_atac_filtered<- HG_10Mb_atac_filtered[,grep("padj_padj|binID", colnames(HG_10Mb_atac_filtered))]
HG_10Mb_atac_filtered[is.na(HG_10Mb_atac_filtered)]<-1

#add missing rows from bins present in the other table

HG_10Mb_filtered$binID[!HG_10Mb_filtered$binID %in% HG_10Mb_atac_filtered$binID]

data.frame(binID=HG_10Mb_filtered$binID[!HG_10Mb_filtered$binID %in% HG_10Mb_atac_filtered$binID],
           t0_wt_mla6_padj_padj=1, t16_wt_mla6_padj_padj=1, t20_wt_mla6_padj_padj=1,
           t32_wt_mla6_padj_padj=1)

HG_10Mb_atac_filtered<- rbind(HG_10Mb_atac_filtered, data.frame(binID=HG_10Mb_filtered$binID[!HG_10Mb_filtered$binID %in% HG_10Mb_atac_filtered$binID],
                                                                t0_wt_mla6_padj_padj=1, t16_wt_mla6_padj_padj=1, t20_wt_mla6_padj_padj=1,
                                                                t32_wt_mla6_padj_padj=1))

HG_10Mb_atac_filtered$binID[!HG_10Mb_atac_filtered$binID %in% HG_10Mb_filtered$binID]


HG_10Mb_filtered<- rbind(HG_10Mb_filtered, data.frame(binID=HG_10Mb_atac_filtered$binID[!HG_10Mb_atac_filtered$binID %in% HG_10Mb_filtered$binID],
                                                                t0_wt_mla6_padj_padj=1, t16_wt_mla6_padj_padj=1, t20_wt_mla6_padj_padj=1,
                                                                t32_wt_mla6_padj_padj=1))

HG_10Mb_filtered<- HG_10Mb_filtered[order(HG_10Mb_filtered$binID),]

HG_10Mb_atac_filtered<- HG_10Mb_atac_filtered[order(HG_10Mb_atac_filtered$binID),]

vec1<- as.numeric(as.vector(as.matrix(HG_10Mb_filtered[,-1])))
vec2<- as.numeric(as.vector(as.matrix(HG_10Mb_atac_filtered[,-1])))

cor<- cor(vec1, vec2, method = "pearson")
cor_df<- data.frame(DE=vec1, DA=vec2)

ggplot(cor_df, aes(x=log(vec1), y=log(vec2))) +
  geom_point()

plot(vec1, vec2)

#make heatmap with bins and significant non significant? thr 0.05
#one per comp and an average 0,1

# Example binary matrices
mat1 <- HG_10Mb_filtered[,2:5]
rownames(mat1)<- HG_10Mb_filtered$binID
mat1<- ifelse(mat1<0.05, 1,0)

mat2 <- HG_10Mb_atac_filtered[,2:5]
rownames(mat2)<- HG_10Mb_atac_filtered$binID
mat2<- ifelse(mat2<0.05, 1,0)
# Generate agreement matrix (1 if equal, 0 if different)
agreement_matrix <- ifelse(mat1 == mat2, 1, 0)

# Visualize with a heatmap
library(ggplot2)
library(reshape2)

# Convert matrix to long format for ggplot
df <- melt(agreement_matrix)

# Heatmap plot
ggplot(df, aes(Var2, Var1, fill=value)) +
  geom_tile() +
  scale_fill_gradient(low="red", high="blue") + # Corrected scale
  theme_minimal() +
  labs(title="Binary Agreement Heatmap", fill="Agreement")


# Recode agreement matrix
recoded_matrix <- ifelse(mat1 == 1 & mat2 == 1, "DE-DA", 
                         ifelse(mat1 == 0 & mat2 == 0, "NotDE-NotDA", 
                                ifelse(mat1 == 0 & mat2 == 1, "NotDE-DA", 
                                       "DE-NotDA")))

# Visualize heatmap
library(ggplot2)
library(reshape2)
library(viridis)

df <- melt(recoded_matrix)
df$chr<- substr(df$Var1, 1,2)
df$bin<- substr(df$Var1, 4,6)

ggplot(df, aes(Var2, Var1, fill=value)) +
  geom_tile() +
  scale_fill_manual(values=c("DE-DA"="darkgreen", "NotDE-NotDA"="gray90", "NotDE-DA"="darkorange", "DE-NotDA"="turquoise")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(title="DE vs DA genes", fill="Category")


ggplot(df, aes(Var2, Var1, fill=value)) +
  geom_tile() +
  scale_fill_manual(values=c("DE-DA"="darkgreen", "NotDE-NotDA"="gray90", "NotDE-DA"="darkorange", "DE-NotDA"="turquoise")) +
  scale_x_discrete(labels=c("t0_wt_mla6_padj_padj"="0",  "t16_wt_mla6_padj_padj"="16", "t20_wt_mla6_padj_padj"="20", "t32_wt_mla6_padj_padj"="32")) +  # Custom labels
  labs(title="DE vs DA genes", fill="Category", y="Location", x="Hours after inoculation")

#"t0_wt_mla6_padj_padj"="0",  "t16_wt_mla6_padj_padj"="16", "t20_wt_mla6_padj_padj"="20", "t32_wt_mla6_padj_padj"="32"


pdf(paste0("~/iowa_state/lab/Chromatin_accessibility_paper/DE_vs_DA_bins_10Mb.pdf"), width = 5, height = 12, fonts = "ArialMT", pointsize = 30)

ggplot(df, aes(Var2, Var1, fill=value)) +
  geom_tile() +
  scale_fill_manual(values=c("DE-DA"="darkgreen", "NotDE-NotDA"="gray90", "NotDE-DA"="darkorange", "DE-NotDA"="turquoise")) +
  scale_x_discrete(labels=c("t0_wt_mla6_padj_padj"="0",  "t16_wt_mla6_padj_padj"="16", "t20_wt_mla6_padj_padj"="20", "t32_wt_mla6_padj_padj"="32")) +  # Custom labels
  labs(title="DE vs DA genes", fill="Category", y="Location", x="Hours after inoculation")

ggplot(df, aes(Var2, Var1, fill=value)) +
  geom_tile() +
  scale_y_discrete(limits=rev) +
  scale_fill_manual(values=c("DE-DA"="darkgreen", "NotDE-NotDA"="gray90", "NotDE-DA"="darkorange", "DE-NotDA"="turquoise")) +
  scale_x_discrete(labels=c("t0_wt_mla6_padj_padj"="0",  "t16_wt_mla6_padj_padj"="16", "t20_wt_mla6_padj_padj"="20", "t32_wt_mla6_padj_padj"="32")) +  # Custom labels
  labs(title="DE vs DA genes", fill="Category", y="Location", x="Hours after inoculation")


ggplot(df, aes(Var2, bin, fill=value)) +
  geom_tile() +
  scale_y_discrete(limits=rev) +
  facet_grid(~chr) +
  scale_fill_manual(values=c("DE-DA"="darkgreen", "NotDE-NotDA"="gray90", "NotDE-DA"="darkorange", "DE-NotDA"="turquoise")) +
  scale_x_discrete(labels=c("t0_wt_mla6_padj_padj"="0",  "t16_wt_mla6_padj_padj"="16", "t20_wt_mla6_padj_padj"="20", "t32_wt_mla6_padj_padj"="32")) +  # Custom labels
  labs(title="DE vs DA genes", fill="Category", y="Bin Location", x="Hours after inoculation")

dev.off()


pdf(paste0("~/iowa_state/lab/Chromatin_accessibility_paper/DE_vs_DA_bins_facet_10Mb.pdf"), width = 8, height = 10, fonts = "ArialMT", pointsize = 30)

ggplot(df, aes(Var2, bin, fill=value)) +
  geom_tile() +
  scale_y_discrete(limits=rev) +
  facet_grid(~chr) +
  scale_fill_manual(values=c("DE-DA"="darkgreen", "NotDE-NotDA"="gray90", "NotDE-DA"="darkorange", "DE-NotDA"="turquoise")) +
  scale_x_discrete(labels=c("t0_wt_mla6_padj_padj"="0",  "t16_wt_mla6_padj_padj"="16", "t20_wt_mla6_padj_padj"="20", "t32_wt_mla6_padj_padj"="32")) +  # Custom labels
  labs(title="DE vs DA genes", fill="Category", y="Bin Location (10Mb)", x="Hours after inoculation")

ggplot(df, aes(Var2, bin, fill=value)) +
  geom_tile() +
  scale_y_discrete(limits=rev) +
  facet_wrap(~chr,ncol = 1) +
  scale_fill_manual(values=c("DE-DA"="darkgreen", "NotDE-NotDA"="gray90", "NotDE-DA"="darkorange", "DE-NotDA"="turquoise")) +
  scale_x_discrete(labels=c("t0_wt_mla6_padj_padj"="0",  "t16_wt_mla6_padj_padj"="16", "t20_wt_mla6_padj_padj"="20", "t32_wt_mla6_padj_padj"="32")) +  # Custom labels
  labs(title="DE vs DA genes", fill="Category", y="Bin Location (10Mb)", x="Hours after inoculation")


dev.off()





#characterize DA genes wt vs mla6


