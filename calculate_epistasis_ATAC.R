#define main operations to calculate epistasis


#start  y defining the colnames expected by epistasis funcitons
FC_padj_table<- read.csv("~/iowa_state/lab/RNAseq/DESeq2 HvR3/hv_R3_genes_deseq2_results_pairwise_genotype_logfc_padj_tax_sp.csv")
colnames(FC_padj_table)
#gene, t_g1_g2_padj, t_g1_g2_log2FoldChange

#atacseq folders: diff.g1t1_vs_g2t2, in this case t1=t2
#file names: diff.g1t1_vs_g2t2.results_allpeaks_annotation

#extract gene, region, coordinates, padj, log2fc
#padj Na =1


# Load required libraries
library(dplyr)

# Define the path to the parent directory containing folders
parent_dir <- "~/iowa_state/lab/ATACSeq/Analysis results feb 2025/Differential_analysis/Differential_analysis/"

# Get the list of folders matching the naming pattern
folders <- list.dirs(parent_dir, full.names = TRUE, recursive = FALSE)
pattern <- "^diff\\.([a-zA-Z0-9]+)t(\\d+)_vs_([a-zA-Z0-9]+)t(\\d+)$"

#folder<- folders[37]

table_promoters<- data.frame()
table_regions<- data.frame()
table_any_regions<- data.frame()


table_promoters_g<- data.frame()
table_regions_g<- data.frame()
table_any_regions_g<- data.frame()

# Loop through the folders
for (folder in folders) {
  if (grepl(pattern, basename(folder))) {
    # Extract g1, g2, t1, t2 values
    matches <- regmatches(basename(folder), regexec(pattern, basename(folder)))[[1]]
    g1 <- matches[2]
    t1 <- matches[3]
    g2 <- matches[4]
    t2 <- matches[5]
    
    # Read the required CSV file
    file_path <- file.path(folder, paste0("diff.", g1, "t",t1, "_vs_", g2, "t",t2, ".results_allpeaks_annotation.csv"))
    if (file.exists(file_path)) {
      data <- read.csv(file_path)
      
      # Extract necessary columns
      filtered_data <- data %>%
        select(Nearest.PromoterID, Annotation, padj, log2FoldChange) 
      filtered_data$Annotation<- sub("\\s.*", "", filtered_data$Annotation)
      
      filtered_data_promoter<- filtered_data[filtered_data$Annotation %in% c("TTS","promoter-TSS"), c(1,4,3)]
      
      filtered_data_promoter<-filtered_data_promoter %>% group_by(Nearest.PromoterID) %>% 
        slice_min(padj) %>% slice_max(abs(log2FoldChange))
      
      filtered_data_any_region<- filtered_data %>% group_by(Nearest.PromoterID) %>% 
        slice_min(padj) %>% slice_max(abs(log2FoldChange))
      
      filtered_data_any_region$gene<- filtered_data_any_region$Nearest.PromoterID
      
      filtered_data_any_region<- filtered_data_any_region[, c("gene", "log2FoldChange", "padj")]
      
      filtered_data<- filtered_data %>% group_by(Nearest.PromoterID, Annotation) %>% 
        slice_min(padj) %>% slice_max(abs(log2FoldChange))
      
      filtered_data$gene<- paste0(filtered_data$Nearest.PromoterID, "_",filtered_data$Annotation)
      
      
      #filtered_data2<- filtered_data %>% group_by(Nearest.PromoterID) %>% 
      #  slice_min(padj) %>% slice_max(abs(log2FoldChange))
      # Handle the case when t1 == t2
      if (t1 == t2) {
        # Create the first modified table for t1 = t2
        
        table1_t1_eq_t2_tss<- filtered_data_promoter
        colnames(table1_t1_eq_t2_tss) <- c("gene", paste0("t", t1, "_",g1, "_",g2, "_log2FoldChange"), paste0("t", t1, "_",g1, "_",g2, "_padj"))
        if(ncol(table_promoters)>0){
          table_promoters<- merge(table_promoters, table1_t1_eq_t2_tss, by="gene", all=T)
        } else{
          table_promoters<- table1_t1_eq_t2_tss
        }
        
        
        # Create the second modified table for t1 = t2
        table1_t1_eq_t2_regions<- filtered_data[, c("gene", "log2FoldChange", "padj")]
        colnames(table1_t1_eq_t2_regions) <- c("gene", paste0("t", t1, "_",g1, "_",g2, "_log2FoldChange"), paste0("t", t1, "_",g1, "_",g2, "_padj"))
        if(ncol(table_regions)>0){
          table_regions<- merge(table_regions, table1_t1_eq_t2_regions, by="gene", all=T)
        } else{
          table_regions<- table1_t1_eq_t2_regions
        }
        
        
        
        # Create the third modified table for t1 = t2
        
        table1_t1_eq_t2_any_region<- filtered_data_any_region
        colnames(table1_t1_eq_t2_any_region) <- c("gene", paste0("t", t1, "_",g1, "_",g2, "_log2FoldChange"), paste0("t", t1, "_",g1, "_",g2, "_padj"))
        if(ncol(table_any_regions)>0){
          table_any_regions<- merge(table_any_regions, table1_t1_eq_t2_any_region, by="gene", all=T)
        } else{
          table_any_regions<- table1_t1_eq_t2_any_region
        }
        
        
      }
      
      # Handle the case when g1 == g2
      if (g1 == g2) {
        # Create the first modified table for t1 = t2
        table1_g1_eq_g2_tss<- filtered_data_promoter
        colnames(table1_g1_eq_g2_tss) <- c("gene", paste0("t", t1, "_t",t2, "_",g2, "_log2FoldChange"), paste0("t", t1, "_t",t2, "_",g2, "_padj"))
        if(ncol(table_promoters_g)>0){
          table_promoters_g<- merge(table_promoters_g, table1_g1_eq_g2_tss, by="gene", all=T)
        } else{
          table_promoters_g<- table1_g1_eq_g2_tss
        }
        
        
        # Create the second modified table for t1 = t2
        
        table1_g1_eq_g2_regions<- filtered_data[, c("gene", "log2FoldChange", "padj")]
        colnames(table1_g1_eq_g2_regions) <- c("gene", paste0("t", t1, "_t",t2, "_",g2, "_log2FoldChange"), paste0("t", t1, "_t",t2, "_",g2, "_padj"))
        if(ncol(table_regions_g)>0){
          table_regions_g<- merge(table_regions_g, table1_g1_eq_g2_regions, by="gene", all=T)
        } else{
          table_regions_g<- table1_g1_eq_g2_regions
        }
        
        # Create the third modified table for t1 = t2
        
        table1_g1_eq_g2_any_region<- filtered_data_any_region[, c("gene", "log2FoldChange", "padj")]
        colnames(table1_g1_eq_g2_any_region) <- c("gene", paste0("t", t1, "_",g1, "_",g2, "_log2FoldChange"), paste0("t", t1, "_",g1, "_",g2, "_padj"))
        if(ncol(table_any_regions_g)>0){
          table_any_regions_g<- merge(table_any_regions_g, table1_g1_eq_g2_any_region, by="gene", all=T)
        } else{
          table_any_regions_g<- table1_g1_eq_g2_any_region
        }
        
        
      }
      
    }
  }
}

#save tables
#table_promoters[]
write.csv(table_promoters, "~/iowa_state/lab/Chromatin_accessibility_paper/epistatic_inputs/table_promoters_epistasis.csv", row.names = F)
write.csv(table_regions, "~/iowa_state/lab/Chromatin_accessibility_paper/epistatic_inputs/table_regions_epistasis.csv", row.names = F)
write.csv(table_any_regions, "~/iowa_state/lab/Chromatin_accessibility_paper/epistatic_inputs/table_any_regions_epistasis.csv", row.names = F)

write.csv(table_promoters_g, "~/iowa_state/lab/Chromatin_accessibility_paper/epistatic_inputs/table_promoters_g_timecourse.csv", row.names = F)
write.csv(table_regions_g, "~/iowa_state/lab/Chromatin_accessibility_paper/epistatic_inputs/table_regions_g_timecourse.csv", row.names = F)
write.csv(table_any_regions_g, "~/iowa_state/lab/Chromatin_accessibility_paper/epistatic_inputs/table_any_regions_g_timecourse.csv", row.names = F)

table_promoters<- read.csv("~/iowa_state/lab/Chromatin_accessibility_paper/epistatic_inputs/table_promoters_epistasis.csv")
table_regions<- read.csv("~/iowa_state/lab/Chromatin_accessibility_paper/epistatic_inputs/table_regions_epistasis.csv")
table_any_regions<- read.csv("~/iowa_state/lab/Chromatin_accessibility_paper/epistatic_inputs/table_any_regions_epistasis.csv")


#make table with stats as region type and distance to TSS, peak size

# Load required libraries
library(dplyr)

# Define the path to the parent directory containing folders
parent_dir <- "~/iowa_state/lab/Chromatin_accessibility_paper/ATACSeq/Analysis results feb 2025/Differential_analysis/Differential_analysis/"

# Get the list of folders matching the naming pattern
folders <- list.dirs(parent_dir, full.names = TRUE, recursive = FALSE)
pattern <- "^diff\\.([a-zA-Z0-9]+)t(\\d+)_vs_([a-zA-Z0-9]+)t(\\d+)$"
example<- read.csv("~/iowa_state/lab/Chromatin_accessibility_paper/ATACSeq/Analysis results feb 2025/Differential_analysis/Differential_analysis/diff.bln1t0_vs_bln1t16/diff.bln1t0_vs_bln1t16.results_allpeaks_annotation.csv")

#folder<- folders[37]
table_all_regions_meta_padj005<- data.frame()
table_all_regions_meta_padj001<- data.frame()
#colnames 
#peak	baseMean	log2FoldChange	lfcSE	stat	pvalue	padj	group	PeakID	Chr	Start	End	Strand	Peak.Score	Focus.Ratio.Region.Size	Annotation	Detailed.Annotation	Distance.to.TSS	Nearest.PromoterID	Entrez.ID	Nearest.Unigene	Nearest.Refseq	Nearest.Ensembl	Gene.Name	Gene.Alias	Gene.Description	Gene.Type

# Loop through the folders
for (folder in folders) {
  if (grepl(pattern, basename(folder))) {
    # Extract g1, g2, t1, t2 values
    matches <- regmatches(basename(folder), regexec(pattern, basename(folder)))[[1]]
    g1 <- matches[2]
    t1 <- matches[3]
    g2 <- matches[4]
    t2 <- matches[5]
    
    # Read the required CSV file
    file_path <- file.path(folder, paste0("diff.", g1, "t",t1, "_vs_", g2, "t",t2, ".results_allpeaks_annotation.csv"))
    if (file.exists(file_path)) {
      data <- read.csv(file_path)
      
      # Extract necessary columns
      filtered_data <- data %>%
        select(Nearest.PromoterID, Annotation, padj, log2FoldChange) 
      filtered_data$Annotation<- sub("\\s.*", "", filtered_data$Annotation)
      
      filtered_data_padj005<- filtered_data[filtered_data$padj<0.005 & !is.na(filtered_data$padj),]
      if(nrow(filtered_data_padj005)>0){
        table_padj005<- data.frame(table(filtered_data_padj005$Annotation))
        colnames(table_padj005) <- c("Annotation", paste0("t", t1, "_", t2,"_",g1, "_",g2, "_padj"))
        if(ncol(table_all_regions_meta_padj005)>0){
          table_all_regions_meta_padj005<- merge(table_all_regions_meta_padj005, table_padj005, by="Annotation", all=T)
          } else{
          table_all_regions_meta_padj005<- table_padj005
          }
        
      }
      
      filtered_data_padj001<- filtered_data[filtered_data$padj<0.001 & !is.na(filtered_data$padj),]
      if(nrow(filtered_data_padj001)>0){
        table_padj001<- data.frame(table(filtered_data_padj001$Annotation))
        colnames(table_padj001) <- c("Annotation", paste0("t", t1, "_", t2,"_",g1, "_",g2, "_padj"))
        if(ncol(table_all_regions_meta_padj001)>0){
          table_all_regions_meta_padj001<- merge(table_all_regions_meta_padj001, table_padj001, by="Annotation", all=T)
          } else{
          table_all_regions_meta_padj001<- table_padj001
          }
      }
      
      
    }
  }
}


#save tables
#table_promoters[]
write.csv(table_all_regions_meta_padj001, "~/iowa_state/lab/Chromatin_accessibility_paper/ATACSeq/stats/table_all_regions_meta_padj001.csv", row.names = F)

write.csv(table_all_regions_meta_padj005, "~/iowa_state/lab/Chromatin_accessibility_paper/ATACSeq/stats/table_all_regions_meta_padj005.csv", row.names = F)


#distance to TSS

# Load required libraries
library(dplyr)

# Define the path to the parent directory containing folders
parent_dir <- "~/iowa_state/lab/Chromatin_accessibility_paper/ATACSeq/Analysis results feb 2025/Differential_analysis/Differential_analysis/"

# Get the list of folders matching the naming pattern
folders <- list.dirs(parent_dir, full.names = TRUE, recursive = FALSE)
pattern <- "^diff\\.([a-zA-Z0-9]+)t(\\d+)_vs_([a-zA-Z0-9]+)t(\\d+)$"

data<- read.csv("~/iowa_state/lab/Chromatin_accessibility_paper/ATACSeq/Analysis results feb 2025/Differential_analysis/Differential_analysis/diff.bln1t0_vs_bln1t16/diff.bln1t0_vs_bln1t16.results_allpeaks_annotation.csv")

#Distance.to.TSS
#folder<- folders[37]
#peak	baseMean	log2FoldChange	lfcSE	stat	pvalue	padj	group	PeakID	Chr	Start	End	Strand	Peak.Score	Focus.Ratio.Region.Size	Annotation	Detailed.Annotation	Distance.to.TSS	Nearest.PromoterID	Entrez.ID	Nearest.Unigene	Nearest.Refseq	Nearest.Ensembl	Gene.Name	Gene.Alias	Gene.Description	Gene.Type

table_promoters<- data.frame()

# Loop through the folders
for (folder in folders) {
  if (grepl(pattern, basename(folder))) {
    # Extract g1, g2, t1, t2 values
    matches <- regmatches(basename(folder), regexec(pattern, basename(folder)))[[1]]
    g1 <- matches[2]
    t1 <- matches[3]
    g2 <- matches[4]
    t2 <- matches[5]
    
    # Read the required CSV file
    file_path <- file.path(folder, paste0("diff.", g1, "t",t1, "_vs_", g2, "t",t2, ".results_allpeaks_annotation.csv"))
    if (file.exists(file_path)) {
      data <- read.csv(file_path)
      
      # Extract necessary columns
      filtered_data <- data %>%
        select(Nearest.PromoterID, Annotation, padj, log2FoldChange, Distance.to.TSS) 
      filtered_data$Annotation<- sub("\\s.*", "", filtered_data$Annotation)
      
      filtered_data_promoter<- filtered_data[filtered_data$Annotation %in% c("TTS","promoter-TSS"), ]
      
      filtered_data_promoter<-filtered_data_promoter %>% group_by(Nearest.PromoterID) %>% 
        slice_min(padj) %>% slice_max(abs(log2FoldChange))
      
      #filtered_data2<- filtered_data %>% group_by(Nearest.PromoterID) %>% 
      #  slice_min(padj) %>% slice_max(abs(log2FoldChange))
      # Handle the case when t1 == t2
      if (t1 == t2) {
        # Create the first modified table for t1 = t2
        
        table1_t1_eq_t2_tss<- filtered_data_promoter[,c(1,3,5)]
        colnames(table1_t1_eq_t2_tss) <- c("gene", paste0("t", t1, "_",g1, "_",g2, "_padj"), paste0("t", t1, "_",g1, "_",g2, "_distanceTSS"))
        if(ncol(table_promoters)>0){
          table_promoters<- merge(table_promoters, table1_t1_eq_t2_tss, by="gene", all=T)
        } else{
          table_promoters<- table1_t1_eq_t2_tss
        }
        
      }
      
    }
  }
}

#save tables
#table_promoters[]
write.csv(table_promoters, "~/iowa_state/lab/Chromatin_accessibility_paper/epistatic_inputs/table_promoters_distanceTSS.csv", row.names = F)







#pbr1 HORVU.MOREX.r3.3HG0320310

#time=c(0, 16, 20, 24, 32, 48)
#gen=c("wt", "mla6", "rar3","bln1", "dm")
#combi <-combn(c("wt", "mla6", "rar3","bln1", "dm"),2)
#combi <-combn(c("wt", "mla6", "rar3","bln1", "dm"),2)

time=c(0, 16, 20, 32)
gen=c("wt", "mla6", "bln1", "dm")
combi <-combn(c("wt", "mla6","bln1", "dm"),2)
comb<- apply(combi, 2, FUN=function(x)paste0(x[1],"_", x[2]))
FC_padj_table<- read.csv("~/iowa_state/lab/Chromatin_accessibility_paper/epistatic_inputs/table_any_regions_epistasis.csv")

#FC_padj_table<- read.csv("~/iowa_state/lab/Chromatin_accessibility_paper/epistatic_inputs/table_promoters_epistasis.csv")
#FC_padj_table<- read.csv("~/iowa_state/lab/RNAseq/DESeq2 HvR3/hv_R3_genes_deseq2_results_pairwise_genotype_logfc_padj_tax_sp.csv")

#FC_padj_table<- read.csv("~/iowa_state/lab/RNAseq/RNASeq_paper/vvz manuscript/analysis/hv_R2_genes_deseq2_results_pairwise_genotype_logfc_padj_tax_sp.csv")
#FC_padj_table<-read.csv("~/iowa_state/lab/RNAseq/DESeq2 analysis HvV2/hv_V2_genes_deseq2_results_pairwise_genotype_FC_padj_tax_sp.csv", stringsAsFactors = F)
FC_padj_table[,grep("FoldChange", colnames(FC_padj_table))][is.na(FC_padj_table[,grep("FoldChange", colnames(FC_padj_table))])]<- 0
FC_padj_table[,grep("padj", colnames(FC_padj_table))][is.na(FC_padj_table[,grep("padj", colnames(FC_padj_table))])]<- 1


# deciding on a threshold for significant comparisons
library(tidyverse)
thr_vals<- c(0.05, 0.01, 0.005, 0.003,0.0025, 0.001, 0.0008)
sig_tests<- data.frame(org=c(rep("Hv", 24)), test=colnames(FC_padj_table)[grep("padj", colnames(FC_padj_table))])
for (names in colnames(FC_padj_table)[grep("padj", colnames(FC_padj_table))]) {
  for(thr in thr_vals){
    sig_tests[sig_tests$test==names, paste0("n_",thr)]<- sum(FC_padj_table[,names]<thr, na.rm = T)*thr
    
  }
}

sig_tests %>% group_by(org) %>% summarise(mean_0.05=mean(n_0.05), mean_0.01=mean(n_0.01), mean_0.005=mean(n_0.005),
                                          mean_0.003=mean(n_0.003), mean_0.0025=mean(n_0.0025), mean_0.001=mean(n_0.001),
                                          mean_0.0008=mean(`n_8e-04`))


# A tibble: 1 × 8
# org   mean_0.05 mean_0.01 mean_0.005 mean_0.003 mean_0.0025 mean_0.001 mean_0.0008
# <chr>     <dbl>     <dbl>      <dbl>      <dbl>       <dbl>      <dbl>       <dbl>
#   1 Hv         27.4      4.09       1.82       1.01       0.821      0.287       0.224

#create a function that compares using the 2 genotypes, and p value, report with FC and fix by num and den, and if DE or !DE

DE_list<- function(DE_table, num, den, padj, DE=T, FC=0, comb, times){
  library(textclean)
  #paste0(num, "_", den) %in% comb
  if(  sum(grepl(paste0(num, "_", den), colnames(DE_table)))>0){
    DE_table_mod<- DE_table[, c(1, grep(paste(paste0(num, "_", den), collapse = "|"), colnames(DE_table)))]
  }else{
    DE_table_mod<- DE_table[, c(1, grep(paste(paste0(den, "_", num), collapse = "|"), colnames(DE_table)))]
    DE_table_mod[, grep("FoldChange", colnames(DE_table_mod))]<- -DE_table_mod[, grep("FoldChange", colnames(DE_table_mod))]
    colnames(DE_table_mod)<- swap(colnames(DE_table_mod), num, den)
  }
  #DE_table_mod<- DE_table_mod[complete.cases(DE_table_mod),]
  result<- list()
  for(t in times){
    if(DE==T){
      if(FC>0){
        result[[paste0("t",t)]]<- DE_table_mod[DE_table_mod[,paste0("t",t,"_", num, "_", den,"_padj")]< padj &
                                                 DE_table_mod[,paste0("t",t,"_", num, "_", den,"_log2FoldChange")] > 0, c(1, grep(paste0("t",t,"_", num, "_", den), colnames(DE_table_mod)))]
      }else if(FC<0){
        result[[paste0("t",t)]]<- DE_table_mod[DE_table_mod[,paste0("t",t,"_", num, "_", den,"_padj")]< padj &
                                                 DE_table_mod[,paste0("t",t,"_", num, "_", den,"_log2FoldChange")] < 0, c(1, grep(paste0("t",t,"_", num, "_", den), colnames(DE_table_mod)))]
      }else{
        result[[paste0("t",t)]]<- DE_table_mod[DE_table_mod[,paste0("t",t,"_", num, "_", den,"_padj")]< padj, c(1, grep(paste0("t",t,"_", num, "_", den), colnames(DE_table_mod)))]
      }
    }else{
      result[[paste0("t",t)]]<- DE_table_mod[DE_table_mod[,paste0("t",t,"_", num, "_", den,"_padj")]> padj, c(1, grep(paste0("t",t,"_", num, "_", den), colnames(DE_table_mod)))]
      #result[[paste0("t",t)]]<- result[[paste0("t",t)]][!duplicated(result[[paste0("t",t)]]$gene),]
    }
    #result[[paste0("t",t)]]<- result[[paste0("t",t)]][!is.na(result[[paste0("t",t)]]$gene),]
  }
  return(result)
}

#z<- result[["t20"]]

set_operations<- function(list1, list2, operation="intersection"){
  result<- list()
  if(operation=="intersection"){
    for(t in names(list1)){
      result[[t]]<- list1[[t]][list1[[t]]$gene %in% list2[[t]]$gene, ]
      if(sum(!colnames(list2[[t]]) %in% colnames(list1[[t]]))>0){
        result[[t]]<- cbind(result[[t]], list2[[t]][match(result[[t]]$gene, list2[[t]]$gene),-grep(paste(colnames(list1[[t]]), collapse = "|"), colnames(list2[[t]]))])
      }
    }
  }else if(operation=="minus"){
    for(t in names(list1)){
      result[[t]]<- list1[[t]][!list1[[t]]$gene %in% list2[[t]]$gene, ]
    }
  }else if(operation=="union"){
    library(dplyr)
    for(t in names(list1)){
      common_genes<- list1[[t]][list1[[t]]$gene %in% list2[[t]]$gene, "gene"]
      result[[t]]<- bind_rows(merge(list1[[t]][list1[[t]]$gene %in% common_genes, ], list2[[t]][match(common_genes, list2[[t]]$gene),], by=colnames(list1[[t]])[colnames(list1[[t]]) %in% colnames(list2[[t]])]),
                          list1[[t]][!list1[[t]]$gene %in% common_genes, ], list2[[t]][!list2[[t]]$gene %in% common_genes, ])
    }
  }
  return(result)
}

get_epistatic<-function(triple_union_list, thr=1, thr_prop=0.2, mode="dev"){
  result<- triple_union_list
  result_notEpi<-list()
  for(t in names(result)){
    result[[t]][,grep("FoldChange", colnames(result[[t]]))][is.na(result[[t]][,grep("FoldChange", colnames(result[[t]]))])]<- 0
    result[[t]]$expected<-result[[t]][,4]+result[[t]][,6]
    result[[t]]$deviation<-result[[t]][,2]-result[[t]][,4]-result[[t]][,6]
    result[[t]]$dev_prop<-abs(result[[t]]$deviation/(result[[t]][,2]+0.1))
    if(mode=="dev"){
      result_notEpi[[t]]<- result[[t]][abs(result[[t]]$deviation)<thr,]
      result[[t]]<- result[[t]][abs(result[[t]]$deviation)>thr,]
    }else{
      result_notEpi[[t]]<- result[[t]][result[[t]]$dev_prop<thr_prop,]
      result[[t]]<- result[[t]][result[[t]]$dev_prop>thr_prop,]
    }
  }
  return(list(result,result_notEpi))
}

#for barley
#horvu_DE_table<- FC_padj_table[grep("HORVU",FC_padj_table$gene), ]
horvu_DE_table<- FC_padj_table
sum(FC_padj_table$t0_dm_bln1_padj<0.05)

#Mla6 and Bln1
padj<- 0.001
padj<- 0.005 #selected with expeted number of significant genes

wt_mla6_DE<- DE_list(DE_table=horvu_DE_table, num="wt", den="mla6", padj=padj, DE=T, comb=comb, times=time)
wt_bln1_DE<- DE_list(DE_table=horvu_DE_table, num="wt", den="bln1", padj=padj, DE=T, comb=comb, times=time)
wt_dm_DE<- DE_list(DE_table=horvu_DE_table, num="wt", den="dm", padj=padj, DE=T, comb=comb, times=time)

mla6_bln1_notDE<- DE_list(DE_table=horvu_DE_table, num="mla6", den="bln1", padj=padj, DE=F, comb=comb, times=time)
mla6_dm_notDE<- DE_list(DE_table=horvu_DE_table, num="mla6", den="dm", padj=padj, DE=F, comb=comb, times=time)
bln1_dm_notDE<- DE_list(DE_table=horvu_DE_table, num="bln1", den="dm", padj=padj, DE=F, comb=comb, times=time)

#shared_eff_Mla6_Bln1<- set_operations(list1=wt_mla6_DE, list2=wt_bln1_DE, operation="intersection")
#shared_eff_Mla6_Bln1<- set_operations(list1=shared_eff_Mla6_Bln1, list2=mla6_bln1_notDE, operation="intersection")
epi_Mla6_Bln1_DE<-set_operations(list1=wt_dm_DE, list2=wt_mla6_DE, operation="union")
epi_Mla6_Bln1_DE<-set_operations(list1=epi_Mla6_Bln1_DE, list2=wt_bln1_DE, operation="union")
#change this so only er look at DE genes in wt dm comparisson
#epi_Mla6_Bln1<-set_operations(list1=epi_Mla6_Bln1, list2=wt_bln1_DE, operation="intersection")

epi_Mla6_Bln1_all_prop<- get_epistatic(epi_Mla6_Bln1_DE, thr_prop = 0.5, mode = "prop")
a1<-epi_Mla6_Bln1_all_prop[[1]][["t16"]]
epi_Mla6_Bln1_prop<-epi_Mla6_Bln1_all_prop[[1]]
additive_Mla6_Bln1_prop<-epi_Mla6_Bln1_all_prop[[2]]

epi_Mla6_Bln1_all_dev<- get_epistatic(epi_Mla6_Bln1_DE, thr = 1)
epi_Mla6_Bln1<-epi_Mla6_Bln1_all_dev[[1]]
additive_Mla6_Bln1<-epi_Mla6_Bln1_all_dev[[2]]


epi_Mla6_Bln1<-epi_Mla6_Bln1_all_prop[[1]]
additive_Mla6_Bln1<-epi_Mla6_Bln1_all_prop[[2]]

a<-epi_Mla6_Bln1[["t16"]]
b<-additive_Mla6_Bln1[["t20"]]

mla6_dm_DE_FC_plus<- DE_list(DE_table=horvu_DE_table, num="mla6", den="dm", padj=padj, DE=T, FC=1, comb=comb, times=time)
mla6_dm_DE_FC_minus<- DE_list(DE_table=horvu_DE_table, num="mla6", den="dm", padj=padj, DE=T, FC=-1, comb=comb, times=time)

a1<-mla6_dm_DE_FC_plus[["t20"]]
b1<-mla6_dm_DE_FC_minus[["t20"]]

bln1_dm_DE_FC_plus<- DE_list(DE_table=horvu_DE_table, num="bln1", den="dm", padj=padj, DE=T, FC=1, comb=comb, times=time)
bln1_dm_DE_FC_minus<- DE_list(DE_table=horvu_DE_table, num="bln1", den="dm", padj=padj, DE=T, FC=-1, comb=comb, times=time)

positive<- set_operations(list1=epi_Mla6_Bln1, list2=bln1_dm_DE_FC_minus, operation="intersection")
positive<- set_operations(list1=positive, list2=mla6_dm_DE_FC_minus, operation="intersection")

apos<-positive[["t16"]]

pseudo_masked<- set_operations(list1=epi_Mla6_Bln1, list2=bln1_dm_DE_FC_plus, operation="intersection")
pseudo_masked<- set_operations(list1=pseudo_masked, list2=mla6_dm_DE_FC_minus, operation="intersection")

pseudo_masked2<- set_operations(list1=epi_Mla6_Bln1, list2=bln1_dm_DE_FC_minus, operation="intersection")
pseudo_masked2<- set_operations(list1=pseudo_masked2, list2=mla6_dm_DE_FC_plus, operation="intersection")
pseudo_masked<- set_operations(list1=pseudo_masked, list2=pseudo_masked2, operation="union")


d<-bln1_dm_DE_FC_plus[["t20"]]
bPM<-pseudo_masked[["t20"]]

negative<- set_operations(list1=epi_Mla6_Bln1, list2=bln1_dm_DE_FC_plus, operation="intersection")
negative<- set_operations(list1=negative, list2=mla6_dm_DE_FC_plus, operation="intersection")

c<-negative[["t20"]]

#masked<- set_operations(list1=wt_dm_DE, list2=wt_mla6_DE, operation="intersection")
#masked<- set_operations(list1=masked, list2=shared_eff_Mla6_Bln1, operation="minus")
masked<- set_operations(list1=epi_Mla6_Bln1, list2=mla6_dm_notDE, operation="intersection")

e<-masked[["t20"]]

#suppression<- set_operations(list1=wt_bln1_DE, list2=wt_dm_DE, operation="intersection")
#suppression<- set_operations(list1=suppression, list2=shared_eff_Mla6_Bln1, operation="minus")
suppression<- set_operations(list1=epi_Mla6_Bln1, list2=bln1_dm_notDE, operation="intersection")

f<-suppression[["t20"]]

sum(e$gene %in% f$gene)
g<-e$gene[e$gene %in% f$gene]

#If we want to separate further symmetric interaction and remove those from masked and suppression
symmetric<- set_operations(list1=masked, list2=suppression, operation="intersection")
masked<- set_operations(list1=masked, list2=symmetric, operation="minus")
suppression<- set_operations(list1=suppression, list2=symmetric, operation="minus")

#intersections
#make a function to classify the genes with the epistasis patters by timepoint
#find if one gene has more than one pattern

combine_lists<- function(oldlist){
  library(data.table)
  new_list<- list()
  for (name in names(oldlist)) {
    if(nrow(oldlist[[name]]>0)){
      new_list[[name]]<- data.frame(gene=oldlist[[name]]$gene, time=name, epistasis=deparse(substitute(oldlist)))
    }
  }
  combined<- rbindlist(new_list)
  return(combined)
}


library(tidyverse)
addi_Mla6_Bln1<- combine_lists(additive_Mla6_Bln1)
pos<- combine_lists(positive)
pseudo<- combine_lists(pseudo_masked)
neg<- combine_lists(negative)
masked_epi<- combine_lists(masked)
supp<- combine_lists(suppression)
symmetric_Mla6_Bln1<- combine_lists(symmetric)

#combined_tables_Mla6_Bln1<- rbind(addi_Mla6_Bln1, pos, pseudo, neg, masked_epi, supp)
combined_tables_Mla6_Bln1<- rbind(addi_Mla6_Bln1, pos, pseudo, neg, masked_epi, supp, symmetric_Mla6_Bln1)
combined_tables_Mla6_Bln1<- combined_tables_Mla6_Bln1 %>% group_by(gene, time) %>% summarise(epistasis = paste0(epistasis, collapse = ";"))
spread_tables_Mla6_Bln1<- spread(combined_tables_Mla6_Bln1, key = time, value = epistasis)

summary_barley_Mla6_Bln1<- sapply(X = spread_tables_Mla6_Bln1[,-1], FUN = table)


#spread_tables_Mla6_Bln1[spread_tables_Mla6_Bln1=="additive_Mla6_Bln1"]<-""
spread_tables_Mla6_Bln1$consensus<- apply(spread_tables_Mla6_Bln1[,2:5], 1, 
                                          FUN=function(x){y=table(unlist(x))
                                          consensus=names(y)[which(y==max(y))]
                                          ifelse("additive_Mla6_Bln1" %in% consensus & length(consensus)>1,
                                                 paste(consensus[!consensus=="additive_Mla6_Bln1"], collapse = ";"),
                                                 paste(consensus, collapse = ";"))})

unique(spread_tables_Mla6_Bln1$consensus)
spread_tables_Mla6_Bln1<- spread_tables_Mla6_Bln1[spread_tables_Mla6_Bln1$consensus!="",]
table(spread_tables_Mla6_Bln1$consensus)
#spread_tables_Mla6_Bln1[spread_tables_Mla6_Bln1$consensus=="", "consensus"]<- "additive_Mla6_Bln1"
#spread_tables_Mla6_Bln1[is.na(spread_tables_Mla6_Bln1)]<- "additive_Mla6_Bln1"

#save the lists
#deparse(substitute(inter_common_epi_mla6_Bln1))
save_lists<- function(list, path){
  for (name in names(list)) {
    if(nrow(list[[name]]>0)){
    write.csv(list[[name]], paste0(path, name,"_",deparse(substitute(list)),".csv"), row.names = F)
    }}
  saveRDS(list, paste0(path, deparse(substitute(list)),".RDS"))
}

#HVR3

save_path<- "~/iowa_state/lab/Chromatin_accessibility_paper/epistatic_inputs/padj005" 
save_path<- "~/iowa_state/lab/Chromatin_accessibility_paper/epistatic_inputs/any_region_padj005/" 

save_lists(list=epi_Mla6_Bln1, path=save_path)
save_lists(list=additive_Mla6_Bln1, path=save_path)
save_lists(list=positive, path=save_path)
save_lists(list=pseudo_masked, path=save_path)
save_lists(list=negative, path=save_path)
save_lists(list=masked, path=save_path)
save_lists(list=suppression, path=save_path)
save_lists(list=symmetric, path=save_path)

write.csv(spread_tables_Mla6_Bln1, paste0(save_path, "spread_tables_Mla6_Bln1.csv"))
write.csv(summary_barley_Mla6_Bln1[2:4], paste0(save_path, "summary_barley_Mla6_Bln1.csv"))

summary_barley_Mla6_Bln1[1]



########################################
#position analysis


#create the lists for chr location
#directory<- "~/iowa_state/lab/RNAseq/RNASeq_paper/vvz manuscript/analysis/chr analysis/R2/Genelists for positioning/"
directory<- "~/iowa_state/lab/Chromatin_accessibility_paper/epistatic_inputs/chr_analysis/gene_lists/"
directory<- "~/iowa_state/lab/Chromatin_accessibility_paper/epistatic_inputs/padj005/chr_analysis/gene_lists/"
directory<- "~/iowa_state/lab/Chromatin_accessibility_paper/epistatic_inputs/any_region_padj005/chr_analysis/gene_lists/"

#annotation<- read.table("~/iowa_state/lab/RNAseq/DESeq2 HvR2/annotation_HvR2.txt", sep=';', stringsAsFactors = F, quote = "\"")[,1:2]
annotation<- read.csv("~/iowa_state/lab/genome annotation files/tritex R3/HR3_annotation.csv")[,c(3,6)]
colnames(annotation)<-c("gene","description")

spread_tables_Mla6_Bln1<- read.csv("~/iowa_state/lab/Chromatin_accessibility_paper/epistatic_inputs/spread_tables_Mla6_Bln1.csv", row.names = 1)
spread_tables_Mla6_Bln1<- read.csv("~/iowa_state/lab/Chromatin_accessibility_paper/epistatic_inputs/padj005/spread_tables_Mla6_Bln1.csv", row.names = 1)
spread_tables_Mla6_Bln1<- read.csv("~/iowa_state/lab/Chromatin_accessibility_paper/epistatic_inputs/any_region_padj005/spread_tables_Mla6_Bln1.csv", row.names = 1)

epistatic_spread_tables_Mla6_Bln1<- spread_tables_Mla6_Bln1[spread_tables_Mla6_Bln1$consensus!="additive_Mla6_Bln1", ]
epistatic_spread_tables_Mla6_Bln1<- merge(epistatic_spread_tables_Mla6_Bln1, annotation, by="gene", all.x=T)
for(type in unique(epistatic_spread_tables_Mla6_Bln1$consensus)){
  table<- epistatic_spread_tables_Mla6_Bln1[epistatic_spread_tables_Mla6_Bln1$consensus==type, c(1,7,2:6)]
  write.csv(table, paste0(directory, "Barley_Mla6_Bln1_", type, "_Genelist.csv"))
}

write.csv(paste0("Barley_Mla6_Bln1_",unique(epistatic_spread_tables_Mla6_Bln1$consensus)),
          "~/iowa_state/lab/Chromatin_accessibility_paper/epistatic_inputs/chr_analysis/trts_barley.csv")

write.csv(paste0("Barley_Mla6_Bln1_",unique(epistatic_spread_tables_Mla6_Bln1$consensus)),
          "~/iowa_state/lab/Chromatin_accessibility_paper/epistatic_inputs/padj005/chr_analysis/trts_barley.csv")

write.csv(paste0("Barley_Mla6_Bln1_",unique(epistatic_spread_tables_Mla6_Bln1$consensus)),
          "~/iowa_state/lab/Chromatin_accessibility_paper/epistatic_inputs/any_region_padj005/chr_analysis/trts_barley.csv")


#HG test
HG_test<- function(trts, bin_table, overlap_table){
  #table trts by bins
  hypergeom_test <- data.frame(matrix(ncol = 1+2*length(trts), nrow = nrow(bin_table))) 
  hypergeom_test[,1]<- bin_table[,1]
  colnames(hypergeom_test)<- c("binID", unlist(sapply(trts, function(x) c(paste0(x,"_pvalue"), paste0(x,"_padj")))))
  
  for(m in 1:length(trts)){
    epi_gene_list<- read.csv(paste0(folder,"gene_lists/",trts[m], "_Genelist.csv"), 
                             stringsAsFactors = F, row.names = 1)[,"gene"]
    overlap_table_filtered<- overlap_table[overlap_table$seqname %in% epi_gene_list, c("seqname", colnames(bin_table)[1])]
    DE_genes_bins<-data.frame(table(overlap_table_filtered[,colnames(bin_table)[1]]))
    overlap_table_filtered<- merge(overlap_table_filtered, DE_genes_bins, by.x=colnames(bin_table)[1], by.y="Var1", all.x=T)
    DE_bins<- bin_table[unlist(bin_table[,1]) %in% DE_genes_bins$Var1,]
    
    #total DE in bin, total DE trt, sum not DE in bins listed in trt, genes in bin
    hypergeom_test[match(DE_genes_bins$Var1, unlist(hypergeom_test$binID)), 
                   paste0(trts[m],"_pvalue")]<- phyper(DE_genes_bins$Freq, length(epi_gene_list), 
                                                       sum(DE_bins$Freq - DE_genes_bins$Freq), 
                                                       DE_bins$Freq, lower.tail = FALSE)
    
    hypergeom_test[unlist(hypergeom_test$binID) %in% DE_genes_bins$Var1,paste0(trts[m],"_padj")]<- p.adjust(hypergeom_test[unlist(hypergeom_test$binID) %in% DE_genes_bins$Var1,paste0(trts[m],"_pvalue")], method = "BH")
  }
  hypergeom_test<- hypergeom_test[rowSums(is.na(hypergeom_test[,-1])) != ncol(hypergeom_test[,-1]), ]
  return(hypergeom_test)
}

#test for all bin tables
folder<- "~/iowa_state/lab/Chromatin_accessibility_paper/epistatic_inputs/chr_analysis/"
folder<- "~/iowa_state/lab/Chromatin_accessibility_paper/epistatic_inputs/padj005/chr_analysis/"
folder<- "~/iowa_state/lab/Chromatin_accessibility_paper/epistatic_inputs/any_region_padj005/chr_analysis/"

chrs = list.files(path = paste0(folder,"gene_lists/"), pattern = "csv")
#trts<- read.csv("~/iowa_state/lab/RNAseq/RNASeq_paper/vvz manuscript/analysis/chr analysis/R2/trts_barley.csv", row.names = 1)
trts<- read.csv("~/iowa_state/lab/Chromatin_accessibility_paper/epistatic_inputs/chr_analysis/trts_barley.csv", row.names = 1)
trts<- read.csv("~/iowa_state/lab/Chromatin_accessibility_paper/epistatic_inputs/padj005/chr_analysis/trts_barley.csv", row.names = 1)
trts<- read.csv("~/iowa_state/lab/Chromatin_accessibility_paper/epistatic_inputs/any_region_padj005/chr_analysis/trts_barley.csv", row.names = 1)

trts<- trts$x
#trts<- gsub('.{0,13}$', '', chrs) 
#overlap_table<- read.csv(paste0(folder, "Binned_genes.csv")
#overlap_table<- read.csv(paste0(folder, "Binned_genes_HR3.csv")
complete_overlap_df<- read.csv(paste0("~/iowa_state/lab/RNAseq/RNASeq_paper/vvz manuscript/analysis/chr analysis/R3/", "Binned_genes_HR3.csv"))

overlap_table<- complete_overlap_df

bins_1Mb<- read.csv( paste0("~/iowa_state/lab/RNAseq/RNASeq_paper/vvz manuscript/analysis/chr analysis/R3/", "bins_1Mb.csv"))
bins_10Mb<-read.csv(paste0("~/iowa_state/lab/RNAseq/RNASeq_paper/vvz manuscript/analysis/chr analysis/R3/", "bins_10Mb.csv"))
bins_100Mb<-read.csv( paste0("~/iowa_state/lab/RNAseq/RNASeq_paper/vvz manuscript/analysis/chr analysis/R3/", "bins_100Mb.csv"))
bins_chr<-read.csv( paste0("~/iowa_state/lab/RNAseq/RNASeq_paper/vvz manuscript/analysis/chr analysis/R3/", "bins_chrs.csv"))


HG_1Mb<- HG_test(trts, bins_1Mb, overlap_table)
write.csv(HG_1Mb, paste0(folder, "hypergeometric_test_bins_1Mb.csv"), row.names = F)

HG_10Mb<- HG_test(trts, bins_10Mb, overlap_table)
write.csv(HG_10Mb, paste0(folder, "hypergeometric_test_bins_10Mb.csv"))

HG_100Mb<- HG_test(trts, bins_100Mb, overlap_table)
write.csv(HG_100Mb, paste0(folder, "hypergeometric_test_bins_100Mb.csv"))

HG_chrs<- HG_test(trts, bins_chr, overlap_table)
write.csv(HG_chrs, paste0(folder, "hypergeometric_test_bins_chrs.csv"))


#Make the plots
#https://www.bioconductor.org/packages/release/bioc/vignettes/ggbio/inst/doc/ggbio.pdf
#BiocManager::install("ggbio")

library(ggbio)
#Fig 4a position
library(reshape2)
#define the significant hotspots and plot coloring by pattern
#folder<- "~/iowa_state/lab/RNAseq/RNASeq_paper/vvz manuscript/analysis/chr analysis/R3/"
folder<- "~/iowa_state/lab/Chromatin_accessibility_paper/epistatic_inputs/chr_analysis/"
folder<- "~/iowa_state/lab/Chromatin_accessibility_paper/epistatic_inputs/padj005/chr_analysis/"
folder<- "~/iowa_state/lab/Chromatin_accessibility_paper/epistatic_inputs/any_region_padj005/chr_analysis/"

# HG_1Mb<- read.csv(paste0(folder, "hypergeometric_test_bins_1Mb.csv"))
# HG_10Mb<- read.csv(paste0(folder, "hypergeometric_test_bins_10Mb.csv"), row.names = 1)
HG_100Mb<- read.csv(paste0(folder, "hypergeometric_test_bins_100Mb.csv"), row.names = 1)
#HG_chrs<- read.csv(paste0(folder, "hypergeometric_test_bins_chrs.csv"), row.names = 1)

# type=c("additive_Mla6_Bln1", "symmetric", "masked", "suppression", "pseudo_masked", "positive", "negative",
#        "unique_eff_Mla6", "unique_eff_Rar3", "equal_eff_Mla6_Rar3", "predominant_Mla6", "predominant_Rar3")
type=c("additive_Mla6_Bln1", "symmetric", "masked", "suppression", "pseudo_masked", "positive", "negative")
trts<- paste0("Barley_Mla6_Bln1_",type, "_padj")

HG_100Mb_mla6_bln1<- HG_100Mb[, c(1, grep(paste(trts, collapse = "|"), colnames(HG_100Mb)))]
#thr 0.005
HG_100Mb_mla6_bln1<- HG_100Mb_mla6_bln1[rowSums(HG_100Mb_mla6_bln1[,-1] <0.005, na.rm = T)>0,
                                        c(T,colSums(HG_100Mb_mla6_bln1[,-1] <0.005, na.rm = T)>0)]
sig_trts<- sapply(colnames(HG_100Mb_mla6_bln1)[-1], FUN=function(x)strsplit(x,"_")[[1]][4])
colnames(HG_100Mb_mla6_bln1)[-1]<-sig_trts
HG_100Mb_mla6_bln1<- melt(HG_100Mb_mla6_bln1, 1, na.rm = T)
HG_100Mb_mla6_bln1<- HG_100Mb_mla6_bln1[HG_100Mb_mla6_bln1$value< 0.005,]

#binned_genes<-read.csv(paste0(folder, "Binned_genes_HR3.csv"))
binned_genes<- read.csv(paste0("~/iowa_state/lab/RNAseq/RNASeq_paper/vvz manuscript/analysis/chr analysis/R3/", "Binned_genes_HR3.csv"))

#spread_tables_Mla6_Bln1<- read.csv("~/iowa_state/lab/RNAseq/RNASeq_paper/vvz manuscript/analysis/epistasis_lists/barley/HVR3/spread_tables_Mla6_Bln1.csv", row.names = 1)
#spread_tables_Mla6_Bln1<- read.csv("~/iowa_state/lab/Chromatin_accessibility_paper/epistatic_inputs/spread_tables_Mla6_Bln1.csv", row.names = 1)
#spread_tables_Mla6_Bln1<- read.csv("~/iowa_state/lab/Chromatin_accessibility_paper/epistatic_inputs/spread_tables_Mla6_Bln1.csv", row.names = 1)

spread_tables_Mla6_Bln1<- read.csv("~/iowa_state/lab/Chromatin_accessibility_paper/epistatic_inputs/spread_tables_Mla6_Bln1.csv", row.names = 1)
spread_tables_Mla6_Bln1<- read.csv("~/iowa_state/lab/Chromatin_accessibility_paper/epistatic_inputs/padj005/spread_tables_Mla6_Bln1.csv", row.names = 1)
spread_tables_Mla6_Bln1<- read.csv("~/iowa_state/lab/Chromatin_accessibility_paper/epistatic_inputs/any_region_padj005/spread_tables_Mla6_Bln1.csv", row.names = 1)

spread_tables_Mla6_Bln1<- spread_tables_Mla6_Bln1[spread_tables_Mla6_Bln1$consensus %in% sig_trts, c("gene","consensus")]
spread_tables_Mla6_Bln1<- merge(spread_tables_Mla6_Bln1, binned_genes[,c(5,14)], by.x="gene", by.y="seqname", all.x=T)
sig_genes_Mla6_Bln1<- merge(spread_tables_Mla6_Bln1, HG_100Mb_mla6_bln1, by.x=c("consensus","bin_100Mb"), by.y=c("variable", "binID"))

write.csv(sig_genes_Mla6_Bln1, paste0(folder, "Barley_chr_position_significant.csv"))

library(tidyverse)
library(GenomicRanges)
horvu <- read.csv(paste0("~/iowa_state/lab/RNAseq/RNASeq_paper/vvz manuscript/analysis/chr analysis/R3/", "Horvu_table.csv"),row.names = 1, stringsAsFactors = F)
len<- horvu %>% group_by(V1) %>% summarise(lenght= max(V2, V3))
names(horvu)<- (c("chr","start","end","strand","seqname"))
Seqinfo(as.character(len$V1), seqlengths=len$lenght, isCircular=NA, genome=NA)
horvu_ranges<- makeGRangesFromDataFrame(horvu,keep.extra.columns=TRUE,ignore.strand=FALSE,
                                        seqinfo=Seqinfo(as.character(len$V1), seqlengths=len$lenght, isCircular=NA, genome=NA),
                                        seqnames.field="chr",start.field="start",end.field="end",
                                        strand.field="strand",starts.in.df.are.0based=FALSE)

#Make the plots
#https://www.bioconductor.org/packages/release/bioc/vignettes/ggbio/inst/doc/ggbio.pdf
#BiocManager::install("ggbio")

library(ggbio)
#Now a loop
p <- autoplot(seqinfo(horvu_ranges), layout = "karyogram")
chrs<- sig_trts
colors  <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
names(colors)<- c("additive_Mla6_Bln1", "symmetric", "masked", "suppression", "pseudo_masked", "positive", "negative")

for(m in 1:length(chrs)){
  name<- chrs[m]
  coords<- sig_genes_Mla6_Bln1[sig_genes_Mla6_Bln1$consensus==chrs[m],]
  coords<- merge(horvu, coords, by.x="seqname",by.y="gene")
  coords_ann<-makeGRangesFromDataFrame(coords, keep.extra.columns=T, seqnames.field= "chr",start.field="start",
                                       seqinfo=Seqinfo(as.character(len$V1)), end.field="end", strand.field="strand")
  p<-p + layout_karyogram(coords_ann, color=colors[chrs[m]])
}

pdf(paste0(folder, "Barley_chr_position_significant.pdf"), width = 5, height = 3, fonts = "ArialMT", pointsize = 52)
print(p)
dev.off()

