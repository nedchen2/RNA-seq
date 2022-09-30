
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(gridExtra))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(viridis))
suppressPackageStartupMessages(library(ComplexHeatmap))
suppressPackageStartupMessages(library(circlize))
suppressPackageStartupMessages(library(pheatmap))

path1 <- "LavendustinA-vs-DMSO-all.gene.xls"
path2 <- "Reversine-vs-DMSO-all.gene.xls"

#提取组名1
groupname <- gsub("\\.(txt|xls)$", "", gsub("-all.gene.xls$", "", path1))#opt$input
Case_names <- unlist(strsplit(groupname,"-vs-"))[1]
Control_names <- unlist(strsplit(groupname,"-vs-"))[2]

#提取组名2
groupname2 <- gsub("\\.(txt|xls)$", "", gsub("-all.gene.xls$", "", path2))#opt$input
Case_names2 <- unlist(strsplit(groupname2,"-vs-"))[1]
Control_names2 <- unlist(strsplit(groupname2,"-vs-"))[2]


one <- read.delim(normalizePath(path1), header=T, sep="\t", check.names=F, quote="")
zero <- read.delim(normalizePath(path2), header=T, sep="\t", check.names=F, quote="")

one = plyr::rename(one, c("log2FoldChange"="logFC_A"))
one = plyr::rename(one, c("q-value"="FDR_A"))
zero = plyr::rename(zero, c("log2FoldChange"="logFC_B"))
zero = plyr::rename(zero, c("q-value"="FDR_B"))

two <- inner_join(one, zero, by=c("gene_id" = "gene_id"))

two <- subset(two,select=c("gene_id","logFC_A","logFC_B","FDR_A","FDR_B"))
#标记显著性（默认 p < 0.05）
two[which(two$FDR_A < 0.05 & two$FDR_B < 0.05),'type1'] <- 'sign'
two[which(two$FDR_A >= 0.05 | two$FDR_B >= 0.05),'type1'] <- 'no'

#标记差异倍数（默认 |log2FC| >= 1）
two[which(two$logFC_A <= -1 & two$logFC_B <= -1),'type2'] <- 'a_down.b_down'
two[which(two$logFC_A >= 1 & two$logFC_B <= -1),'type2'] <- 'a_up.b_down'
two[which(two$logFC_A <= -1 & two$logFC_B >= 1),'type2'] <- 'a_down.b_up'
two[which(two$logFC_A >= 1 & two$logFC_B >= 1),'type2'] <- 'a_up_b_up'
two[is.na(two$type2),'type2'] <- 'no'

rownames(two)=two[,1]

showname_q<-paste0("Sig.","(",groupname,")")
showname_q2<-paste0("Sig.","(",groupname2,")")
showname_FC<-paste0("log2FC.","(",groupname,")")
showname_FC2<-paste0("log2FC.","(",groupname2,")")

two = plyr::rename(two, c("FDR_A"=showname_q))
two = plyr::rename(two, c("FDR_B"=showname_q2))
two = plyr::rename(two, c("logFC_A"=showname_FC))
two = plyr::rename(two, c("logFC_B"=showname_FC2))

Three <- subset(two,type1=='sign'& type2!='no',select=c(showname_FC,showname_FC2))
Three <- as.matrix(Three)

Four <- subset(two,type1=='sign'& type2!='no',select=c(showname_q,showname_q2))
Four <- as.matrix(Four)

getComplexHeatmap <- function (fpkm_matrix,palette,legend_name,width) {
  p=Heatmap(fpkm_matrix,show_column_names = TRUE,show_row_names = FALSE, 
            row_names_gp =  gpar(fontsize = 2, fontface="italic", col="black"),
            #clustering_distance_columns = "pearson",
            clustering_distance_rows = "euclidean",
#            row_split = annotation_row$Regulation,
#            column_split = annotation_col$Group,
#            cluster_rows = opt$rowcluster,
            cluster_columns = F,
#            top_annotation = annotation_col_final,
#            left_annotation = annotation_row_final,
            #row_labels = "right",
            row_names_side = "right",
            width = unit(width, "cm"),
            column_names_rot = 45,
            #row
            column_title_rot = 0,
            row_title_rot = 0,
            column_dend_side = "top",
            row_dend_side = "left",
            column_title = legend_name, 
            #            column_title = "gene_id",
            col = palette,
            #column_title_side = "Top",
            row_title_side = "left",
            column_title_gp = gpar(fontsize=8, fontface="bold"), 
            row_title_gp = gpar(fontsize=1, fontface="bold"),
            heatmap_legend_param=list(title= legend_name,legend_direction="vertical",title_rot = 90)#labels_rot = 0
  )
  return (p)
}

palette <- colorRamp2(c(-4, 0, 4), c("RoyalBlue2", "White", "Red2"))
palette2 <- colorRamp2(c(0, 0.05), c("Black", "white"))

p1=getComplexHeatmap(Three,palette, "FC Color Key", 5)
p2=getComplexHeatmap(Four,palette2, "Sig", 3)
ht_list=p1+p2

pdf(paste0("./","diff_FC_SIG",".","output.pdf"),width=8,height=9, onefile = F)
draw(ht_list, 
     ht_gap = unit(1, "cm"))
dev.off()

png(paste0("./","diff_FC_SIG",".","output.png"),width=8,height=9)
draw(ht_list, 
     ht_gap = unit(1, "cm"))
dev.off()