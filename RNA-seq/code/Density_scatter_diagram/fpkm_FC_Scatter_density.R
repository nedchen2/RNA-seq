#这个脚本尚未完善，请注意


suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(gridExtra))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(viridis))


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

#如果要筛选显著，则需要这个代码
#标记显著性（默认 p < 0.05）
#two[which(two$FDR_A < 0.05 & two$FDR_B < 0.05),'type1'] <- 'sign'
#two[which(two$FDR_A >= 0.05 | two$FDR_B >= 0.05),'type1'] <- 'no'

#Three <- subset(two,type1=='sign',select=c("logFC_A","logFC_B"))

#标记差异倍数（默认 |log2FC| >= 1）
#two[which(two$logFC_A <= -1 & two$logFC_B <= -1),'type2'] <- 'a_down.b_down'
#two[which(two$logFC_A >= 1 & two$logFC_B <= -1),'type2'] <- 'a_up.b_down'
#two[which(two$logFC_A <= -1 & two$logFC_B >= 1),'type2'] <- 'a_down.b_up'
#two[which(two$logFC_A >= 1 & two$logFC_B >= 1),'type2'] <- 'a_up_b_up'
#two[is.na(two$type2),'type2'] <- 'no'

#rownames(two)=two[,1]

#showname_q<-paste0("Sig.","(",groupname,")")
#showname_q2<-paste0("Sig.","(",groupname2,")")

#two = plyr::rename(two, c("FDR_A"=showname_q))
#two = plyr::rename(two, c("FDR_B"=showname_q2))

showname_FC<-paste0("log2FC.","(",groupname,")")
showname_FC2<-paste0("log2FC.","(",groupname2,")")
two = plyr::rename(two, c("logFC_A"=showname_FC))
two = plyr::rename(two, c("logFC_B"=showname_FC2))

rownames(two)=two[,1]

Three <- subset(two,type1=='sign',select=c(showname_FC,showname_FC2))

p <- ggplot(Three, aes(x=logFC_A, y=logFC_B) )  + geom_bin2d(bins = 150) +
  scale_fill_continuous(type = "viridis") + theme_bw() +   geom_vline(xintercept = c(0, 0), lty = 2) +
  geom_hline(yintercept = c(0, 0), lty = 2) + xlab(showname_FC) + 
  ylab(showname_FC2)
p

a<-cor.test(Three$logFC_A,Three$logFC_B,method = "spearman")

cor.test(Three$logFC_A,Three$logFC_B,method = "spearman",exact = FALSE)
l <- list(a = as.numeric(format(a$estimate, digits = 4)),
          b = as.numeric(0.05))

if (a$p.value < 0.05){
  eq <- substitute(italic(Cor)~"="~a~","~italic(P)~"<"~b, l)
}else{
  eq <- substitute(italic(Cor)~"="~a~","~italic(P)~">"~b, l)
}

p1 <- p + labs(x = paste0('log2fpkm ',"(",groupname,")"), y = paste0('log2fpkm ',"(",groupname2,')'),subtitle = as.expression(eq))
p1

ggsave(paste0("./",'DensityScatter-FC2-log2fpkm.pdf'), p1, width = 7, height = 7)
ggsave(paste0("./",'DensityScatter-FC2-log2fpkm.png'), p1, width = 7, height = 7)
print(paste0( "Scatter", ".png(pdf) is OK"))
        



