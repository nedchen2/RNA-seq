#!/usr/bin/env Rscript
# 用于输出三张go分类图
library("optparse")
option_list = list(
	make_option(c("-i", "--input"), type="character", default=NULL, help="input file name", metavar="character"),
	make_option(c("-m", "--mark"), type="character", default=NULL, help="select Up, Total or Down", metavar="character"),
	make_option(c("-o", "--outpath"), type="character", default=NULL, help="outfile directory", metavar="character")
);


opt_parser = OptionParser(option_list=option_list,epilogue = "Rscript top***.r -i enrichment-go-Group1-vs-Group2-Down.txt -m Down -o outdir/");
opt = parse_args(opt_parser);
if(is.null(opt$input) | is.null(opt$outpath) | is.null(opt$mark)){
	print_help(opt_parser)
	stop("--input --outpath --mark must be supplied", call.=FALSE)
}
if(!file.exists(opt$outpath)){dir.create(opt$outpath,recursive = T)}
opt$outpath<-gsub("/$", "", opt$outpath)

library(ggplot2)
library(stringr)
library(grid)
library(RColorBrewer)

groupname <- gsub("\\.(txt|xls)$", "", gsub("^enrichment-go-", "", basename(opt$input)))
if(grepl("-Down$", groupname, ignore.case = F, perl = F,fixed = F, useBytes = F)){
        groupname <- gsub("-Down$", "(Down)", groupname)
}
if(grepl("-Up$", groupname, ignore.case = F, perl = F,fixed = F, useBytes = F)){
        groupname <- gsub("-Up$", "(Up)", groupname)
}
if(grepl("-Total$", groupname, ignore.case = F, perl = F,fixed = F, useBytes = F)){
        groupname <- gsub("-Total$", "(Total)", groupname)
}

limit70 <- function(s) {
	k <- as.character(s)
	if(str_length(s)>70){k <- sub("[^ ]+$", "...", substr(k,1,67))}
	return(k)	
}

top10 <- function(i) { return(i[order(head(i, 10)["ListHits"]),]) }

d <- read.delim(opt$input, sep="\t", header=T, quote="")
d <- d[which(d[,"ListHits"]>2),]
d <- d[order(d[,"ListHits"],decreasing=T),]

#stopifnot(nrow(d)>0)
if(nrow(d)==0){
    print("d items = 0, program exit!")
    q()
}

dp <- rbind(top10(d[d["Category"]=="biological_process", ]),
	top10(d[d["Category"]=="cellular_component", ]),
	top10(d[d["Category"]=="molecular_function", ]))
write.table(dp[,c(1,2,3,4,8,10,11)], paste0(opt$outpath, "/GO.top.", opt$mark, ".xls"), sep="\t", quote=FALSE,
	col.names=TRUE, row.names=FALSE)

d1<-dp[which(dp[,3]=="biological_process"),]
d1_s<-d1[order(d1[,4],decreasing=T),]
d1_s["Term"] <- apply(d1_s["Term"], 1, limit70)
d1_s$Term <- factor(d1_s$Term, levels=d1_s$Term)

d2<-dp[which(dp[,3]=="cellular_component"),]
d2_s<-d2[order(d2[,4],decreasing=T),]
d2_s["Term"] <- apply(d2_s["Term"], 1, limit70)
d2_s$Term <- factor(d2_s$Term, levels=d2_s$Term)
d3<-dp[which(dp[,3]=="molecular_function"),]
d3_s<-d3[order(d3[,4],decreasing=T),]
d3_s["Term"] <- apply(d3_s["Term"], 1, limit70)
d3_s$Term <- factor(d3_s$Term, levels=d3_s$Term)

#d_l["Term"] <- apply(d_l["Term"], 1, limit70)
#d_l$Term <- factor(d_l$Term, levels=d_l$Term)
p=ggplot(data=d1_s, aes(x=Term, y=ListHits,width=0.6,space=0.6,cex.main=3,cex.lab=2)) +
  geom_bar(stat="identity",position=position_dodge(0.7),width=0.5,fill="#336699") +
  labs(x="", y="Number of gene", title = paste0(groupname,": ","(Biological process)","Top 10 GO Term")) +
  theme_bw() + 
  theme(axis.text.x=element_text(angle = 45, hjust=1, size=14))+
  theme(legend.position=c(-0.10,0.8), legend.key.width=unit(1, "lines"))+
  theme(legend.text=element_text(size=10))+
  theme(text=element_text(size=14)) +
  theme(plot.title = element_text(hjust = 0.5, vjust=4, size=18, family = "ArialMT"))+
  theme(plot.margin=unit(c(2,2,2,2), "lines"))+
  theme(panel.grid =element_blank()) 
ggsave(paste0(opt$outpath, "/GO(Biological_process).top.", "pdf"), height=10, width=18, plot=p)
ggsave(paste0(opt$outpath, "/GO(Biological_process).top.", "png"), type="cairo-png", height=10, width=18, plot=p)
print(paste0(opt$outpath, "/GO(Biological_process).top.", "png(pdf) is OK"));


#d_l<-rbind(d1,d2,d3)
#d_l[which(dp$Category=="biological_process"), "color"] <- "#4DAF4A"
#d_l[which(dp$Category=="cellular_component"), "color"] <- "#377EB8"
#d_l[which(dp$Category=="molecular_function"), "color"] <- "#E41A1C"


p2=ggplot(data=d2_s, aes(x=Term, y=ListHits,width=0.6,space=0.6,cex.main=3,cex.lab=2)) +
  geom_bar(stat="identity",position=position_dodge(0.7),width=0.5,fill="#336699") +
  labs(x="", y="Number of gene", title = paste0(groupname,": ","(Cellular component)","Top 10 GO Term")) +
  theme_bw() +   
  theme(axis.text.x=element_text(angle = 45, hjust=1, size=14))+
  theme(legend.position=c(-0.10,0.8), legend.key.width=unit(1, "lines"))+
  theme(legend.text=element_text(size=10))+
  theme(text=element_text(size=14)) +
  theme(plot.title = element_text(hjust = 0.5, vjust=4, size=18, family = "ArialMT"))+
  theme(plot.margin=unit(c(2,2,2,2), "lines"))+
  theme(panel.grid =element_blank()) 

ggsave(paste0(opt$outpath, "/GO(Cellular_component).top.", "pdf"), height=10, width=18, plot=p2)
ggsave(paste0(opt$outpath, "/GO(Cellular_component).top.", "png"), type="cairo-png", height=10, width=18, plot=p2)
print(paste0(opt$outpath, "/GO(Cellular_component).top.", "png(pdf) is OK"));

p3=ggplot(data=d3_s, aes(x=Term, y=ListHits,width=0.6,space=0.6,cex.main=3,cex.lab=2)) +
  geom_bar(stat="identity",position=position_dodge(0.7),width=0.5,fill="#336699") +
  labs(x="", y="Number of gene", title = paste0(groupname,": ","(Molecular function)","Top 10 GO Term")) +
  theme_bw() +
  theme(axis.text.x=element_text(angle = 45, hjust=1, size=14))+
  theme(legend.position=c(-0.10,0.8), legend.key.width=unit(1, "lines"))+
  theme(legend.text=element_text(size=10))+
  theme(text=element_text(size=14)) +
  theme(plot.title = element_text(hjust = 0.5, vjust=4, size=18, family = "ArialMT"))+
  theme(plot.margin=unit(c(2,2,2,2), "lines"))+
  theme(panel.grid =element_blank()) 

ggsave(paste0(opt$outpath, "/GO(Molecular_function).top.", "pdf"), height=10, width=18, plot=p3)
ggsave(paste0(opt$outpath, "/GO(Molecular_function).top.", "png"), type="cairo-png", height=10, width=18, plot=p3)
print(paste0(opt$outpath, "/GO(Molecular_function).top.", "png(pdf) is OK"));

