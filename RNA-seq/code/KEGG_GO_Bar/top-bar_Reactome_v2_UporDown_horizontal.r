#!/usr/bin/env Rscript
# 平行条形图
library("optparse")
option_list = list(
	make_option(c("-i", "--input"), type="character", default=NULL, help="input file name", metavar="character"),
	make_option(c("-m", "--mark"), type="character", default=NULL, help="select Up, Total or Down", metavar="character"),
	make_option(c("-o", "--outpath"), type="character", default=NULL, help="outfile directory", metavar="character"),
	make_option(c("-n", "--Number"), type="character", default=10, help="The number of terms you want to conserve", metavar="character")
);
opt_parser = OptionParser(option_list=option_list,epilogue = "Rscript top10X3_GO.r -i enrichment-kegg-Group1-vs-Group2-Down.txt -m Down -o outdir/  Rscript top10X3_kegg.r -i enrichment-kegg-H-vs-C-Fake.xls -m Total -o outdir/
");
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

groupname <- gsub("\\.(txt|xls)$", "", gsub("^enrichment-kegg-", "", basename(opt$input)))
if(grepl("-Down$", groupname, ignore.case = F, perl = F,fixed = F, useBytes = F)){
        groupname <- gsub("-Down$", "(Down)", groupname)
}
if(grepl("-Up$", groupname, ignore.case = F, perl = F,fixed = F, useBytes = F)){
        groupname <- gsub("-Up$", "(Up)", groupname)
}
if(grepl("-Total$", groupname, ignore.case = F, perl = F,fixed = F, useBytes = F)){
        groupname <- gsub("-Total$", "(Total)", groupname)
}

Number <-as.numeric(opt$Number)

print ("We will do top ")
print (Number)

limit70 <- function(s) {
	k <- as.character(s)
	if(str_length(s)>70){k <- sub("[^ ]+$", "...", substr(k,1,67))}
	return(k)	
}

top10 <- function(i) { return(i[order(head(i, Number)["ListHits"]),]) }

d <- read.delim(opt$input, sep="\t", header=T, quote="")
d <- d[which(d[,"ListHits"]>0),]

d_l <- head(d[order(d$ListHits,decreasing=T),],Number)

d_l <- top10(d_l)
#stopifnot(nrow(d)>0)
if(nrow(d)==0){
    print("d items = 0, program exit!")
    q()
}
write.table(d_l[,c(1,2,3,4,8)], paste0(opt$outpath, "/Reactome.top.", opt$mark, ".xls"), sep="\t", quote=FALSE,
	col.names=TRUE, row.names=FALSE)


#d_l["Term"] <- apply(d_l["Term"], 1, limit70)
d_l$Term <- factor(d_l$Term, levels=d_l$Term)

p=ggplot(data=d_l, aes(x=Term, y=ListHits, width=0.6, fill=ListHits,space=0.6,cex.main=3,cex.lab=2)) +

  geom_bar(stat="identity",position=position_dodge(0.7),width=0.5) +
  labs(x="", y='ProteinCounts', fill='ProteinCounts', title = paste0(groupname,": ","Top ",opt$Number," Reactome Term")) +
  theme_bw() + scale_colour_gradient(low = "green", high = "red") +
  theme(axis.text.x=element_text(hjust=1, size=14,color=d_l$color))+
  theme(legend.position=c(-0.10,0.8), legend.key.width=unit(1, "lines"))+
  theme(legend.text=element_text(size=10))+
  theme(text=element_text(size=14)) +
  theme(plot.title = element_text(hjust = 0.5, vjust=4, size=18, family = "ArialMT"))+
  theme(plot.margin=unit(c(2,10,2,10), "lines"))+
  theme(panel.grid =element_blank())+
  coord_flip()+
  theme(legend.position="right")


ggsave(paste0(opt$outpath, "/Reactome.top.", opt$mark, ".pdf"), height=6, width=18, plot=p)
ggsave(paste0(opt$outpath, "/Reactome.top.", opt$mark, ".png"), type="cairo-png", height=6, width=18, plot=p)
print(paste0(opt$outpath, "/Reactome.top.", opt$mark, ".png(pdf) is OK"));
