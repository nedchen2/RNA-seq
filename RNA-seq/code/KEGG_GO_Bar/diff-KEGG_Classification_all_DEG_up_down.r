#!/usr/bin/env Rscript
library("optparse")
#library("ggthemes")

option_list = list(
	make_option(c("-i", "--input"), type="character", default=NULL, help="input file name", metavar="character"),
	make_option(c("-n", "--group"), type="character", default=NULL, help="group name, eg: Group2-vs-Group1", metavar="character"),
	make_option(c("-o", "--output"), type="character", default=NULL, help="output file prefix", metavar="character")
);
opt_parser = OptionParser(option_list=option_list, epilogue = "Rscript diff-KEGG_Classification_all_DEG_up_down.r -i KEGG_Classification_Up_Down.xls,KEGG_Classification_Up_Down.xls -n Group2-vs-Group1 -o out/KEGG_Classification_Up_Down");
opt = parse_args(opt_parser);
if (is.null(opt$input) | is.null(opt$output) | is.null(opt$group)){
	print_help(opt_parser)
	stop("--input --output --group must be supplied", call.=FALSE)
}

library(ggplot2)
data <- read.delim(opt$input, sep="\t", header=T, quote="")

uniq<-data$Classification_level2[!duplicated(data$Classification_level2)]
data$Classification_level2 <- factor(data$Classification_level2, levels=uniq)
mylabel=paste(data$Classification_level1, data$Classification_level2, sep="--")
#col=vector(length=length(data$Type))
#col[data$Type=="Down"]="red"
#col[data$Type=="Up"]="blue"
data$Type <- factor(data$Type, levels=c("48h","96h"));

cols=c("48h"="#377EB8", "96h"="#E41A1C")
p=ggplot(data=data, aes(x=Classification_level2, y=Percentage, width=0.8, fill=Type, space=0)) +
  coord_flip() +
  geom_bar(stat="identity",position="dodge") +
  geom_text(aes(label=Gene_number),position = position_dodge(width = 0.8), hjust=-0.5,size=3.5) + 
  scale_x_discrete(breaks=data$Classification_level2, labels=mylabel)+
  guides(fill=guide_legend(reverse=TRUE))+
  scale_fill_manual(values=cols)+
  xlab("")+ylab("Percent of Genes(%)")+
  labs(title=paste0("KEGG Pathway Classification", "(", as.character(opt$group), ")"),fill="Group") +
  theme(plot.title = element_text(hjust = 0.5, size=25,face = "bold")) +
  theme_classic() +
  theme(axis.text.y=element_text(size=11,color="black")) + 
  theme(axis.text.x=element_text(hjust=1, size=11,color="black")) +
  theme(legend.text=element_text(size=11))+
  theme(panel.grid =element_blank()) 

ggsave(paste0(opt$output, ".pdf"), height=10, width=14, plot=p)
ggsave(paste0(opt$output, ".png"), type="cairo-png", height=10, width=14, plot=p)
print(paste0(opt$output, ".png(pdf) is OK"));

