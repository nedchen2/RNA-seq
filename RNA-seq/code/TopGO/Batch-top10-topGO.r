#!/usr/bin/env Rscript
library("optparse")

#########安装###########
##if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")

#BiocManager::install("topGO")


option_list = list(
	make_option(c("-d", "--diff"), type="character", default=NULL, help="diff result file", metavar="character"),
	make_option(c("-m", "--mark"), type="character", default=NULL, help="select Up, Total or Down [optional]", metavar="character"),
	make_option(c("-j", "--gobg"), type="character", default=NULL, help="go backgroud file", metavar="character"),
	make_option(c("-o", "--outpath"), type="character", default=NULL, help="outfile directory", metavar="character")
);
opt_parser = OptionParser(option_list=option_list,epilogue = "Rscript topGO.r -d diff-Group2-vs-Group1-Total.xls -m Total -j go.backgroud.xls -o outdir/");
opt = parse_args(opt_parser);
if (is.null(opt$diff) | is.null(opt$gobg) | is.null(opt$outpath)){
	print_help(opt_parser)
	stop("--diff --gobg --outpath must be supplied", call.=FALSE)
}

groupname <- gsub("\\.(txt|xls)$", "", gsub(".split_result.xls$", "", basename(opt$diff)))

outputpath <- paste0(opt$outpath,"/",groupname,"/")

if(!dir.exists(outputpath)){dir.create(outputpath,recursive = T)}
#opt$outpath<-gsub("/$", "", opt$outpath)

library(topGO)
diff_genes<-read.delim(opt$diff, header=T, sep="\t", quote="")
geneID2GO<-readMappings(opt$gobg)
#在两个文件中都有的，被当作完整基因列表
interesting_genes=factor(diff_genes[,1])
geneNames<- names(geneID2GO)
geneList <- factor(as.integer (geneNames %in% interesting_genes))
names(geneList)=geneNames

#生成GO数据
GOdata <- new("topGOdata", ontology = "MF", allGenes = geneList,annot = annFUN.gene2GO, gene2GO = geneID2GO)
#进行富集分析
resultFisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
allRes <- GenTable(GOdata,classicFisher = resultFisher, orderBy = "classicFisher", topNodes = 20)
write.table(allRes,file=paste0(outputpath,"/","MF.topGO_RESULT.xls"),col.names=T,row.names=F,quote = FALSE,sep='\t')

anno_gene <-max(allRes[,4])
if (anno_gene>2){
	png(paste0(outputpath, "/topGO_MF", opt$mark, ".png"), width=3000, height=3000, res=300)
	showSigOfNodes(GOdata, score(resultFisher), firstSigNodes = 10, useInfo = 'all')
	dev.off()

	pdf(paste0(outputpath, "/topGO_MF", opt$mark, ".pdf"))
	showSigOfNodes(GOdata, score(resultFisher), firstSigNodes = 10, useInfo = 'all')
	dev.off()
	print(paste0(outputpath, "/topGO_MF", opt$mark, ".pdf(png) is OK"));
}

GOdata <- new("topGOdata", ontology = "BP", allGenes = geneList,annot = annFUN.gene2GO, gene2GO = geneID2GO)
resultFisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
#筛选并且出表格
allRes <- GenTable(GOdata,classicFisher = resultFisher, orderBy = "classicFisher", topNodes = 20)
write.table(allRes,file=paste0(outputpath,"/","BP.topGO_RESULT.xls"),row.names=T,quote = FALSE,sep='\t')


anno_gene <-max(allRes[,4])
if (anno_gene>2){
	png(paste0(outputpath, "/topGO_BP", opt$mark, ".png"), width=3000, height=3000, res=300)
	showSigOfNodes(GOdata, score(resultFisher), firstSigNodes = 10, useInfo = 'all')
	dev.off()

	pdf(paste0(outputpath, "/topGO_BP", opt$mark, ".pdf"))
	showSigOfNodes(GOdata, score(resultFisher), firstSigNodes = 10, useInfo = 'all')
	dev.off()
	print(paste0(outputpath, "/topGO_BP", opt$mark, ".pdf(png) is OK"));
}

GOdata <- new("topGOdata", ontology = "CC", allGenes = geneList,annot = annFUN.gene2GO, gene2GO = geneID2GO)
resultFisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
allRes <- GenTable(GOdata,classicFisher = resultFisher, orderBy = "classicFisher", topNodes = 20)
write.table(allRes,file=paste0(outputpath,"/","CC.topGO_RESULT.xls"),row.names=T,quote = FALSE,sep='\t')


anno_gene <-max(allRes[,4])
if (anno_gene>2){
	png(paste0(outputpath, "/topGO_CC", opt$mark, ".png"), width=3000, height=3000, res=300)
	showSigOfNodes(GOdata, score(resultFisher), firstSigNodes = 10, useInfo = 'all')
	dev.off()

	pdf(paste0(outputpath, "/topGO_CC", opt$mark, ".pdf"))
	showSigOfNodes(GOdata, score(resultFisher), firstSigNodes = 10, useInfo = 'all')
	dev.off()
	print(paste0(outputpath,"/topGO_CC", opt$mark, ".pdf(png) is OK"));
}
