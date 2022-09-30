#==========import library==========
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(tidyr))

#############################################
usage = "\
 根据某一列，拆分表达量表格
"
cat(usage)

#==========parameter import==========
suppressPackageStartupMessages(library(optparse))
option_list = list(
    make_option(c("-e", "--expression"), type = "character", default = NULL,
        help = "Expression matrix file name(force). "),
    make_option(c("-o", "--outputdir"), type = "character", default = "./Result",
        help = "output directory for split results ,[ default: ./Result] . "),
    make_option(c("-c", "--col_name"), type = "character", default = "Cluster",
        help = "the col_name which you want to split,[ default: Cluster] . ")
);
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

#==========parameter==========

#############读取表达文件###########
if (! is.null(opt$expression)) {
    expression <- read.delim(normalizePath(opt$expression), header=T, sep="\t", check.names=F, quote="")
    
}else {
    print_help(opt_parser)
    stop("expression matrix file is not supplied, heatmap for replicate groups  will be skipped\n", call. = FALSE)
}

#######路径##########
if ( is.null(opt$outputdir) ){
    output_dir = "Heatmap"
}else{
    if ( file.exists(opt$outputdir) ){
        output_dir = opt$outputdir
    }else{
        output_dir = opt$outputdir
        dir.create(output_dir)
    }
}

col <- expression[1:nrow(expression),opt$col_name]
print ("本次的目标列：")
print (head(col))
if ( is.numeric(col) ){
print ("数值型")
expression <- mutate(expression,
  Cluster_Name = paste0("cluster", expression$Cluster)) #神一样的函数
term_list <- unique(expression$Cluster_Name)
#term_list <- levels(expression["Cluster_Name"])
for (i in 1:length(term_list)) {
  print (paste0("我们现在正在处理",term_list[i]))
####处理表达量文件
  tmp=expression[expression["Cluster_Name"]==term_list[i],]
  row.names(tmp)<-tmp$gene_id
  write.table(tmp,file=paste0(output_dir,"/",term_list[i],".","split_result.txt"),row.names=F,quote = FALSE,sep='\t')
  print ("================================================")
   }
}else{
print ("非数值型")
term_list <- expression[,opt$col_name][!duplicated(expression[,opt$col_name])]
print (term_list)
##############根据总表进行拆分#############
for (i in 1:length(term_list)) {
  print (paste0("我们现在正在处理",term_list[i]))
####处理表达量文件
  tmp=expression[expression[,opt$col_name]==term_list[i],]
  row.names(tmp)<-tmp$gene_id
  groupname1 <- gsub("\\.(txt|xls)$", "", gsub("(-diff-pval-0.05-FC-2.mRNA.xls)$", "", basename(opt$expression)))
  groupname <- gsub("/", "_", gsub(" ", "_", term_list[i])) #去除特殊符号
  groupname <- str_replace_all(groupname, "[^[:alnum:]]", "_")
  write.table(expression,file=paste0(output_dir,"/",groupname1,".","Total",".split_result.xls"),row.names=F,quote = FALSE,sep='\t')
  write.table(tmp,file=paste0(output_dir,"/",groupname1,".",groupname,".split_result.xls"),row.names=F,quote = FALSE,sep='\t')
  print ("================================================")
  }
}