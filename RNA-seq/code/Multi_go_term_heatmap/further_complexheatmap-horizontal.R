#==========import library==========
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(gridExtra))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(viridis))
suppressPackageStartupMessages(library(ComplexHeatmap))
suppressPackageStartupMessages(library(circlize))
suppressPackageStartupMessages(library(pheatmap))


#############################################
usage = "\
usage:
    *X基因在所有样本中数值相等将报错 , 需调整输入文件！！*
    /home/hanmin/anaconda3/envs/mro3.5/bin/Rscript further_heatmap.R  -e [fpkm or data matrix file]  -d [phenodata matrix file] -o ./Heatmap
input grouping file format:
    *sample_group.xls*:
    Sample      Group
    Sample_A    group2
    Sample_B    group2
input phenotype file format:
    *anno.xls*:
    gene_id    Celltype
    gene1      Endothelial
    gene2      Fibroblasts
"
cat(usage)

#==========parameter import==========
suppressPackageStartupMessages(library(optparse))
option_list = list(
    make_option(c("-e", "--expression"), type = "character", default = NULL,
        help = "Expression matrix file name(force). "),
    make_option(c("-d", "--group"), type = "character", default = NULL,
        help = "Group information of samplenames.  e.g. sample_group.xls. " ),
    make_option(c("-t", "--title"), type = "character", default = NULL,
        help = "Graphic title information: Group_A-vs-Group_B or \"Group A vs Group B\". "),
    make_option(c("-s", "--showgenes"), type = "logical", default = "F",
        help = "Whether to show gene id, [default: F] . " ),
    make_option(c("-r", "--rowcluster"), type = "logical", default = "T",
        help = "Whether rows (genes) are clustered : T or F , [default: T] . "),
    make_option(c("-l", "--colcluster"), type = "logical", default = "T",
        help = "Whether the columns (samples) are clustered : T or F , [default: T] . "),
    make_option(c("-a", "--angle"), type = "double",  default = 90,
        help = "angle of the column labels, right now one can choose only from few predefined : 0, 45, 90, 270 and 315, [default 90]", metavar = "double"),
    make_option(c("-f", "--fontsize"), type = "double",default = NULL, 
        help = "Base fontsize for the plot, [default %default]", metavar = "double"),
    make_option(c("-p", "--phenotype"), type = "character", default = NULL,
        help = "Phenotype information of genes.  e.g. phenotype.xls. " ),
    make_option(c("-c", "--colors"), type = "character", default = "redwhiteblue",
        help = "colors choise for Heatmap picture : redwhiteblue, redblackgreen ,yellowblackblue , [default: redwhiteblue]. "),
    make_option(c("-g", "--height"), type = "double",  default = NULL,
        help = "Height limit [default %default]", metavar = "double"),
    make_option(c("-k", "--width"), type = "double",  default = NULL,
        help = "Width limit [default %default]", metavar = "double"),
    make_option(c("-o", "--outputdir"), type = "character", default = "./Heatmap",
        help = "output directory for Heatmap results ,[ default: ./Heatmap] . "),
	make_option(c("-u", "--cutrow"), type = "double", default = "1",
        help = "How many sections do you want? ,[ default:1] . ")
);
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

#==========parameter==========

#############读取表达文件###########
if (! is.null(opt$expression)) {
    expression <- read.table(normalizePath(opt$expression), header=T, sep="\t", check.names=F, quote="")
}else {
    print_help(opt_parser)
    stop("expression matrix file is not supplied, heatmap for replicate groups  will be skipped\n", call. = FALSE)
}

#############组名注释##############
if ( !is.null(opt$group) ){
    annotation_col <- read.delim(normalizePath(opt$group), header=T, sep="\t", row.names=1,check.names=F, quote="")
}else {
    annotation_col = NA
}


##############基因注释################
if ( !is.null(opt$phenotype) ){
    annotation_row <- read.delim(normalizePath(opt$phenotype), header=T, sep="\t", row.names=1,check.names=F, quote="")
}else {
    annotation_row = NA
}

##############标题#############
if ( !is.null(opt$title) ){
    title = opt$title
}else {
    title = "Heatmap"
}

##############角度#################
angle_list = unlist( strsplit( "0,45,90,270,315", ",", perl = T) )
if ( !opt$angle %in% angle_list ){stop("You can choose only from few predefined : 0, 45, 90, 270 and 315")}
if ( !is.null(opt$angle) ){
    angle = opt$angle
}else {
    angle = 90
}

##########颜色################
if ( !is.null(opt$colors) ){
    if ( opt$colors == "redwhiteblue" ){
        palette <- colorRampPalette(c("RoyalBlue2", "White", "Red2"))(n=256)
    }else if ( opt$colors == "redblackgreen" ){
        palette <- colorRampPalette(c("Green", "Black", "Red"))(n=256)
    }else if ( opt$colors == "yellowblackblue" ){
        palette <- colorRampPalette(c("Blue", "Black", "Yellow"))(n=256)
    }
}else {
    palette <- colorRampPalette(c("RoyalBlue2", "White", "Red2"))(n=256)
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

#==========Annotation================

groupname <- gsub("\\.(txt|xls)$", "", gsub("(-annotation.xls)$", "", basename(opt$expression)))
if(grepl("-Down$", groupname, ignore.case = F, perl = F,fixed = F, useBytes = F)){
  groupname <- gsub("-Down$", "(Down)", groupname)
}
if(grepl("-Up$", groupname, ignore.case = F, perl = F,fixed = F, useBytes = F)){
  groupname <- gsub("-Up$", "(Up)", groupname)
}
if(grepl("-Total$", groupname, ignore.case = F, perl = F,fixed = F, useBytes = F)){
  groupname <- gsub("-Total$", "(Total)", groupname)
}

#生成不重复的term列表
#term_list <- expression$term[!duplicated(expression$term)]

term_list <- levels(expression$term)

#term_list <- unique(expression$term)
#palette <- colorRampPalette(viridis(3))(n=256)
palette <- colorRamp2(c(-1, 0, 1), c("navy","white","firebrick3"))
#palette <- colorRamp2(c(-2, 0, 2), c(viridis(3)))

########定义画图脚本##########
limit30 <- function(s) {
  k <- as.character(s)
  if(str_length(s)>30){k <- sub("[^ ]+$", "...", substr(k,1,27))}
  return(k)	
}

GetComplexHeatmap <- function (fpkm_matrix,term_name) {
  p=Heatmap(fpkm_matrix,show_column_names = FALSE,show_row_names = TRUE,name = "Color key",
            #clustering_distance_columns = "pearson",
            #clustering_distance_rows = "euclidean",
            cluster_rows = opt$rowcluster,
            cluster_columns = T,
            column_title_rot = 90,
            row_title_rot = 270,
            column_dend_side = "bottom",
            #column_title = limit30(term_name), 
            column_title = term_name,
            col = palette,
            #column_title_side = "Top",
            row_title_side = "right",
            column_title_gp = gpar(fontsize=8, fontface="bold"), 
            row_title_gp = gpar(fontsize=1, fontface="bold"),
            heatmap_legend_param=list(title= "Color Key", legend_direction="horizontal",labels_rot = 90,title_rot = 90)
            )
  return (p)
}

GetClusterMatrix <- function (fpkm,pic) {
  if ( (opt$rowcluster == T) & (opt$colcluster == T ) ){
    met_reorder = fpkm[row_order(pic),column_order(pic)]
    
  }else if( (opt$rowcluster == F) & (opt$colcluster == T ) ) {
    met_reorder = fpkm[,column_order(pic)]
  }
  return (met_reorder)
}

##########生成该条目的表达量文件##########
#在热图列表中始终有一个主热图，用于控制全局的行顺序，其他所有的热图会根据主热图的配置自动进行调整，调整方式为：

#不对行进行聚类，按照主热图的行顺序排列#
#删除行标题
#如果主热图进行了分割，也会进行分割
#主热图的高度就是所有热图的高度


##############根据总表进行拆分#############
for (i in 1:length(term_list)) {
  print (paste0("我们现在正在处理",term_list[i]))
####处理表达量文件
  tmp=expression[expression["term"]==term_list[i],][,-2]
  row.names(tmp)<-tmp$gene_id
  tmp=scale(t(tmp[,-1]))
  assign(paste("df_fpkm", i, sep = ""), tmp)#批量创建数据框并赋值  
  assign(paste("df_term", i, sep = ""), term_list[i])
  assign(paste("pic", i, sep = ""), GetComplexHeatmap(tmp,
                                                        term_list[i]))
  #批量画图add_heatmap(pic1,pic2)
  #tmp <- get(paste("df", i, sep="")) #得到我们需要的数据
  print ("================================================")
}

for (i in 1:length(term_list)) {
  tmp=get(paste("df_fpkm", i, sep=""))
  tmp_pic=get(paste("pic", i, sep=""))
  assign(paste("mat",i,sep = ""),GetClusterMatrix(tmp,tmp_pic)) 
}

#mat1
#print ("================================================")
#df_fpkm1[row_order(pic1),column_order(pic1)]
#print ("================================================")
#df_fpkm1[row_order(get(paste("pic",1,sep=""))),column_order(get(paste("pic", 1, sep="")))]

#####################出图######################################
pdf(paste0(output_dir,"/",groupname,".",opt$rowcluster,"-output.pdf"),width=11,height=8)
a=get(paste("pic",1,sep = ""))
for (i in 2:length(term_list)) {
  a=a+get(paste("pic",i,sep = "")) #本质上就是add_Heatmap
}
draw(a,main_heatmap = 1,
     #ht_gap=unit(c(3, 1), "mm"),
    gap=unit(1,"mm"), padding = unit(c(2, 10, 8, 2), "mm"),#col_dend_side = "bott",
    merge_legend = TRUE, heatmap_legend_side = "top",
     row_title = "Heatmaps of Gene Expression for Multiple GO Term", row_title_gp = gpar(fontsize = 16, fontface="bold"))
dev.off()


#########################出第一张图########################

##################输出聚类之后的表格########

b=get(paste("mat",1,sep = ""))
tmp_df=cbind(gene_id=row.names(b), b)
row.names(tmp_df)=NULL

for (i in 2:length(term_list)) {
  m2=get(paste("mat",i,sep = ""))
  m3=cbind(gene_id=row.names(m2), m2)
  row.names(m3)=NULL
  tmp_df=cbind(tmp_df,m3)
}


write.table(tmp_df,file=paste0(output_dir,"/",groupname,".",opt$rowcluster,".genes_reorder.cluster_result.xls"),row.names=F,quote = FALSE,sep='\t')
write.table(t(tmp_df),file=paste0(output_dir,"/",groupname,".",opt$rowcluster,".genes_reorder.after_transformation.xls"),row.names=T,quote = FALSE,col.names =F,sep='\t')





##########################还可以用split来拆分#######################
#列子
#Heatmap(df, name = "mtcars", col = mycol, split = mtcars$cyl )


#ComplexHeatmap
###可以在顶部展示注释信息和注释图
#添加点图
#annotation = data.frame(value = rnorm(10))
#value = 1:10
#ha = HeatmapAnnotation(df = annotation, points = anno_points(value), 
#                       annotation_height = c(1, 2))
# 添加barplot注释信息
#ha = HeatmapAnnotation(barplot = anno_barplot(1:12, which = "row", bar_width=0.4, gp= gpar(fill="red")), which = "row")
#add_heatmap(ht1, ha) ####可以添加注释也可以添加热图
# 添加boxplot注释信息
#ha = HeatmapAnnotation(boxplot = anno_boxplot(matrix(rnorm(60), nrow=12), which = "row", border = F, gp= gpar(fill="blue")), which = "row")
#add_heatmap(ht2, ha)
# 添加histogram注释信息
#ha = HeatmapAnnotation(histogram = anno_histogram(matrix(rnorm(48), nrow=12), which = "row", gp= gpar(fill="red")), which = "row")
#add_heatmap(ht2, ha)
# 添加density注释信息
#ha = HeatmapAnnotation(density = anno_density(matrix(rnorm(48), nrow=12), which = "row", type="heatmap"), which = "row")
#add_heatmap(ht2, ha)
#ht1=Heatmap(expression,heatmap_legend_param=list(title= "Color Key", title_position = "topcenter", legend_height=unit(8,"cm"), legend_direction="vertical"),column_title = "ComplexHeatmap", column_title_gp = gpar(fontsize = 20, fontface = "bold", #col="red"),row_dend_side = "right",show_row_names = TRUE,show_column_names = FALSE, show_row_dend = TRUE,km = 2, row_title_gp = gpar(col=rainbow(4)), row_names_gp = gpar(col="red", fontsize=20))#top_annotation=


#cell_fun = NULL, #cell_fun：自定义在cell中增加绘图项的函数。7个参数：i(row index,矩阵中的行index）, j(column index，矩阵中的列index), x,y(热图体区中中间点的坐标）,width,height(cell的宽度和高度）,fill(cell的填充颜色） 

