#meta-bishe

#准备itol注释文件####
# Tong Zhou, Kuidong Xu, Feng Zhao et. al. itol.toolkit accelerates working with 
# iTOL (Interactive Tree Of Life) by an automated generation of annotation files, 
# Bioinformatics, 2023;, btad339, https://doi.org/10.1093/bioinformatics/btad339

#参考学习https://tongzhou2017.github.io/itol.toolkit/articles/Get_Start.html
# load package
library(itol.toolkit)
library(data.table) # 用于元数据文件读取
library(ape) # 用于树文件读写
library(stringr) # 用于字符串处理
setwd('D:/wenjian/毕业设计/大论文/数据/meta/itol')
tree <- "bin.unrooted.tree" # 树分枝名存在多余信息会导致注释文件匹配失败
# 修改树分枝名
phylo <- ape::read.tree(tree)
#phylo$tip.label <- stringr::str_remove(phylo$tip.label, "_.*$") # 去除多余信息
#tree <- "tree_short_label.nwk"
#ape::write.tree(phylo,tree) # 导出树，实际上传该树进行可视化


phylo_depth <- vcv.phylo(phylo)
phylo_depth_diag <- diag(phylo_depth) %>% as.data.frame() %>% setNames('phylo_depth')


# 创建本次分析的数据及样式仓库
hub <- create_hub(tree)

# 加载元数据
# data_file <- "iTOL_annotation.txt"
# data <- data.table::fread(data_file)
# data <- data %>% filter(ID != "Seq30_") # 去除多余信息

#读取物种信息
setwd("D:/wenjian/毕业设计/大论文/数据/meta/MAGs")
library(readxl)
tax <- read.delim('gtdbtk.bac120.summary.tsv')
tax <- tax[,1:2]
max_cols <- max(sapply(strsplit(tax$classification, ";"), length))# 确定最大的拆分列数
col_names <- paste0("col", 1:max_cols)# 生成列名
library(tidyr)
taxonomy <- tax %>%
  separate(classification, into = col_names, sep = ";", fill = "right")# 拆分字符串到列

#读取丰度信息
mags <- read.delim('bins_abundance_tpm.tsv')

#合并物种与丰度
mags <- merge(mags,taxonomy,by.x = 'Genome',by.y = 'user_genome')

setwd('D:/wenjian/毕业设计/大论文/数据/meta/itol')
#添加属信息（分枝重命名）
library(dplyr)
# df_genera <- mags %>% select(Genome, col6)
# unit_2 <- create_unit(data = df_genera, 
#                       key = "rep_Li2022jhm_4_text_1", 
#                       type = "DATASET_TEXT", 
#                       size_factor = 1,
#                       rotation= 180,
#                       position = -1,
#                       color = "#000000",
#                       tree = tree)
#unit_2@common_themes$basic_theme$margin <- -500
df_phycolor <- mags %>% select(Genome, col2)
unit_1 <- create_unit(data = df_phycolor, 
                      key = "tree_color", 
                      type = "TREE_COLORS", 
                      subtype = "c", 
                     # line_type = c(rep("normal",4),"dashed"),
                     # size_factor = 5, 
                      tree = tree)
write_unit(unit_1)

#分枝重命名
df_genera <- mags %>% select(Genome, col6)
unit_3 <- create_unit(data = df_genera,
                      key = "genera_labels",
                      type = "LABELS",
                      tree = tree)
write_unit(unit_3)
#按门设置颜色
df_phycolor <- mags %>% select(Genome, col2)
set.seed(666)
unit_4 <- create_unit(data = df_phycolor, 
                      key = "rep_Li2022jhm_3_range", 
                      type = "TREE_COLORS", 
                      subtype = "range", 
                      color = "wesanderson",
                      tree = tree)
#绘制基因丰度热图
genes <- read.delim('D:/wenjian/毕业设计/大论文/数据/meta/MAGs/tongji2.txt')

unit_5 <- create_unit(data = genes,
                      key = 'genes_abundance_heatmap',
                      type = 'DATASET_HEATMAP',
                      tree = tree
                      )
write_unit(unit_5)

#绘制质粒丰度热图
plasmid <- read.delim('D:/wenjian/毕业设计/大论文/数据/meta/itol/result.txt')

unit_5 <- create_unit(data = plasmid,
                      key = 'plasmid_abundance',
                      type = 'DATASET_MULTIBAR',
                      tree = tree
)
write_unit(unit_5)

#绘制MAGs丰度
df_magabun <- mags[,-17:-23]
unit_30 <- create_unit(data = df_magabun,
                       key = "mag_abundance", 
                       type = "DATASET_MULTIBAR", 
                       tree = tree)
write_unit(unit_30)

#添加MAGs的大小
#读取数据
tax_mag <- read.delim('D:/wenjian/毕业设计/大论文/数据/meta/MAGs/gtdbtk.bac120.summary.tsv')[,1:2]
com_cont <- read.csv('D:/wenjian/毕业设计/大论文/数据/meta/MAGs/genomeInformation.csv')
#tax_inf <- merge(tax_mag,com_cont,by.x = 'user_genome',by.y = 'genome')
com_cont$genome <- gsub(".fa", "", com_cont$genome)
tax_inf <- com_cont[com_cont$genome %in% tax_mag$user_genome,]

df_maglenth <- tax_inf %>% select(genome,length)
unit_10 <- create_unit(data = df_maglenth,
                       key = "mag_lenth", 
                       type = "DATASET_MULTIBAR", 
                       tree = tree)
write_unit(unit_10)


hub <- hub + unit_2 + unit_4
write_hub(hub, getwd())


#MAGs的种类及丰度分析####
setwd('D:/wenjian/毕业设计/大论文/数据/meta/MAGs')
library(readxl)
tax <- read.delim('gtdbtk.bac120.summary.tsv')
tax <- tax[,1:2]
max_cols <- max(sapply(strsplit(tax$classification, ";"), length))# 确定最大的拆分列数
col_names <- paste0("col", 1:max_cols)# 生成列名
library(tidyr)
taxonomy <- tax %>%
  separate(classification, into = col_names, sep = ";", fill = "right")# 拆分字符串到列

#读取丰度信息
mags <- read.delim('bins_abundance_tpm.tsv')
#合并物种与丰度
mags <- merge(mags,taxonomy,by.x = 'Genome',by.y = 'user_genome')
men <- select(mags,-Genome,-col1,-col3,-col4,-col5,-col6,-col7)
men <- men %>% 
  group_by(col2) %>% 
  summarise_all(sum) %>% 
  as.data.frame() #将重复门合并

#reads层面ARG种类及丰度分析####
#读取数据
setwd("D:/wenjian/毕业设计/大论文/数据/meta/ARG")
type_arg <- read.delim('normalized_cell.type.txt',row.names = 1)
colnames(type_arg) <- c('CG11','CG21','Cd11',
                        'Cd21','Cr11','Cr21',
                        'Cu11','Cu21','PAP11',
                        'PAP21','PNP11','PNP21',
                        'Seed','oh11','oh21')
type_arg10 <- type_arg[c('bacitracin','multidrug','sulfonamide',
                       'aminoglycoside','macrolide-lincosamide-streptogramin',
                       'beta_lactam','polymyxin'),]
type_arg10$type <- rownames(type_arg10)


library(reshape2)
type_long <- melt(type_arg10)
type_long <- type_long[order(as.character(type_long$variable)),]
type_long$category <- c(rep('heavy_metal',times=14),
                        rep('no',times=14),
                        rep('heavy_metal',times=28),
                        rep('organic',times=42),
                        rep('no',times=7)
                        )
library(dplyr)
type_long <- type_long

type_long2 <- type_long[type_long$variable=='Cd11'|
                         type_long$variable=='Cu11'|
                         type_long$variable=='Cr11'|
                          type_long$variable=='oh11'|
                          type_long$variable=='PAP11'|
                          type_long$variable=='PNP11'
                         # type_long$variable=='CG11'|
                         # type_long$variable=='Seed'
                        ,]
type_long2$category <-factor(type_long2$category,ordered=TRUE,
                             levels=c("heavy_metal","organic")) #修改因子水平 

type_long3 <- type_long[type_long$variable=='Cd21'|
                          type_long$variable=='Cu21'|
                          type_long$variable=='Cr21'|
                          type_long$variable=='oh21'|
                          type_long$variable=='PAP21'|
                          type_long$variable=='PNP21'
                         # type_long$variable=='CG21'|
                         # type_long$variable=='Seed'
                        ,]
type_long3$category <-factor(type_long3$category,ordered=TRUE,
                             levels=c("heavy_metal","organic","no")) #修改因子水平 

# #有些type丰度太高，需要挑出来画
library(ggplot2)
ggplot(type_long3,aes(x=type,y=value,fill=category))+
  geom_boxplot(
    position=position_dodge(0.9)                 ##因为分组比较，需设组间距
  )+
  stat_boxplot(mapping=aes(x = reorder(type, -value, FUN = median),y=value),
               geom ="errorbar",      ##添加箱子的bar为最大、小值
               width=0.3,
               position=position_dodge(0.9)   ###这里的0.9要与上面的0.9对应
               )+     ##bar宽度和组间距
  scale_y_continuous(limits=c(0,0.3))+      ##修改y轴的范围
  scale_fill_manual(values=c("cyan3","indianred1","#8f888b"))+ #设置分组的颜色
  labs(x="",y="Copies of ARG per cell (60 d)")+        ##修改坐标轴和图例的文本
  theme_classic()+
  theme(axis.text.x = element_text(angle = 15,vjust = 0.85,hjust = 0.75))+ #横轴旋转
  theme(panel.background=element_rect(fill='transparent'),   ##去掉底层阴影
        panel.grid=element_blank(),                          ##去掉网格线
        #panel.border=element_blank(),                        ##去掉图的边界线
        panel.border=element_rect(fill=NA,        ##设置图的边界线
                                  color='black',
                                  size = 0.5        ##设置边框粗细
                                  ),
        axis.line=element_line(colour="black"),     ##设置坐标轴的线条
        axis.title=element_text(face="bold",size = 13),      ##设置坐标轴文本字体
        # axis.title.x=element_blank(),                        ##删除x轴文本               
        axis.text=element_text(color = 'black',size = 10),     ##设置坐标轴标签字体
        #axis.ticks=element_line(color='black'),              ##设置坐标轴刻度线
        # legend.position=c(0.15, 0.9),                        ##设置图例的位置
        # legend.direction = "horizontal",                      ##设置图列标签的排列方式
        #legend.title=element_text(face = "bold",size=12),    ##设置图例标题字体
        legend.title=element_blank(),                        ##去掉图例标题         
        # legend.text=element_text(face = "bold",size=12),     ##设置图例文本的字体
        legend.background=element_rect(linetype="solid",     ##设置图列的背景和线条颜色
                                       colour ="black"))
##手动显著性检验
try <- type_long
try$value <- try$value
kruskal.test(value~category,data=try)

#arg的α多样性####
setwd("D:/wenjian/毕业设计/大论文/数据/meta/ARG")
type_subarg <- read.delim('normalized_cell.subtype.txt',row.names = 1)
colnames(type_subarg) <- c('CG11','CG21','Cd11',
                           'Cd21','Cr11','Cr21',
                           'Cu11','Cu21','PAP11',
                           'PAP21','PNP11','PNP21',
                           'Seed','oh11','oh21')

#计算α多样性指数并绘制小提琴图
library(vegan)
richness <- colSums(type_subarg > 0)
shannon <- diversity(t(type_subarg),"shannon")
simpson <- t(type_subarg) %>% diversity("simpson")
fisher <- t(type_subarg) %>% diversity("invsimpson")
diversityindex <- data.frame(SampleID = rownames(t(type_subarg)),
                             ShannonIndex = shannon,
                             simpson = simpson,
                             fisher = fisher,
                             richness = richness)
#读取metadata
md <- read.delim('D:/wenjian/毕业设计/大论文/数据/meta/metadata.txt')
#合并
md <- merge(md,diversityindex,by.x = 'sample',by.y = 'SampleID')
#调整因子顺序，以便绘图时是想要的顺序
md$rongyu <- factor(md$rongyu,levels = c("y","n"))
md$category <- factor(md$category,levels = c("heavy_metal","organic","no"))
#计算显著性差异
library(tidyverse)
library(ggh4x)
library(rstatix)
library(ggpubr)
md_p_val1 <- md %>% group_by(category)%>%
  wilcox_test(richness ~rongyu) %>%
  adjust_pvalue(p.col="p",method="bonferroni") %>%
  add_significance(p.col="p.adj") %>% 
  add_xy_position(x="rongyu",dodge=0.8)

#绘图
md[-15,] %>% ggplot(aes(rongyu,ShannonIndex))+
  # annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = 1,
  #          fill = "cornflowerblue", alpha = .3, color = NA)+
  # annotate("rect", xmin = -Inf, xmax = Inf, ymin = 2, ymax = Inf,
  #          fill = "#FCAE12", alpha = .3, color = NA) +
  geom_violin(aes(fill=category),trim = FALSE,show.legend = F)+
  geom_boxplot(width = 0.2,outliers = FALSE, staplewidth = 0.5) +
  scale_y_continuous(sec.axis = sec_axis(~ ., name = ""))+
  stat_summary(fun=mean,geom="point",col="black",fill="#F98400",
               shape=23,show.legend = F)+
  stat_summary(fun=mean, geom="line", aes(group=category), col="black")+
  stat_pvalue_manual(md_p_val1,label="p.adj.signif",hide.ns=T,
                     tip.length = 0,label.size = 5,color="black")+
  # geom_hline(yintercept = 2,linetype=2)+
  # geom_hline(yintercept = 1,linetype=2)+
  geom_vline(xintercept = 2.6,linetype=2)+
  facet_nested_wrap(. ~ category,
                    strip = strip_nested(background_x = 
                                           elem_list_rect(fill=c("cyan3","indianred1","#8f888b")))) +
  scale_fill_manual(values = c("cyan3","indianred1","#8f888b"))+
  labs(x=NULL,y='Copies of ARG per cell')+
  theme(axis.text.x=element_text(angle = 0,vjust=0.5,hjust=0.5,color="black"),
        axis.text.y=element_text(color="black"),
        plot.background = element_rect(fill="white"), 
        panel.background = element_rect(fill="white"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.spacing.x = unit(0,"cm"),
        plot.margin=unit(c(0.5,0.5,0.5,0.5),unit="cm"),
        axis.line.x=element_line(color="black"),
        axis.line.y.left = element_line(color="grey30"),
        axis.line.y.right = element_line(color="grey30"),
        axis.line.x.bottom = element_line(color="grey30"),
        axis.text.y.right = element_blank(),
        axis.ticks.y.right = element_blank())+
  guides(y = guide_axis(minor.ticks = TRUE))

#arg的β多样性，nmds####
#calculating stress values for NMDS,判断是否能用nmds分析。stress<=0.2,stressplot图中点基本都在线附近
dat <- t(type_subarg)
md <- read.delim('D:/wenjian/毕业设计/大论文/数据/meta/metadata.txt')
md$sample <- order(md$sample,levels = c("CG11","CG21","Cd11","Cd21","Cr11",
                                        "Cr21","Cu11","Cu21","PAP11","PAP21",
                                        "PNP11","PNP21","Seed","oh11","oh21"))

vare.dis<-vegdist(dat)
library(MASS)
vare.mds0<-isoMDS(vare.dis)
stressplot(vare.mds0,vare.dis)
## to calculate ordination
ord<-metaMDS(dat, trace=FALSE)
#saveRDS(ord, file = file.path(result_folder,"ord_v6.rds")) # to save the R object 'ord'

# to draw the ordination result
#ord <- readRDS(file.path(result_folder,"ord_v6.rds")) # to read in the R object 'ord'
par(mar = c(4, 4, 1, 2)) # mar = c(bottom, left, top, right)
pl<-plot(ord,dis="sites",type = 't') # type = 't' for plotting site-labels
col<-as.vector(md$categorycol) #color are respond with location
#pch<-as.vector(dat.graph$pch)#pch点的形状对应是否加了污染物
points(ord,display="sites",cex=2,col=col,
       pch=21, 
       bg=col)#设置点的大小、形状、颜色等
abline(v=0,h=0,lty=2,col='grey')#增加0,0线

#nmds的显著性检验
##stress
stress <- ord[["stress"]]
stressplot(ord)
##分别进行Anosim分析（Analysis of similarities）和PERMANOVA（即adonis）检验分析。
set.seed(123)#设置随机种子
###基于bray-curtis距离进行PERMANOVA分析(即adonis);
adonis <-  adonis2(dat ~ category, md, permutations = 999, method = "bray")
adonis$R2
###基于bray-curtis距离进行anosim分析
anosim = anosim(dat, md$category, permutations = 999, distance = "bray")
R <- anosim$statistic
p <- anosim$signif
plot(anosim)



#reads层面MRG种类及丰度分析####
#读取数据
setwd("D:/wenjian/毕业设计/大论文/数据/meta/MRG")
type_mrg <- read.delim('normalized_cell.level2.txt',row.names = 1)
colnames(type_mrg) <- c('CG11','CG21','Cd11',
                        'Cd21','Cr11','Cr21',
                        'Cu11','Cu21','PAP11',
                        'PAP21','PNP11','PNP21',
                        'Seed','oh11','oh21')
# type_arg10 <- type_arg[c('bacitracin','multidrug','sulfonamide',
#                          'aminoglycoside','macrolide-lincosamide-streptogramin',
#                          'beta_lactam','polymyxin'),]
 type_mrg$type <- rownames(type_mrg)


library(reshape2)
type_long <- melt(type_mrg)
type_long <- type_long[order(as.character(type_long$variable)),]
type_long$category <- c(rep('heavy_metal',times=32),
                        rep('no',times=32),
                        rep('heavy_metal',times=64),
                        rep('organic',times=96),
                        rep('no',times=16)
)
library(dplyr)
type_long <- type_long

type_long2 <- type_long[type_long$variable=='Cd11'|
                          type_long$variable=='Cu11'|
                          type_long$variable=='Cr11'|
                          type_long$variable=='oh11'|
                          type_long$variable=='PAP11'|
                          type_long$variable=='PNP11'|
                          type_long$variable=='CG11'|
                          type_long$variable=='Seed',]
type_long2$category <-factor(type_long2$category,ordered=TRUE,
                             levels=c("heavy_metal","organic","no")) #修改因子水平 

type_long3 <- type_long[type_long$variable=='Cd21'|
                          type_long$variable=='Cu21'|
                          type_long$variable=='Cr21'|
                          type_long$variable=='oh21'|
                          type_long$variable=='PAP21'|
                          type_long$variable=='PNP21'|
                          type_long$variable=='CG21'|
                          type_long$variable=='Seed',]
type_long3$category <-factor(type_long3$category,ordered=TRUE,
                             levels=c("heavy_metal","organic","no")) #修改因子水平 

# #有些type丰度太高，需要挑出来画
library(ggplot2)
ggplot(type_long3,aes(x=type,y=value,fill=category))+
  geom_boxplot(
    position=position_dodge(0.9)                 ##因为分组比较，需设组间距
  )+
  stat_boxplot(mapping=aes(x = reorder(type, -value, FUN = median),y=value),
               geom ="errorbar",      ##添加箱子的bar为最大、小值
               width=0.3,
               position=position_dodge(0.9)   ###这里的0.9要与上面的0.9对应
  )+     ##bar宽度和组间距
#  scale_y_continuous(limits=c(0,0.3))+      ##修改y轴的范围
  scale_fill_manual(values=c("cyan3","indianred1","#8f888b"))+ #设置分组的颜色
  labs(x="",y="Copies of MRG per cell")+        ##修改坐标轴和图例的文本
  theme_classic()+
  theme(axis.text.x = element_text(angle = 15,vjust = 0.85,hjust = 0.75))+ #横轴旋转
  theme(panel.background=element_rect(fill='transparent'),   ##去掉底层阴影
        panel.grid=element_blank(),                          ##去掉网格线
        #panel.border=element_blank(),                        ##去掉图的边界线
        panel.border=element_rect(fill=NA,        ##设置图的边界线
                                  color='black',
                                  size = 1        ##设置边框粗细
        ),
        #axis.line=element_line(colour="black"),     ##设置坐标轴的线条
        axis.title=element_text(face="bold",size = 13),      ##设置坐标轴文本字体
        # axis.title.x=element_blank(),                        ##删除x轴文本               
        axis.text=element_text(color = 'black',size = 10),     ##设置坐标轴标签字体
        #axis.ticks=element_line(color='black'),              ##设置坐标轴刻度线
        # legend.position=c(0.15, 0.9),                        ##设置图例的位置
        # legend.direction = "horizontal",                      ##设置图列标签的排列方式
        #legend.title=element_text(face = "bold",size=12),    ##设置图例标题字体
        legend.title=element_blank(),                        ##去掉图例标题         
        # legend.text=element_text(face = "bold",size=12),     ##设置图例文本的字体
        legend.background=element_rect(linetype="solid",     ##设置图列的背景和线条颜色
                                       colour ="black"))
##手动显著性检验
try <- type_long
try$value <- try$value
kruskal.test(value~category,data=try)




#reads层面aromadeg种类及丰度分析####
#读取数据
setwd("D:/wenjian/毕业设计/大论文/数据/meta/aromadeg")
type_aromadeg <- read.delim('normalized_cell.level4.txt',row.names = 1)
colnames(type_aromadeg) <- c('CG11','CG21','Cd11',
                        'Cd21','Cr11','Cr21',
                        'Cu11','Cu21','PAP11',
                        'PAP21','PNP11','PNP21',
                        'Seed','oh11','oh21')
# type_arg10 <- type_arg[c('bacitracin','multidrug','sulfonamide',
#                          'aminoglycoside','macrolide-lincosamide-streptogramin',
#                          'beta_lactam','polymyxin'),]
type_aromadeg$type <- rownames(type_aromadeg)


library(reshape2)
type_long <- melt(type_aromadeg)
type_long <- type_long[order(as.character(type_long$variable)),]
type_long$category <- c(rep('heavy_metal',times=88),
                        rep('no',times=88),
                        rep('heavy_metal',times=176),
                        rep('organic',times=264),
                        rep('no',times=44)
)
library(dplyr)
type_long <- type_long

type_long2 <- type_long[type_long$variable=='Cd11'|
                          type_long$variable=='Cu11'|
                          type_long$variable=='Cr11'|
                          type_long$variable=='oh11'|
                          type_long$variable=='PAP11'|
                          type_long$variable=='PNP11'|
                          type_long$variable=='CG11'|
                          type_long$variable=='Seed',]
type_long2$category <-factor(type_long2$category,ordered=TRUE,
                             levels=c("heavy_metal","organic","no")) #修改因子水平 

type_long3 <- type_long[type_long$variable=='Cd21'|
                          type_long$variable=='Cu21'|
                          type_long$variable=='Cr21'|
                          type_long$variable=='oh21'|
                          type_long$variable=='PAP21'|
                          type_long$variable=='PNP21'|
                          type_long$variable=='CG21'|
                          type_long$variable=='Seed',]
type_long3$category <-factor(type_long3$category,ordered=TRUE,
                             levels=c("heavy_metal","organic","no")) #修改因子水平 

#有些type丰度差异较大，需要挑出来
top15 <- type_aromadeg$type[order(rowSums(type_aromadeg[,-16]),decreasing = 1)][1:15]
type_long22 <- type_long2 %>% subset(type %in% top15)
type_long33 <- type_long3 %>% subset(type %in% top15)
type_long222 <- type_long2 %>% subset(type=='Non')
type_long333 <- type_long3 %>% subset(type=='Non')

library(ggplot2)
ggplot(type_long333,aes(x=type,y=value,fill=category))+
  geom_boxplot(
    position=position_dodge(0.9)                 ##因为分组比较，需设组间距
  )+
  stat_boxplot(mapping=aes(x = reorder(type, -value, FUN = median),y=value),
               geom ="errorbar",      ##添加箱子的bar为最大、小值
               width=0.3,
               position=position_dodge(0.9)   ###这里的0.9要与上面的0.9对应
  )+     ##bar宽度和组间距
  #  scale_y_continuous(limits=c(0,0.1))+      ##修改y轴的范围
  scale_fill_manual(values=c("cyan3","indianred1","#8f888b"))+ #设置分组的颜色
  labs(x="",y="Copies of AromaDeg per cell")+        ##修改坐标轴和图例的文本
  theme_classic()+
  theme(axis.text.x = element_text(angle = 15,vjust = 0.85,hjust = 0.75))+ #横轴旋转
  theme(panel.background=element_rect(fill='transparent'),   ##去掉底层阴影
        panel.grid=element_blank(),                          ##去掉网格线
        #panel.border=element_blank(),                        ##去掉图的边界线
        panel.border=element_rect(fill=NA,        ##设置图的边界线
                                  color='black',
                                  size = 1        ##设置边框粗细
        ),
        #axis.line=element_line(colour="black"),     ##设置坐标轴的线条
        axis.title=element_text(face="bold",size = 13),      ##设置坐标轴文本字体
        # axis.title.x=element_blank(),                        ##删除x轴文本               
        axis.text=element_text(color = 'black',size = 10),     ##设置坐标轴标签字体
        #axis.ticks=element_line(color='black'),              ##设置坐标轴刻度线
        # legend.position=c(0.15, 0.9),                        ##设置图例的位置
        # legend.direction = "horizontal",                      ##设置图列标签的排列方式
        #legend.title=element_text(face = "bold",size=12),    ##设置图例标题字体
        legend.title=element_blank(),                        ##去掉图例标题   
        legend.position = "none",                           ##删除图例
        # legend.text=element_text(face = "bold",size=12),     ##设置图例文本的字体
        legend.background=element_rect(linetype="solid",     ##设置图列的背景和线条颜色
                                       colour ="black"))
##手动显著性检验
try <- type_long
try$value <- try$value
kruskal.test(value~category,data=try)




#reads层面MGE种类及丰度分析####
#读取数据
setwd("D:/wenjian/毕业设计/大论文/数据/meta/MGE")
type_mge <- read.delim('normalized_cell.level2.txt',row.names = 1)
colnames(type_mge) <- c('CG11','CG21','Cd11',
                             'Cd21','Cr11','Cr21',
                             'Cu11','Cu21','PAP11',
                             'PAP21','PNP11','PNP21',
                             'Seed','oh11','oh21')
# type_arg10 <- type_arg[c('bacitracin','multidrug','sulfonamide',
#                          'aminoglycoside','macrolide-lincosamide-streptogramin',
#                          'beta_lactam','polymyxin'),]
type_mge$type <- rownames(type_mge)


library(reshape2)
type_long <- melt(type_mge)
type_long <- type_long[order(as.character(type_long$variable)),]
type_long$category <- c(rep('heavy_metal',times=10),
                        rep('no',times=10),
                        rep('heavy_metal',times=20),
                        rep('organic',times=30),
                        rep('no',times=5)
)
library(dplyr)
type_long <- type_long
type_long$category <-factor(type_long$category,ordered=TRUE,
                             levels=c("heavy_metal","organic","no")) #修改因子水平 


type_long2 <- type_long[type_long$variable=='Cd11'|
                          type_long$variable=='Cu11'|
                          type_long$variable=='Cr11'|
                          type_long$variable=='oh11'|
                          type_long$variable=='PAP11'|
                          type_long$variable=='PNP11'|
                          type_long$variable=='CG11'|
                          type_long$variable=='Seed',]
type_long2$category <-factor(type_long2$category,ordered=TRUE,
                             levels=c("heavy_metal","organic","no")) #修改因子水平 

type_long3 <- type_long[type_long$variable=='Cd21'|
                          type_long$variable=='Cu21'|
                          type_long$variable=='Cr21'|
                          type_long$variable=='oh21'|
                          type_long$variable=='PAP21'|
                          type_long$variable=='PNP21'|
                          type_long$variable=='CG21'|
                          type_long$variable=='Seed',]
type_long3$category <-factor(type_long3$category,ordered=TRUE,
                             levels=c("heavy_metal","organic","no")) #修改因子水平 

#有些type丰度差异较大，需要挑出来
top15 <- type_aromadeg$type[order(rowSums(type_aromadeg[,-16]),decreasing = 1)][1:15]
type_long22 <- type_long2 %>% subset(type %in% top15)
type_long33 <- type_long3 %>% subset(type %in% top15)
type_long222 <- type_long2 %>% subset(type=='Non')
type_long333 <- type_long3 %>% subset(type=='Non')

library(ggplot2)
ggplot(type_long3,aes(x=type,y=value,fill=category))+
  geom_boxplot(
    position=position_dodge(0.9)                 ##因为分组比较，需设组间距
  )+
  stat_boxplot(mapping=aes(x = reorder(type, -value, FUN = median),y=value),
               geom ="errorbar",      ##添加箱子的bar为最大、小值
               width=0.3,
               position=position_dodge(0.9)   ###这里的0.9要与上面的0.9对应
  )+     ##bar宽度和组间距
  #  scale_y_continuous(limits=c(0,0.1))+      ##修改y轴的范围
  scale_fill_manual(values=c("cyan3","indianred1","#8f888b"))+ #设置分组的颜色
  labs(x="",y="Copies of MGE per cell")+        ##修改坐标轴和图例的文本
  theme_classic()+
  theme(axis.text.x = element_text(angle = 15,vjust = 0.85,hjust = 0.75))+ #横轴旋转
  theme(panel.background=element_rect(fill='transparent'),   ##去掉底层阴影
        panel.grid=element_blank(),                          ##去掉网格线
        #panel.border=element_blank(),                        ##去掉图的边界线
        panel.border=element_rect(fill=NA,        ##设置图的边界线
                                  color='black',
                                  size = 1        ##设置边框粗细
        ),
        #axis.line=element_line(colour="black"),     ##设置坐标轴的线条
        axis.title=element_text(face="bold",size = 13),      ##设置坐标轴文本字体
        # axis.title.x=element_blank(),                        ##删除x轴文本               
        axis.text=element_text(color = 'black',size = 10),     ##设置坐标轴标签字体
        #axis.ticks=element_line(color='black'),              ##设置坐标轴刻度线
        # legend.position=c(0.15, 0.9),                        ##设置图例的位置
        # legend.direction = "horizontal",                      ##设置图列标签的排列方式
        #legend.title=element_text(face = "bold",size=12),    ##设置图例标题字体
        legend.title=element_blank(),                        ##去掉图例标题   
        legend.position = "none",                           ##删除图例
        # legend.text=element_text(face = "bold",size=12),     ##设置图例文本的字体
        legend.background=element_rect(linetype="solid",     ##设置图列的背景和线条颜色
                                       colour ="black"))
##手动显著性检验
try <- type_long
try$value <- try$value
kruskal.test(value~category,data=try)




#对args总丰度绘制小提琴图####
library(tidyverse)
library(ggh4x)
library(rstatix)
library(ggpubr)
#读取args数据
setwd("D:/wenjian/毕业设计/大论文/数据/meta/ARG")
type_arg <- read.delim('normalized_cell.type.txt',row.names = 1)
type_arg <- read.delim('rpkm.type.txt',row.names = 1)
colnames(type_arg) <- c('CG11','CG21','Cd11',
                        'Cd21','Cr11','Cr21',
                        'Cu11','Cu21','PAP11',
                        'PAP21','PNP11','PNP21',
                        'Seed','oh11','oh21')
#读取metadata
md <- read.delim('D:/wenjian/毕业设计/大论文/数据/meta/metadata.txt')
library(dplyr)
md <- md %>% arrange(sample)
#把arg总丰度合并到metadata
abunarg <- type_arg %>% summarise_all(sum) %>% t() %>% as.data.frame()
colnames(abunarg) <- c('argabun')
md <- bind_cols(md,abunarg)
#调整因子顺序，以便绘图时是想要的顺序
md$rongyu <- factor(md$rongyu,levels = c("y","n"))
md$category <- factor(md$category,levels = c("heavy_metal","organic","no"))

md <- read.delim('D:/wenjian/毕业设计/大论文/数据/meta/metadata2.txt')
#计算显著性差异
md_p_val1 <- md %>% group_by(category)%>%
  wilcox_test(argabun ~rongyu) %>%
  adjust_pvalue(p.col="p",method="bonferroni") %>%
  add_significance(p.col="p.adj") %>% 
  add_xy_position(x="rongyu",dodge=0.8)

# 单样本t检验,判断一列数字是否显著高于某个值。
t_test_result <- t.test(md[c(8,10,15),]$argabun, mu = 0.3390, alternative = "greater")

#绘图
md[,] %>% ggplot(aes(rongyu,argabun))+
  # annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = 1,
  #          fill = "cornflowerblue", alpha = .3, color = NA)+
  # annotate("rect", xmin = -Inf, xmax = Inf, ymin = 2, ymax = Inf,
  #          fill = "#FCAE12", alpha = .3, color = NA) +
  geom_violin(aes(fill=category),trim = FALSE,show.legend = F)+
  geom_boxplot(width = 0.2,outliers = FALSE, staplewidth = 0.5) +
  scale_y_continuous(sec.axis = sec_axis(~ ., name = ""))+
  stat_summary(fun=mean,geom="point",col="black",fill="#F98400",
               shape=23,show.legend = F)+
  stat_summary(fun=mean, geom="line", aes(group=category), col="black")+
  stat_pvalue_manual(md_p_val1,label="p.adj.signif",hide.ns=T,
                     tip.length = 0,label.size = 5,color="black")+
  # geom_hline(yintercept = 2,linetype=2)+
  # geom_hline(yintercept = 1,linetype=2)+
  geom_vline(xintercept = 2.6,linetype=2)+
  facet_nested_wrap(. ~ category,
                    strip = strip_nested(background_x = 
                                           elem_list_rect(fill=c("cyan3","indianred1","#8f888b")))) +
  scale_fill_manual(values = c("cyan3","indianred1","#8f888b"))+
  labs(x=NULL,y='Copies of ARG per cell')+
  theme(axis.text.x=element_text(angle = 0,vjust=0.5,hjust=0.5,color="black"),
        axis.text.y=element_text(color="black"),
        plot.background = element_rect(fill="white"), 
        panel.background = element_rect(fill="white"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.spacing.x = unit(0,"cm"),
        plot.margin=unit(c(0.5,0.5,0.5,0.5),unit="cm"),
        axis.line.x=element_line(color="black"),
        axis.line.y.left = element_line(color="grey30"),
        axis.line.y.right = element_line(color="grey30"),
        axis.line.x.bottom = element_line(color="grey30"),
        axis.text.y.right = element_blank(),
        axis.ticks.y.right = element_blank())+
  guides(y = guide_axis(minor.ticks = TRUE))

#对其他基因绘制小提琴图####
library(tidyverse)
library(ggh4x) #facet_nested_wrap
library(rstatix)
library(ggpubr) #stat_pvalue_manual
#读取args数据
setwd("D:/wenjian/毕业设计/大论文/数据/meta/aromadeg")
type_arg <- read.delim('normalized_cell.level4.txt',row.names = 1)
colnames(type_arg) <- c('CG11','CG21','Cd11',
                        'Cd21','Cr11','Cr21',
                        'Cu11','Cu21','PAP11',
                        'PAP21','PNP11','PNP21',
                        'Seed','oh11','oh21')
#读取metadata
md <- read.delim('D:/wenjian/毕业设计/大论文/数据/meta/metadata.txt')
library(dplyr)
md <- md %>% arrange(sample)
#把arg总丰度合并到metadata
abunarg <- type_arg %>% summarise_all(sum) %>% t() %>% as.data.frame()
colnames(abunarg) <- c('argabun')
md <- bind_cols(md,abunarg)
#调整因子顺序，以便绘图时是想要的顺序
md$rongyu <- factor(md$rongyu,levels = c("y","n"))
md$category <- factor(md$category,levels = c("heavy_metal","organic","no"))
#计算显著性差异
md_p_val1 <- md %>% group_by(category)%>%
  wilcox_test(argabun ~rongyu) %>%
  adjust_pvalue(p.col="p",method="bonferroni") %>%
  add_significance(p.col="p.adj") %>% 
  add_xy_position(x="rongyu",dodge=0.8)

md_p_val2 <- md %>% group_by(rongyu)%>%
  wilcox_test(argabun ~category) %>%
  adjust_pvalue(p.col="p",method="bonferroni") %>%
  add_significance(p.col="p.adj") %>% 
  add_xy_position(x="category",dodge=0.8)
#绘图
md[c(-13,-1,-2),] %>% ggplot(aes(rongyu,argabun))+
  # annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = 1.94598,
  #          fill = "cornflowerblue", alpha = .3, color = NA)+
  # annotate("rect", xmin = -Inf, xmax = Inf, ymin = 2, ymax = Inf,
  #          fill = "#FCAE12", alpha = .3, color = NA) +
  geom_violin(aes(fill=category),trim = FALSE,show.legend = F)+
  geom_boxplot(width = 0.2,outliers = FALSE, staplewidth = 0.5) +
  scale_y_continuous(sec.axis = sec_axis(~ ., name = ""))+
  stat_summary(fun=mean,geom="point",col="black",fill="#F98400",
               shape=23,show.legend = F)+
  stat_summary(fun=mean, geom="line", aes(group=category), col="black")+
  stat_pvalue_manual(md_p_val1,label="p.adj.signif",hide.ns=T,
                     tip.length = 0,label.size = 5,color="black")+
  # geom_hline(yintercept = 2,linetype=2)+
  geom_hline(yintercept = 0.795,linetype=2)+
  geom_vline(xintercept = 2.6,linetype=2)+
  facet_nested_wrap(. ~ category,
                    strip = strip_nested(background_x = 
                                           elem_list_rect(fill=c("cyan3","indianred1","#8f888b")))) +
  scale_fill_manual(values = c("cyan3","indianred1","#8f888b"))+
  labs(x=NULL,y='Copies of ADG per cell')+
  theme(axis.text.x=element_text(angle = 0,vjust=0.5,hjust=0.5,color="black"),
        axis.text.y=element_text(color="black"),
        plot.background = element_rect(fill="white"), 
        panel.background = element_rect(fill="white"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.spacing.x = unit(0,"cm"),
        plot.margin=unit(c(0.5,0.5,0.5,0.5),unit="cm"),
        axis.line.x=element_line(color="black"),
        axis.line.y.left = element_line(color="grey30"),
        axis.line.y.right = element_line(color="grey30"),
        axis.line.x.bottom = element_line(color="grey30"),
        axis.text.y.right = element_blank(),
        axis.ticks.y.right = element_blank())+
  guides(y = guide_axis(minor.ticks = TRUE))

# 单样本t检验,判断一列数字是否显著高于某个值。
t_test_result <- t.test(md[c(3,5,7),]$argabun, mu = 1.945, alternative = "greater")








#对ARG subtypes绘制venn图####
library(ggvenn)
library(tidyverse)
#
setwd("D:/wenjian/毕业设计/大论文/数据/meta/ARG")
type_subarg <- read.delim('normalized_cell.subtype.txt',row.names = 1)
colnames(type_subarg) <- c('CG11','CG21','Cd11',
                        'Cd21','Cr11','Cr21',
                        'Cu11','Cu21','PAP11',
                        'PAP21','PNP11','PNP21',
                        'Seed','oh11','oh21')
ss <- t(type_subarg)##转置表

metal30 <- colSums(ss[c(3,5,7),])
organic30 <- colSums(ss[c(9,11,14),])
CK30 <- ss[c(1),]

metal60 <- colSums(ss[c(4,6,8),])
organic60 <- colSums(ss[c(10,12,15),])
CK60 <- ss[c(2),]
seed <- ss[c(13),]
a <- names(ss[1,])[metal30>0] #如果第一行的某一列值大于0，则将这一列的列名赋值给a
b <- names(ss[1,])[organic30>0]#如果第二行的某一列值大于0，则将这一列的列名赋值给a
c <- names(ss[1,])[CK30>0]

d <- names(ss[1,])[metal60>0]
e <- names(ss[1,])[organic60>0]
f <- names(ss[1,])[CK60>0]
s <- names(ss[1,])[seed>0]

patterns <- list(
  radialGradient(c("white","cyan3"), group = F),
  radialGradient(c("white","indianred1"), group = FALSE),
  radialGradient(c("white","#8f888b"), group = FALSE)
  #radialGradient(c("white","#FF8C00"), group = FALSE)
  )

list(metal30=a,organnic30=b) %>% 
  ggvenn(show_percentage = T,show_elements = F,label_sep = ",",
         digits = 1,stroke_color =NA)+
  scale_fill_manual(values = patterns)

list(metal60=d,organnic60=e) %>% 
  ggvenn(show_percentage = T,show_elements = F,label_sep = ",",
         digits = 1,stroke_color =NA)+
  scale_fill_manual(values = patterns)




#对其他基因绘制venn图####
library(ggvenn)
library(tidyverse)
#
setwd("D:/wenjian/毕业设计/大论文/数据/meta/aromadeg")
type_subarg <- read.delim('normalized_cell.level1.txt',row.names = 1)
colnames(type_subarg) <- c('CG11','CG21','Cd11',
                           'Cd21','Cr11','Cr21',
                           'Cu11','Cu21','PAP11',
                           'PAP21','PNP11','PNP21',
                           'Seed','oh11','oh21')
ss <- t(type_subarg)##转置表

metal30 <- colSums(ss[c(3,5,7),])
organic30 <- colSums(ss[c(9,11,14),])
CK30 <- ss[c(1),]

metal60 <- colSums(ss[c(4,6,8),])
organic60 <- colSums(ss[c(10,12,15),])
CK60 <- ss[c(2),]
seed <- ss[c(13),]
a <- names(ss[1,])[metal30>0] #如果第一行的某一列值大于0，则将这一列的列名赋值给a
b <- names(ss[1,])[organic30>0]#如果第二行的某一列值大于0，则将这一列的列名赋值给a
c <- names(ss[1,])[CK30>0]

d <- names(ss[1,])[metal60>0]
e <- names(ss[1,])[organic60>0]
f <- names(ss[1,])[CK60>0]
s <- names(ss[1,])[seed>0]

patterns <- list(
  radialGradient(c("white","cyan3"), group = F),
  radialGradient(c("white","indianred1"), group = FALSE),
  radialGradient(c("white","#8f888b"), group = FALSE),
  radialGradient(c("white","#FF8C00"), group = FALSE)
)

list(metal30=a,organnic30=b) %>% 
  ggvenn(show_percentage = T,show_elements = F,label_sep = ",",
         digits = 1,stroke_color =NA)+
  scale_fill_manual(values = patterns)

list(metal60=d,organnic60=e) %>% 
  ggvenn(show_percentage = T,show_elements = F,label_sep = ",",
         digits = 1,stroke_color =NA)+
  scale_fill_manual(values = patterns)

#MAGs层面对各个MAG中基因存在情况的统计####
#过滤条件：e值小于1e-5
 setwd('D:/wenjian/毕业设计/大论文/数据/meta/MAGs')
# tongjigenes <- read.delim('tongji.txt',row.names = 1)
# tj1_guiyi <- log2(tongjigenes + 1)
# library(pheatmap)
# pheatmap(tj1_guiyi)

#过滤条件：-e 1e-5 --id 80 --query-cover 80,
#amino acid identity ≥ 80% and a query coverage ≥ 80%
tongjigenes2 <- read.delim('tongji2.txt',row.names = 1)
#归一化数据
tj2_guiyi <- log2(tongjigenes2 + 1)
# 自定义 Nature 风格的颜色向量
#nature_colors <- c("#FFFFFF", "#FF0000", "#FFA500", "#FFFF00", "#008000", "#0000FF")

# 筛选所有列之和不为0的行
genes_filtered <- tongjigenes2[rowSums(tongjigenes2) != 0, ]
genesfil_guiyi <- log2(genes_filtered + 1)

# 自定义渐变色向量
library(pheatmap)
custom_colors <- colorRampPalette(c("lightblue", "yellow", "red"))(100)
pheatmap(tj2_guiyi,
         cluster_cols = FALSE,  #取消对列的聚类
         show_rownames=F,
         color = custom_colors,
         border_color = 'lightgrey'
         )


#MAGs丰度分析####
#读取丰度信息
setwd('D:/wenjian/毕业设计/大论文/数据/meta/MAGs')
abun_mag <- read.delim('bins_abundance_tpm.tsv')
colnames(abun_mag) <- c('Genome','Cd11','Cd21','CG11',
                        'CG21','Cr11','Cr21',
                        'Cu11','Cu21','oh11',
                        'oh21','PAP11','PAP21',
                        'PNP11','PNP21','Seed')

#读取物种信息
tax_mag <- read.delim('gtdbtk.bac120.summary.tsv')[,1:2]
library(tidyr)
library(dplyr)
max_cols <- max(sapply(strsplit(tax_mag$classification, ";"), length))# 确定最大的拆分列数
col_names <- paste0("col", 1:max_cols)# 生成列名
taxonomy <- tax_mag %>%
  separate(classification, into = col_names, sep = ";", fill = "right")# 拆分字符串到列

#合并物种与丰度
mags <- merge(abun_mag,taxonomy,by.x = 'Genome',by.y = 'user_genome')

#属级别
mags_shu <- mags[,c(2:16,22)]

#绘图
library(reshape2)
library(ggplot2)
mags_shu_melt <- melt(mags_shu)

# library(RColorBrewer)
# # 使用颜色调色板自动分配颜色
# num_groups <- 257  # 计算组别数量
# colors <- brewer.pal(num_groups, "Set1")  # 使用 Set1 调色板

ggplot(mags_shu_melt,aes(x=variable,y=value,fill=col6))+
  geom_bar(stat = 'identity',position = 'stack')+
  theme(
    legend.position = "none" #删除图例
  )


#MAGs完整度、污染度、基因组大小分析####

setwd('D:/wenjian/毕业设计/大论文/数据/meta/MAGs')
#读取数据
tax_mag <- read.delim('gtdbtk.bac120.summary.tsv')[,1:2]
com_cont <- read.csv('genomeInformation.csv')
#tax_inf <- merge(tax_mag,com_cont,by.x = 'user_genome',by.y = 'genome')
com_cont$genome <- gsub(".fa", "", com_cont$genome)
tax_inf <- com_cont[com_cont$genome %in% tax_mag$user_genome,]
sum(tax_inf$completeness>95)

#完整度大于95的MAG的基因组大小
genome_95_5 <- tax_inf[tax_inf$completeness>95,]

hist(genome_95_5$length, main = "Histogram of MAGs length", xlab = "length", ylab = "Frequency", 
     col = "skyblue", border = "black")


#高质量基因组中携带ARG...基因与不携带那些基因的MAGs的基因组大小比较####
#组合信息，得到高质量MAGs中的基因携带情况
genome_95_5_genes <- tongjigenes2[rownames(tongjigenes2) %in% genome_95_5$genome,]
genome_95_5_genes$genome <- rownames(genome_95_5_genes)
genome_95_5 <- merge(genome_95_5_genes,genome_95_5,by.x = 'genome',by.y = 'genome')
genome_95_5$genesum <- apply(genome_95_5[,2:6],1,sum)

#
hist(genome_95_5[genome_95_5$genesum==0,]$length, main = "Histogram of MAGs length", xlab = "length", ylab = "Frequency", 
     col = "skyblue", border = "black")
hist(genome_95_5[genome_95_5$genesum!=0,]$length, main = "Histogram of MAGs length", xlab = "length", ylab = "Frequency", 
     col = "skyblue", border = "black")

#对otu和几个基因表相互之间的mantel检验####
###arg与otu
arggenes <- read.delim('D:/wenjian/毕业设计/大论文/数据/meta/ARG/normalized_cell.subtype.txt',
                       row.names = 1
                       )
colnames(arggenes) <- c('CG11','CG21','Cd11',
                           'Cd21','Cr11','Cr21',
                           'Cu11','Cu21','PAP11',
                           'PAP21','PNP11','PNP21',
                           'Seed','oh11','oh21')
totu <- t(otu[,-16])
targgenes <- t(arggenes)
##计算距离
#根据otu物种丰度数据，计算样方间的 Bray-curtis 距离
dist.abund <- vegdist(totu, method = 'bray')
#根据arg基因丰度表，计算样方间的欧几里得距离
dist.arg <- vegdist(targgenes, method = 'bray')
#物种丰度和arg的相关性，以 spearman 相关系数为例，9999 次置换检验显著性（Mantel tests 基于随机置换的方法获取 p 值）
abund_arg <- mantel(dist.abund, dist.arg, method = 'spearman', permutations = 9999, na.rm = TRUE)
abund_arg #otu与arg_r:0.2924,p:0.0434

#mrg与otu
mrggenes <- read.delim('D:/wenjian/毕业设计/大论文/数据/meta/MGE/normalized_cell.level3.txt',
                       row.names = 1)
colnames(mrggenes) <- c('CG11','CG21','Cd11',
                        'Cd21','Cr11','Cr21',
                        'Cu11','Cu21','PAP11',
                        'PAP21','PNP11','PNP21',
                        'Seed','oh11','oh21')
tmrggenes <- t(mrggenes)
#根据arg基因丰度表，计算样方间的bray距离
dist.mrg <- vegdist(tmrggenes, method = 'bray')
#物种丰度和arg的相关性，以 spearman 相关系数为例，9999 次置换检验显著性（Mantel tests 基于随机置换的方法获取 p 值）
abund_mrg <- mantel(dist.abund, dist.mrg, method = 'spearman', permutations = 9999, na.rm = TRUE)
abund_mrg #otu与mrg——r:0.3149,p:0.03
#otu与aromadeg——r:0.2425,p:0.0727;otu与MGE——r:0.4529,p:0.0058

#几个基因之间,mrg
mrggenes <- read.delim('D:/wenjian/毕业设计/大论文/数据/meta/MRG/normalized_cell.level3.txt',
                       row.names = 1)
colnames(mrggenes) <- c('CG11','CG21','Cd11',
                        'Cd21','Cr11','Cr21',
                        'Cu11','Cu21','PAP11',
                        'PAP21','PNP11','PNP21',
                        'Seed','oh11','oh21')
tmrggenes <- t(mrggenes)
#根据arg基因丰度表，计算样方间的bray距离
dist.mrg <- vegdist(tmrggenes, method = 'bray')
#Aromadeg
aromagenes <- read.delim('D:/wenjian/毕业设计/大论文/数据/meta/aromadeg/normalized_cell.level3.txt',
                       row.names = 1)
colnames(aromagenes) <- c('CG11','CG21','Cd11',
                        'Cd21','Cr11','Cr21',
                        'Cu11','Cu21','PAP11',
                        'PAP21','PNP11','PNP21',
                        'Seed','oh11','oh21')
taromagenes <- t(aromagenes)
#根据arg基因丰度表，计算样方间的bray距离
dist.aroma <- vegdist(taromagenes, method = 'bray')
#MGE
mgegenes <- read.delim('D:/wenjian/毕业设计/大论文/数据/meta/MGE/normalized_cell.level3.txt',
                       row.names = 1)
colnames(mgegenes) <- c('CG11','CG21','Cd11',
                        'Cd21','Cr11','Cr21',
                        'Cu11','Cu21','PAP11',
                        'PAP21','PNP11','PNP21',
                        'Seed','oh11','oh21')
tmgegenes <- t(mgegenes)
#根据arg基因丰度表，计算样方间的bray距离
dist.mge <- vegdist(tmgegenes, method = 'bray')

#arg与mrg的相关性，以 spearman 相关系数为例，9999 次置换检验显著性（Mantel tests 基于随机置换的方法获取 p 值）
arg_mrg <- mantel(dist.arg, dist.mrg, method = 'spearman', permutations = 9999, na.rm = TRUE)
arg_mrg #arg与mrg——r:0.6904,p:2e-04
#arg与aromadeg的相关性，以 spearman 相关系数为例，9999 次置换检验显著性（Mantel tests 基于随机置换的方法获取 p 值）
arg_aroma <- mantel(dist.arg, dist.aroma, method = 'spearman', permutations = 9999, na.rm = TRUE)
arg_aroma #arg与aroma——r:0.6344,p:1e-04
#arg与mge的相关性，以 spearman 相关系数为例，9999 次置换检验显著性（Mantel tests 基于随机置换的方法获取 p 值）
arg_mge <- mantel(dist.arg, dist.mge, method = 'spearman', permutations = 9999, na.rm = TRUE)
arg_mge #arg与mge——r:0.6562,p:2e-04
#mrg aroma的相关性，以 spearman 相关系数为例，9999 次置换检验显著性（Mantel tests 基于随机置换的方法获取 p 值）
mrg_aroma <- mantel(dist.mrg, dist.aroma, method = 'spearman', permutations = 9999, na.rm = TRUE)
mrg_aroma #mrg与aroma——r:0.6759,p:1e-04
#mrg mge的相关性，以 spearman 相关系数为例，9999 次置换检验显著性（Mantel tests 基于随机置换的方法获取 p 值）
mrg_mge <- mantel(dist.mrg, dist.mge, method = 'spearman', permutations = 9999, na.rm = TRUE)
mrg_mge #mrg与mge——r:0.72,p:1e-04
#aroma mge的相关性，以 spearman 相关系数为例，9999 次置换检验显著性（Mantel tests 基于随机置换的方法获取 p 值）
aroma_mge <- mantel(dist.aroma, dist.mge, method = 'spearman', permutations = 9999, na.rm = TRUE)
aroma_mge #aroma与mge——r:0.7432,p:2e-04

#绘制mantel结果热图
library(readxl)
setwd('D:/wenjian/毕业设计/大论文/数据/meta')
r <- read.csv('mantel.csv',row.names = 1)
p <- read.csv('mantelp.csv',row.names = 1)

library(reshape2)


# 将相关系数矩阵转换为数据框
r_df <- melt(as.matrix(r))

# 将显著性水平矩阵转换为数据框
p_df <- melt(as.matrix(p))

# 合并数据框
df <- cbind(r_df, significance = p_df$value)

# 添加显著性标记
df$signif_mark <- ifelse(df$significance < 0.05 & df$significance > 0.01, "*", 
                         ifelse(df$significance < 0.01 & df$significance > 0.001, "**", 
                                ifelse(df$significance < 0.001, "***", "")))

# 绘制热图
library(ggplot2)
ggplot(df, aes(Var2, Var1, fill = value)) +
  geom_tile(color = "white") +
  geom_text(aes(label = paste0(round(value, 2), signif_mark))) +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0, limit = c(0, 1)) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid = element_blank(),
        axis.title = element_blank()) +
  coord_fixed()



#linkET包中的mantel_test适合与环境因子/较少的变量,进而绘制与相关性热图的组合图
# mantel_test(spec=t(otu[,-16]),env=t(arggenes),
#        seed=12345,
#        method = 'pearson',
#        
#        na_omit=TRUE
#        )













#########################

#对ARGs等基因总丰度变化的比较，只与seed比较####
setwd("D:/wenjian/毕业设计/大论文/数据/meta/ARG")
type_arg <- read.delim('normalized_cell.type.txt',row.names = 1)
type_arg <- read.delim('ppm.type.txt',row.names = 1)
colnames(type_arg) <- c('CG11','CG21','Cd11',
                        'Cd21','Cr11','Cr21',
                        'Cu11','Cu21','PAP11',
                        'PAP21','PNP11','PNP21',
                        'Seed','oh11','oh21')
#读取metadata
md <- read.delim('D:/wenjian/毕业设计/大论文/数据/meta/metadata.txt')
library(dplyr)
md <- md %>% arrange(sample)
#把arg总丰度合并到metadata
abunarg <- type_arg %>% summarise_all(sum) %>% t() %>% as.data.frame()
colnames(abunarg) <- c('argabun')
md <- bind_cols(md,abunarg)
#删除CG11、CG21
md <- md[-c(1,2,13),]
#调整因子顺序，以便绘图时是想要的顺序
md$rongyu <- factor(md$rongyu,levels = c("y","n"))
md$category <- factor(md$category,levels = c("heavy_metal","organic","no"))
#计算显著性差异
library(rstatix) #add_xy_position
md_p_val1 <- md %>% group_by(category)%>%
  wilcox_test(argabun ~rongyu) %>%
  adjust_pvalue(p.col="p",method="bonferroni") %>%
  add_significance(p.col="p.adj") %>% 
  add_xy_position(x="rongyu",dodge=0.8)
#绘图
library(ggplot2)
library(ggpubr) #stat_pvalue_manual
library(ggh4x) #facet_nested_wrap
md %>% ggplot(aes(rongyu,argabun))+
   annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = 0.3390,
            fill = "cornflowerblue", alpha = .3, color = NA)+
  # annotate("rect", xmin = -Inf, xmax = Inf, ymin = 2, ymax = Inf,
  #          fill = "#FCAE12", alpha = .3, color = NA) +
  geom_violin(aes(fill=category),trim = FALSE,show.legend = F)+
  geom_boxplot(width = 0.2,outliers = FALSE, staplewidth = 0.5) +
  scale_y_continuous(sec.axis = sec_axis(~ ., name = ""))+
  stat_summary(fun=mean,geom="point",col="black",fill="#F98400",
               shape=23,show.legend = F)+
  stat_summary(fun=mean, geom="line", aes(group=category), col="black")+
  stat_pvalue_manual(md_p_val1,label="p.adj.signif",#hide.ns=T,
                     tip.length = 0,label.size = 5,color="black")+
  # geom_hline(yintercept = 2,linetype=2)+
   geom_hline(yintercept = 0.3390,linetype=2)+
  geom_vline(xintercept = 2.6,linetype=2)+
  facet_nested_wrap(. ~ category,
                    strip = strip_nested(background_x = 
                                           elem_list_rect(fill=c("cyan3","indianred1","#8f888b")))) +
  scale_fill_manual(values = c("cyan3","indianred1","#8f888b"))+
  labs(x=NULL,y='Copies of ARG per cell')+
  theme(axis.text.x=element_text(angle = 0,vjust=0.5,hjust=0.5,color="black"),
        axis.text.y=element_text(color="black"),
        plot.background = element_rect(fill="white"), 
        panel.background = element_rect(fill="white"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.spacing.x = unit(0,"cm"),
        plot.margin=unit(c(0.5,0.5,0.5,0.5),unit="cm"),
        axis.line.x=element_line(color="black"),
        axis.line.y.left = element_line(color="grey30"),
        axis.line.y.right = element_line(color="grey30"),
        axis.line.x.bottom = element_line(color="grey30"),
        axis.text.y.right = element_blank(),
        axis.ticks.y.right = element_blank())+
  guides(y = guide_axis(minor.ticks = TRUE))

# 单样本t检验,判断一列数字是否显著高于某个值。
t_test_result <- t.test(md[c(7,8,9,10,11,12),]$argabun, mu = 90, alternative = "greater")



#两种处理对arg-subtype层面的影响，绘制ARG基因的丰度热图####
setwd("D:/wenjian/毕业设计/大论文/数据/meta/ARG")
type_subarg <- read.delim('normalized_cell.subtype.txt',row.names = 1)
colnames(type_subarg) <- c('CG11','CG21','Cd11',
                        'Cd21','Cr11','Cr21',
                        'Cu11','Cu21','PAP11',
                        'PAP21','PNP11','PNP21',
                        'Seed','oh11','oh21')
row_max <- apply(type_subarg[,c(9:12,14,15)], MARGIN = 1, max) #使用apply函数获取每一行的最大值,MARGIN = 1指定了函数是沿着行应用的
type_subarg$MaxValue <- row_max #将最大值作为新的列添加到df中
subarg_sorted <- type_subarg[order(-type_subarg$MaxValue), ]# 根据MaxValue列的值对df进行降序排序
subarg_sorted <- head(subarg_sorted, 100)# 选取MaxValue值最大的前5行
guiyi <- log10(subarg_sorted[,c(-1,-2,-16)]*1000000 + 1) #参考Lee ISME做的归一化
library(pheatmap)
pheatmap(t(guiyi),
         cluster_cols = FALSE,
         color = colorRampPalette(c('white','blue','red'))(1000), #热图色块颜色是从蓝到红分为100个等级
         )

# 对数据框按照行名重新排序
guiyi <- guiyi[order(rownames(guiyi)), ]
guiyi <- guiyi[,order(colnames(guiyi)) ]

library(RColorBrewer)
Colors <- brewer.pal(9, "Blues")
pheatmap(t(guiyi),
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         #color = Colors
         color = colorRampPalette(c('azure','lightpink','red'))(1000), #热图色块颜色是从蓝到红分为100个等级
         #color = colorRampPalette(c('red','white','green'))(1000)
         )

library(ggplot2)
pheatmap(t(guiyi),
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         xlim=0.1,
         ylim=NA,
         #color = Colors
         color = colorRampPalette(c('#4079B0','#C65360'))(1000), #热图色块颜色是从蓝到红分为100个等级
         #color = colorRampPalette(c('red','white','green'))(1000)
         )

#pheatmap做多注释热图，包括ARG基因的type和机制####
setwd("D:/wenjian/毕业设计/大论文/数据/meta/ARG")
type_subarg <- read.delim('normalized_cell.subtype.txt',row.names = 1)
#setwd("D:/wenjian/毕业设计/大论文/数据/meta/argmeta/jili2")
#type_subarg <- read.delim('normalized_cell.Subtype.txt',row.names = 1)

colnames(type_subarg) <- c('CG11','CG21','Cd11',
                           'Cd21','Cr11','Cr21',
                           'Cu11','Cu21','PAP11',
                           'PAP21','PNP11','PNP21',
                           'Seed','oh11','oh21')
argname <- rownames(type_subarg)
type <- strsplit(argname,split = '__')
# 将拆分后的结果存储在数据框中
result_type <- data.frame(do.call(rbind, type))
df <- cbind(result_type,type_subarg)

# 使用duplicated函数检查第二列是否有重复值
duplicate_check <- duplicated(df[, 2])
# 打印出包含重复值的行
if (any(duplicate_check)) {
  print("第二列中存在重复值的行：")
  print(df[duplicate_check, ])
} else {
  print("第二列中没有重复值。")
}
#删除重复的行，其实不用
df <- df[df$X2!='amrB',]
rownames(df)<- df$X2

#取七大type
df <- df[df$X1 %in% c('bacitracin','multidrug','sulfonamide',
                       'aminoglycoside','macrolide-lincosamide-streptogramin',
                          'beta_lactam','polymyxin'),]

row_max <- apply(df[,c(-1,-2,-3,-4)], MARGIN = 1, max) #使用apply函数获取每一行的最大值,MARGIN = 1指定了函数是沿着行应用的
df$MaxValue <- row_max #将最大值作为新的列添加到df中
df_sorted <- df[order(-df$MaxValue), ]# 根据MaxValue列的值对df进行降序排序
df_sorted <- head(df_sorted, 100)# 选取MaxValue值最大的前5行
#取前100基因
guiyi <- log10(df_sorted[,c(-1,-2,-3,-4,-18)]*1000000 + 1) #参考Lee ISME做的归一化
guiyi <- guiyi[,order(colnames(guiyi))]
guiyi$type <- df_sorted$X1[1:100]
guiyi <- guiyi[order(guiyi$type),]

#构造行注释数据框
md <- read.delim('D:/wenjian/毕业设计/大论文/数据/meta/metadata.txt')#读取args数据
md <- md[order(md$sample),]
hang <- md[c(-3,-4),c(1,3)]
rownames(hang) <- hang$sample
library(dplyr)
hang <- select(hang,category)

###构造列注释数据框
guiyitype <- as.data.frame(guiyi$type)
rownames(guiyitype) <- rownames(guiyi)[1:100]
#读取arg机理文件
jili <- read.delim('D:/wenjian/毕业设计/大论文/数据库/Short_subdatabase_V3.2.1/Short_subdatabase/4.SARG_v3.2_20220917_Short_subdatabase_structure.txt')
# 使用separate函数将Subtype列根据"__"拆分成两列
library(tidyr)
jili_split <- separate(jili, Subtype, into = c("col1", "col2"), sep = "__")

ji2 <- jili_split[jili_split$col2 %in% rownames(guiyitype),]
ji2 <- select(ji2,col2,Mechanism.group)
# 使用distinct函数根据指定列去除冗余行
ji2_unique <- ji2 %>% distinct(col2, .keep_all = TRUE)
rownames(ji2_unique) <- ji2_unique[,1]

guiyitype$gene <- rownames(guiyitype)
lie <- merge(guiyitype,ji2_unique,by.x = 'gene',by.y = 'col2')

rownames(lie) <- lie$gene
lie <- lie[,-1]
colnames(lie)[1] <- 'type'
#颜色设置
mechcolor <- c("#85B22E","#5F80B4","#E29827","#922927",'#57C3F3','#F3B1A0') 
names(mechcolor) <- c("Antibiotic target alteration",
                                "Antibiotic target replacement",
                                "Efflux pump",
                                "Efflux pump RND family",
                                "Enzymatic inactivation",'') #类型颜色

#Agecolor <- colorRampPalette(c("white","#99CCCC","#66CC99","#339966"))(6)
#names(Agecolor) <- c("30","31","33","34","35","45")

#Sexcolor <- c("red","#016D06") 
#names(Sexcolor) <- c("F","M") #类型颜色

#argtype <- c("#708090",'#68A180','#F3B1A0', '#D6E7A3')
#names(BPcolor) <- c("Immune response","Proteoglycans in cancer","Glycolysis","Endocytosis")

categorycolor <- c("cyan3",'indianred1','#8f888b')
names(categorycolor) <- c("heavy_metal","organic","no")

typecolors <- c("#001F3F", "#807023", "#228B22", "#C0C0C0", "#4B0082", "#FF7F50", "#8A2BE2")
names(typecolors) <- c('bacitracin','multidrug','sulfonamide',
                       'aminoglycoside','macrolide-lincosamide-streptogramin',
                       'beta_lactam','polymyxin')
ann_colors <- list(Mechanism.group=mechcolor,category=categorycolor,
                   type=typecolors) #颜色设置

#绘图
library(pheatmap)
pheatmap(t(guiyi[,-14]),
         color = colorRampPalette(c('azure','lightpink','red'))(1000),
         cluster_cols = T,
         cluster_rows = FALSE,
         scale = 'none',
         #cutree_cols=3,
         annotation_col = lie,
         annotation_row = hang,  # 仅包括行注释列
         annotation_colors = ann_colors
         )


#对几种基因之间做回归拟合或者结构方程####
#读取metadata
md <- read.delim('D:/wenjian/毕业设计/大论文/数据/meta/metadata.txt')#读取args数据
library(dplyr)
md <- md %>% arrange(sample)
#把arg总丰度合并到metadata
setwd("D:/wenjian/毕业设计/大论文/数据/meta/ARG")
type_arg <- read.delim('normalized_cell.type.txt',row.names = 1)
colnames(type_arg) <- c('CG11','CG21','Cd11',
                        'Cd21','Cr11','Cr21',
                        'Cu11','Cu21','PAP11',
                        'PAP21','PNP11','PNP21',
                        'Seed','oh11','oh21')
abunarg <- type_arg %>% summarise_all(sum) %>% t() %>% as.data.frame()
colnames(abunarg) <- c('argabun')
md <- bind_cols(md,abunarg)
#把mrg总丰度合并到metadata
setwd("D:/wenjian/毕业设计/大论文/数据/meta/MRG")
type_mrg <- read.delim('normalized_cell.level2.txt',row.names = 1)
colnames(type_mrg) <- c('CG11','CG21','Cd11',
                        'Cd21','Cr11','Cr21',
                        'Cu11','Cu21','PAP11',
                        'PAP21','PNP11','PNP21',
                        'Seed','oh11','oh21')
abunmrg <- type_mrg %>% summarise_all(sum) %>% t() %>% as.data.frame()
colnames(abunmrg) <- c('mrgabun')
md <- bind_cols(md,abunmrg)
#把aromadeg总丰度合并到metadata
setwd("D:/wenjian/毕业设计/大论文/数据/meta/aromadeg")
type_aro <- read.delim('normalized_cell.level4.txt',row.names = 1)
colnames(type_aro) <- c('CG11','CG21','Cd11',
                        'Cd21','Cr11','Cr21',
                        'Cu11','Cu21','PAP11',
                        'PAP21','PNP11','PNP21',
                        'Seed','oh11','oh21')
abunaro <- type_aro %>% summarise_all(sum) %>% t() %>% as.data.frame()
colnames(abunaro) <- c('aroabun')
md <- bind_cols(md,abunaro)
#把mge总丰度合并到metadata
setwd("D:/wenjian/毕业设计/大论文/数据/meta/MGE")
type_mge <- read.delim('normalized_cell.level2.txt',row.names = 1)
colnames(type_mge) <- c('CG11','CG21','Cd11',
                        'Cd21','Cr11','Cr21',
                        'Cu11','Cu21','PAP11',
                        'PAP21','PNP11','PNP21',
                        'Seed','oh11','oh21')
abunmge <- type_mge %>% summarise_all(sum) %>% t() %>% as.data.frame()
colnames(abunmge) <- c('mgeabun')
md <- bind_cols(md,abunmge)
#回归分析
##t test
#t <- cor.test(md$mrgabun,md$aroabun)##p0.034  r:0.704
##draw plot and line picture
library(ggplot2)
library(ggpmisc)
ggplot(md[c(-1,-2,-13,-5),],aes(x=mrgabun,y=argabun,color=category))+
  geom_point(size=5)+
  geom_smooth(method="lm",se=T
              #,color="steelblue2"
              )+ #method 还有lm,glm,gam,loess,
  stat_poly_eq(use_label(c("eq", "adj.R2", "p.value.label")),
               formula = y ~ x,  parse = TRUE,
               size = 5 #公式字体大小
              # label.x = 0.05,  #位置 ，0-1之间的比例
              # label.y = 0.95
               )+
  scale_color_manual(values = c("heavy_metal" = "cyan3",
                                "organic" = "indianred1"
                                )) +  # 手动设置每个category的颜色
  ylim(0,1.7)+
  theme_classic()+
  #theme(legend.position = "none")+
  labs(x = "Copies of MRG per cell",y = 'Copies of ARG per cell')+
  theme(
    panel.background = element_blank(),  # 设置背景为透明
    panel.grid.major = element_blank(),  # 去掉主要网格线
    panel.grid.minor = element_blank(),  # 去掉次要网格线
    #panel.border = element_blank(),  # 去掉边界线
    panel.border = element_rect(color = "grey20", fill = NA, size = 0.5),  # 添加边界线，颜色为深灰色，填充为透明，大小为0.5
    axis.line = element_line(color = "grey20", size = 0.5),  # 设置坐标轴线条为深灰色，更加细致
    axis.title = element_text(face = "bold", size = 12, color = "grey20"),  # 设置坐标轴标题字体为粗体，大小为12，颜色为深灰色
    axis.text = element_text(size = 10, color = "grey20"),  # 设置坐标轴标签字体大小为10，颜色为深灰色
    axis.ticks = element_line(color = "grey20", size = 0.5),  # 设置坐标轴刻度线颜色为深灰色，更加细致
    #legend.position = "right",  # 设置图例位置在右侧
    legend.position = c(0.2,0.8),
    legend.title = element_text(face = "bold", size = 12, color = "grey20"),  # 设置图例标题字体为粗体，大小为12，颜色为深灰色
    legend.text = element_text(size = 10, color = "grey20"),  # 设置图例文本字体大小为10，颜色为深灰色
    legend.background = element_blank(),  # 设置图例背景为透明
    legend.box.background = element_rect(color = "grey80", size = 0.5),  # 设置图例框背景为浅灰色，边框更加细致
    legend.key = element_blank()  # 去掉图例键的背景
  )

#对ARG风险等级进行作图####
setwd("D:/wenjian/毕业设计/大论文/数据/meta/ARG")
dat <- read.delim('risk.txt')

md <- read.delim('D:/wenjian/毕业设计/大论文/数据/meta/metadata.txt')#读取args数据
library(dplyr)
md <- md %>% arrange(sample)

#
risk <- t(dat)
colnames(risk) <- risk[1,]
risk <- risk[-1,]
md <- cbind(md,risk)
#这里是30d+60d合起来的
library(reshape2)
risklevel <- melt(md[,c(1,3,8,9,10,11)],id.vars = c('category','sample'))
risklevel$value <- as.numeric(risklevel$value)
risklevel <- risklevel[risklevel$variable!='level4',]

#30d的
risklevel <- melt(md[md$days==60|md$days==0,c(1,3,8,9,10,11)],id.vars = c('category','sample'))
risklevel$value <- as.numeric(risklevel$value)
risklevel <- risklevel[risklevel$variable!='level4',]

##plot
library(ggplot2)
library(ggbreak)
library(ggpubr)
library(ggsignif)

ggplot(risklevel,aes(x=variable,y=value,fill=category))+
  geom_boxplot(aes(x=factor(variable,levels = c('level1','level2','level3','level4'))),
               outlier.shape = NA
               )+
  scale_fill_manual(values=c("cyan3","#8f888b","indianred1"))+
 # scale_y_break(c(5, 40),scales = 0.3,space = 0.01)+
  coord_cartesian(ylim = c(0, 2)) +
  labs(x="Risk level",y="the realative proportion of ARGs (%)")+ 
  stat_compare_means(aes(group=category),)+
  theme(
    panel.background = element_blank(),  # 设置背景为透明
    panel.grid.major = element_blank(),  # 去掉主要网格线
    panel.grid.minor = element_blank(),  # 去掉次要网格线
    #panel.border = element_blank(),  # 去掉边界线
    panel.border = element_rect(color = "grey20", fill = NA, size = 0.5),  # 添加边界线，颜色为深灰色，填充为透明，大小为0.5
    axis.line = element_line(color = "grey20", size = 0.5),  # 设置坐标轴线条为深灰色，更加细致
    axis.title = element_text(face = "bold", size = 12, color = "grey20"),  # 设置坐标轴标题字体为粗体，大小为12，颜色为深灰色
    axis.text = element_text(size = 10, color = "grey20"),  # 设置坐标轴标签字体大小为10，颜色为深灰色
    axis.ticks = element_line(color = "grey20", size = 0.5),  # 设置坐标轴刻度线颜色为深灰色，更加细致
    legend.position = "right",  # 设置图例位置在右侧
    legend.title = element_text(face = "bold", size = 12, color = "grey20"),  # 设置图例标题字体为粗体，大小为12，颜色为深灰色
    legend.text = element_text(size = 10, color = "grey20"),  # 设置图例文本字体大小为10，颜色为深灰色
    legend.background = element_blank(),  # 设置图例背景为透明
    legend.box.background = element_rect(color = "grey80", size = 0.5),  # 设置图例框背景为浅灰色，边框更加细致
    legend.key = element_blank()  # 去掉图例键的背景
    )


kruskal.test(data=risklevel[risklevel$variable=="level3",],value ~category)

#对ARG机理分析####
setwd("D:/wenjian/毕业设计/大论文/数据/meta/argmeta/jili2")
mech_group <- read.delim('normalized_cell.Mechanism.group.txt',row.names = 1)
colnames(mech_group) <- c('CG11','CG21','Cd11',
                        'Cd21','Cr11','Cr21',
                        'Cu11','Cu21','PAP11',
                        'PAP21','PNP11','PNP21',
                        'Seed','oh11','oh21')
library(reshape2)
mech_group2 <- t(mech_group)
#metadata
md <- read.delim('D:/wenjian/毕业设计/大论文/数据/meta/metadata.txt')#读取args数据
library(dplyr)
md <- md %>% arrange(sample)

# 将df2中的Value2列添加到df1中；注意，在向一个数值的数据框中添加字符型列，会把那些数字转变成字符类型
mechgroup_combined <- cbind(mech_group2, category = md$category,days=md$days) %>% as.data.frame()
# 删除CG11 CG21
mechgroup_combined <- mechgroup_combined[c(-1,-2,-13),]
#转换成长表，保留category days
long_mechgroup <- melt(mechgroup_combined,id.vars = c('category','days'))
#melted_df <- gather(mechgroup_combined, key = "variable", value = "value", -category,-days)
#把long_mechgroup中的字符型数字改为数值型
long_mechgroup$value<-as.numeric(long_mechgroup$value)

#绘图
library(ggplot2)
ggplot(long_mechgroup[long_mechgroup$days==60,],aes(x=variable,y=value,fill=category))+
  geom_boxplot(
    position=position_dodge(0.9)                 ##因为分组比较，需设组间距
  )+
  stat_boxplot(
    mapping=aes(x = reorder(variable, -value, FUN = median),y=value),
    geom ="errorbar",      ##添加箱子的bar为最大、小值
    width=0.3,
    position=position_dodge(0.9)   ###这里的0.9要与上面的0.9对应
  )+
  scale_y_continuous(limits=c(0,1))+      ##修改y轴的范围
  scale_fill_manual(values=c("cyan3","indianred1","#8f888b"))+ #设置分组的颜色
  labs(x="",y="Copies of ARG per cell (30 d)")+        ##修改坐标轴和图例的文本
  theme(
    panel.background = element_blank(),  # 设置背景为透明
    panel.grid.major = element_blank(),  # 去掉主要网格线
    panel.grid.minor = element_blank(),  # 去掉次要网格线
    #panel.border = element_blank(),  # 去掉边界线
    panel.border = element_rect(color = "grey20", fill = NA, size = 0.5),  # 添加边界线，颜色为深灰色，填充为透明，大小为0.5
    axis.line = element_line(color = "grey20", size = 0.5),  # 设置坐标轴线条为深灰色，更加细致
    axis.title = element_text(face = "bold", size = 12, color = "grey20"),  # 设置坐标轴标题字体为粗体，大小为12，颜色为深灰色
    axis.text = element_text(size = 10, color = "grey20"),  # 设置坐标轴标签字体大小为10，颜色为深灰色
    axis.text.x = element_text(angle = 15,vjust = 0.85,hjust = 0.75), #横轴旋转
    axis.ticks = element_line(color = "grey20", size = 0.5),  # 设置坐标轴刻度线颜色为深灰色，更加细致
    legend.position = "right",  # 设置图例位置在右侧
    legend.title = element_text(face = "bold", size = 12, color = "grey20"),  # 设置图例标题字体为粗体，大小为12，颜色为深灰色
    legend.text = element_text(size = 10, color = "grey20"),  # 设置图例文本字体大小为10，颜色为深灰色
    legend.background = element_blank(),  # 设置图例背景为透明
    legend.box.background = element_rect(color = "grey80", size = 0.5),  # 设置图例框背景为浅灰色，边框更加细致
    legend.key = element_blank()  # 去掉图例键的背景
    )

##手动显著性检验
try <- type_long
try$value <- try$value
kruskal.test(value~category,data=try)




#############################
#对三种基因arg、mrg、aroma的丰度表构建网络图####
aroma <- read.delim('D:/wenjian/毕业设计/大论文/数据/meta/aromadeg/normalized_cell.level1.txt',sep = '\t')
sarg <- read.csv('D:/wenjian/毕业设计/大论文/数据/meta/ARG/normalized_cell.subtype.txt',sep = '\t')
mrg <- read.csv('D:/wenjian/毕业设计/大论文/数据/meta/MRG/normalized_cell.level3.txt',sep = '\t')
## 合并数据框
colnames(aroma)[1] <- 'gene'
colnames(sarg)[1] <- 'gene'
colnames(mrg)[1] <- 'gene'
merged_df <- rbind(aroma,sarg,mrg)

## NA  > 0
library(tidyverse)
merged_df <- merged_df %>% mutate_all(~replace(., is.na(.), 0))
rownames(merged_df) <- merged_df[,1]
merged_df <- merged_df[,-1]
library(igraph)
library(Hmisc) #提供rcorr函数
library(ggplot2)
library(ggpubr)
library(ggsignif)
#将shu转成矩阵
cooccur_allshu <- as.matrix(t(merged_df))

#此处由于点不多，不进行低丰度/低出现率属的过滤

#计算相关性系数；
sp.cor<- rcorr(cooccur_allshu,type="spearman") #这里计算的是所有物种之间的
#提取r、p值矩阵；
r.cor<-sp.cor$r
p.cor<-sp.cor$P
#使用Benjamini-Hochberg("FDR-BH")法进行多重检验校正；
p.adj <- p.adjust(p.cor, method="BH")
#确定物种间存在相互作用关系的阈值，将相关性R矩阵内不符合的数据转换为0
r.cor[p.cor>=0.05|r.cor<0.6] = 0
#对角线处的1不计
diag(r.cor) <- 0

r.cor[1:1383,1:1383] <- 0
r.cor[1384:1858,1384:1858] <- 0
r.cor[1859:2062,1859:2062] <- 0

r.cor[1:1383,1858:2062] <- 0
r.cor[1858:2062,1:1383] <- 0

r.cor[1:1383,1384:1858]<-0 #只留下mrg与arg的关系
r.cor[1384:1858,1:1383]<-0 #只留下mrg与arg的关系

#r.cor[1384:1858,1859:2062]<-0 #只留下aroma与arg的关系
#r.cor[1859:2062,1384:1858]<-0 #只留下aroma与arg的关系

#使用邻接矩阵（即相关系数矩阵）创建网络；
net_allshu<-graph.adjacency(r.cor,weight=T,mode="undirected")
#去掉冗余的边（multiple edges、loop edges）；
net_allshu<-simplify(net_allshu)
#提取权重
df_weight = E(net_allshu)$weight
# 设定边的宽度，这里我们将相关系数与边宽进行关联
E(net_allshu)$width = abs(df_weight)*5
#生成网络图的结点标签（OTU id）和degree属性；degree表示一个点的中心性
V(net_allshu)$label <- V(net_allshu)$name
V(net_allshu)$degree <- degree(net_allshu)

V(net_allshu)$type <- c(rep('aroma',times=1383),rep('arg',times=475),
                        rep('mrg',times=204)
)

#查看网络图的对象结构;
print(net_allshu)

#将网络图导出为"graphml"、"gml"格式，方便导入Gephi中使用；
write_graph(net_allshu, "arg_mrg.graphml", format="graphml")

########################
#分析质粒####
#读取质粒数据
setwd("D:/wenjian/毕业设计/大论文/数据/meta")
plasmid <- read.delim('plasmid_result.txt')

#读取每个mags信息
#读取物种信息
setwd("D:/wenjian/毕业设计/大论文/数据/meta/MAGs")
library(readxl)
tax <- read.delim('gtdbtk.bac120.summary.tsv')
tax <- tax[,1:2]
max_cols <- max(sapply(strsplit(tax$classification, ";"), length))# 确定最大的拆分列数
col_names <- paste0("col", 1:max_cols)# 生成列名
library(tidyr)
taxonomy <- tax %>%
  separate(classification, into = col_names, sep = ";", fill = "right")# 拆分字符串到列

#把质粒信息合并
# 使用merge函数进行左连接，并将B中不存在的值设为0
merged_data <- merge(taxonomy, plasmid, by.x = "user_genome",by.y = 'File_Name',
                     all.x = TRUE)
merged_data[is.na(merged_data)] <- 0  # 将NA值替换为0

#变形菌门与其他门之间的比较
# 根据"genome"列分组
merged_data$men <- NA
merged_data[merged_data$col2 == "p__Pseudomonadota",]$men <- "Pseudomonadota"
merged_data[merged_data$col2 != "p__Pseudomonadota",]$men <- "other phylum"

# 绘制箱线图比较不同组之间的数值差异

library(ggplot2)
library(ggsignif)#geom_signif
library(ggpubr)#stat_pvalue_manual stat_compare_means
ggplot(merged_data, aes(x = men, y = Count, fill=men)) +
  geom_boxplot(outlier.shape = NA) +  # 隐藏箱线图中的异常值点
  geom_jitter(width = 0.2, alpha = 0.7)+  # 显示异常值点
  stat_compare_means(comparisons = list(c("Pseudomonadota", "other phylum")),
                     method = "wilcox.test", test = "p.adjust",label = 'p.signif')+
  labs(x=NULL,y='Number of plasmids (per MAG)')+
  scale_fill_manual(values = c("#8f888b","#81A88D"))+
  theme(
    panel.background = element_blank(),  # 设置背景为透明
    panel.grid.major = element_blank(),  # 去掉主要网格线
    panel.grid.minor = element_blank(),  # 去掉次要网格线
    #panel.border = element_blank(),  # 去掉边界线
    panel.border = element_rect(color = "grey20", fill = NA, size = 0.5),  # 添加边界线，颜色为深灰色，填充为透明，大小为0.5
    axis.line = element_line(color = "grey20", size = 0.5),  # 设置坐标轴线条为深灰色，更加细致
    axis.title = element_text(face = "bold", size = 12, color = "grey20"),  # 设置坐标轴标题字体为粗体，大小为12，颜色为深灰色
    axis.text = element_text(size = 10, color = "grey20"),  # 设置坐标轴标签字体大小为10，颜色为深灰色
    axis.text.x = element_text(angle = 0,vjust = 1,hjust = 0.5), #横轴旋转
    axis.ticks = element_line(color = "grey20", size = 0.5),  # 设置坐标轴刻度线颜色为深灰色，更加细致
    legend.position = c(0.8,0.7),  # 设置图例位置在右侧
    legend.title = element_text(face = "bold", size = 12, color = "grey20"),  # 设置图例标题字体为粗体，大小为12，颜色为深灰色
    legend.text = element_text(size = 10, color = "grey20"),  # 设置图例文本字体大小为10，颜色为深灰色
    legend.background = element_blank(),  # 设置图例背景为透明
    legend.box.background = element_rect(color = "grey80", size = 0.5),  # 设置图例框背景为浅灰色，边框更加细致
    legend.key = element_blank()  # 去掉图例键的背景
  )



#变形菌门中，伽马变形菌纲与其他纲之间的比较
# 根据"genome"列分组
merged_data2 <- merged_data[merged_data$col2 == "p__Pseudomonadota",]
merged_data2$gang <- NA
merged_data2[merged_data2$col3 == "c__Gammaproteobacteria",]$gang <- "Gammaproteobacteria"
merged_data2[merged_data2$col3 != "c__Gammaproteobacteria",]$gang <- "other class"
ggplot(merged_data2, aes(x = gang, y = Count, fill=gang)) +
  geom_boxplot(outlier.shape = NA) +  # 隐藏箱线图中的异常值点
  geom_jitter(width = 0.2, alpha = 0.7)+  # 显示异常值点
  stat_compare_means(comparisons = list(c("Gammaproteobacteria", "other class")),
                     method = "wilcox.test", test = "p.adjust",label = 'p.signif')+
  labs(x=NULL,y='Number of plasmids (per MAG)')+
  scale_fill_manual(values = c("#FCE6D5","#7294D4"))+
  theme(
    panel.background = element_blank(),  # 设置背景为透明
    panel.grid.major = element_blank(),  # 去掉主要网格线
    panel.grid.minor = element_blank(),  # 去掉次要网格线
    #panel.border = element_blank(),  # 去掉边界线
    panel.border = element_rect(color = "grey20", fill = NA, size = 0.5),  # 添加边界线，颜色为深灰色，填充为透明，大小为0.5
    axis.line = element_line(color = "grey20", size = 0.5),  # 设置坐标轴线条为深灰色，更加细致
    axis.title = element_text(face = "bold", size = 12, color = "grey20"),  # 设置坐标轴标题字体为粗体，大小为12，颜色为深灰色
    axis.text = element_text(size = 10, color = "grey20"),  # 设置坐标轴标签字体大小为10，颜色为深灰色
    axis.text.x = element_text(angle = 0,vjust = 1,hjust = 0.5), #横轴旋转
    axis.ticks = element_line(color = "grey20", size = 0.5),  # 设置坐标轴刻度线颜色为深灰色，更加细致
    legend.position = c(0.8,0.7),  # 设置图例位置在右侧
    legend.title = element_text(face = "bold", size = 12, color = "grey20"),  # 设置图例标题字体为粗体，大小为12，颜色为深灰色
    legend.text = element_text(size = 10, color = "grey20"),  # 设置图例文本字体大小为10，颜色为深灰色
    legend.background = element_blank(),  # 设置图例背景为透明
    legend.box.background = element_rect(color = "grey80", size = 0.5),  # 设置图例框背景为浅灰色，边框更加细致
    legend.key = element_blank()  # 去掉图例键的背景
  )

#计算MAGs中质粒非0比例
#伽马变形菌纲
per1 <- sum(merged_data[merged_data$col3 == "c__Gammaproteobacteria",]$Count != 0)/nrow(merged_data[merged_data$col3 == "c__Gammaproteobacteria",])*100
#变形菌门其他纲
per2 <- sum(merged_data2[merged_data2$col3 != "c__Gammaproteobacteria",]$Count != 0)/nrow(merged_data2[merged_data2$col3 != "c__Gammaproteobacteria",])*100
#其他门
per3 <- sum(merged_data[merged_data$col2 != "p__Pseudomonadota",]$Count != 0)/nrow(merged_data[merged_data$col2 != "p__Pseudomonadota",])*100

bili <- data.frame(name=c('Gammaproteobacteria','non-Gammaproteobacteria','other phylum'),
                   percentage =c(51.32,38.46,25.58)
                   )
ggplot(bili) +
  geom_bar(aes(x = reorder(name,-percentage), y = percentage),
           stat = "identity", position = "dodge",
           show.legend = T,alpha = .9,linewidth = 0.5,
           fill='skyblue'
  )+
  scale_y_continuous(limits = c(0,60),expand = c(0,0))+
  labs( x=NULL,y="Proportion of plasmids present in MAGs (%)") +
  theme(
    panel.background = element_blank(),  # 设置背景为透明
    panel.grid.major = element_blank(),  # 去掉主要网格线
    panel.grid.minor = element_blank(),  # 去掉次要网格线
    #panel.border = element_blank(),  # 去掉边界线
    panel.border = element_rect(color = "grey20", fill = NA, size = 0.5),  # 添加边界线，颜色为深灰色，填充为透明，大小为0.5
    axis.line = element_line(color = "grey20", size = 0.5),  # 设置坐标轴线条为深灰色，更加细致
    axis.title = element_text(face = "bold", size = 12, color = "grey20"),  # 设置坐标轴标题字体为粗体，大小为12，颜色为深灰色
    axis.text = element_text(size = 10, color = "grey20"),  # 设置坐标轴标签字体大小为10，颜色为深灰色
    axis.text.x = element_text(angle = 0,vjust = 1,hjust = 0.5), #横轴旋转
    axis.ticks = element_line(color = "grey20", size = 0.5),  # 设置坐标轴刻度线颜色为深灰色，更加细致
    legend.position = c(0.8,0.7),  # 设置图例位置在右侧
    legend.title = element_text(face = "bold", size = 12, color = "grey20"),  # 设置图例标题字体为粗体，大小为12，颜色为深灰色
    legend.text = element_text(size = 10, color = "grey20"),  # 设置图例文本字体大小为10，颜色为深灰色
    legend.background = element_blank(),  # 设置图例背景为透明
    legend.box.background = element_rect(color = "grey80", size = 0.5),  # 设置图例框背景为浅灰色，边框更加细致
    legend.key = element_blank()  # 去掉图例键的背景
  )+
  guides(y = guide_axis(minor.ticks = TRUE))

#####################补充分析，reviewers#########################

#R1 对各个处理组绘制ARG总丰度############
setwd("D:/wenjian/毕业设计/大论文/数据/meta/ARG")

type_arg <- read.delim('normalized_cell.type.txt',row.names = 1)
type_arg <- read.delim('rpkm.type.txt',row.names = 1)
colnames(type_arg) <- c('CK11','CK21','Cd11',
                        'Cd21','Cr11','Cr21',
                        'Cu11','Cu21','PAP11',
                        'PAP21','PNP11','PNP21',
                        'Seed','phOH11','phOH21')
#读取metadata
md <- read.delim('D:/wenjian/毕业设计/大论文/数据/meta/metadata.txt')
library(dplyr)
md <- md %>% arrange(sample)
#把arg总丰度合并到metadata
abunarg <- type_arg %>% summarise_all(sum) %>% t() %>% as.data.frame()
colnames(abunarg) <- c('argabun')
md <- bind_cols(md,abunarg)
#调整因子顺序，以便绘图时是想要的顺序
md$rongyu <- factor(md$rongyu,levels = c("y","n"))
md$category <- factor(md$category,levels = c("heavy_metal","organic","no"))

md <- read.delim('D:/wenjian/毕业设计/大论文/数据/meta/1reviewers/metadata22.txt')
md$sample2 <- c('Cd11',
                'Cd21','Cr11','Cr21',
                'Cu11','Cu21','PAP11',
                'PAP21','PNP11','PNP21',
                'Seed','CK11','CK21','phOH11','phOH21')

library(ggplot2)
ARG_total_abundance <- 
  ggplot(md, aes(x = sample2, y = argabun,fill = category)) +
  geom_col(aes(y = argabun)) +
  geom_vline(xintercept = c(7.5,14.5),linetype=2)+
  scale_fill_manual(values = c("heavy_metal"="cyan3","organic"="indianred1",
                               "p:Acidobacteriota"="#dadada","p:Actinobacteriota"="#fbf398",
                               "bbfasi"="orchid","p:Planctomycetota"="#e77381",
                               "p:Bacteroidota"="#9b8191",'no'= "#8f888b"),)+
  scale_x_discrete(limits=c('Cd11','Cr11','Cu11','phOH11','PNP11','PAP11','CK11',
                            'Cd21','Cr21','Cu21','phOH21','PNP21','PAP21','CK21',
                            'Seed'))+
  labs(x=NULL,y='Copies of ARG per cell')+
  theme_minimal()+  # 使用简洁主题
  theme(axis.text.x=element_text(angle = 30,vjust=0.5,hjust=0.5,color="black"),
        axis.text.y=element_text(color="black"),
        plot.background = element_rect(fill="white"), 
        panel.background = element_rect(fill="white"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.spacing.x = unit(0,"cm"),
        legend.position = 'none',
        plot.margin=unit(c(0.5,0.5,0.5,0.5),unit="cm"),
        axis.line.x=element_line(color="black"),
        axis.line.y.left = element_line(color="grey30"),
        axis.line.y.right = element_line(color="grey30"),
        axis.line.x.bottom = element_line(color="grey30"),
        axis.text.y.right = element_blank(),
        axis.ticks.y.right = element_blank())+
  guides(y = guide_axis(minor.ticks = TRUE))
ARG_total_abundance
setwd("D:/wenjian/毕业设计/大论文/数据/meta/1reviewers")
ggsave("ARG_total_abundance.pdf",ARG_total_abundance,width = 4.5,height = 4,units = 'in')

#R2 分各个type画各个处理组的ARG丰度，用热图画
setwd("D:/wenjian/毕业设计/大论文/数据/meta/ARG")

type_arg <- read.delim('normalized_cell.type.txt',row.names = 1)
#type_arg <- read.delim('rpkm.type.txt',row.names = 1)
colnames(type_arg) <- c('CK11','CK21','Cd11',
                        'Cd21','Cr11','Cr21',
                        'Cu11','Cu21','PAP11',
                        'PAP21','PNP11','PNP21',
                        'Seed','phOH11','phOH21')

library(dplyr)
setwd("D:/wenjian/毕业设计/大论文/数据/meta/1reviewers")
arg_ck2 <- read.delim("2normalized_cell.type.txt",row.names = 1) %>% as.data.frame()
arg_ck2 <- select(arg_ck2,B814_final_pure_reads)

type_arg$argtype <- rownames(type_arg)
arg_ck2$argtype <- rownames(arg_ck2)
arg_new <- merge(type_arg,arg_ck2,by = 'argtype')

arg_ck1 <- read.delim("normalized_cell.type.txt",row.names = 1) %>% as.data.frame()
arg_ck1 <- select(arg_ck1,seed2_final_pure_reads)
arg_ck1$argtype <- rownames(arg_ck1)
arg_new <- merge(arg_new,arg_ck1,by = 'argtype')

arg_new <- arg_new[,c(-2,-3)]
colnames(arg_new)[15:16] <- c("CK21",'CK11')
rownames(arg_new) <- arg_new[,1]
arg_new <- arg_new[,-1]

# 计算每列的总和
column_sums <- apply(arg_new, 2, sum)
# 转换为百分比
df_percentage <- sweep(arg_new, 2, column_sums, "/")

library(reshape2)
df_percentage$type <- rownames(df_percentage)
df_long <- melt(df_percentage)

library(ggplot2)
colorsss <- c(
     "#17becf",  # 青色
    "#ff7f0e",  # 橙色
 # "#F5E320",  # 红色
 # "#CC281C",  # 绿色
  "#D1CB7E",  # 紫色
  "#5AAAB2",  # 棕色
  "#E77F5A",  # 粉色
  "#A0E6C5",  # 灰色
  "#70A47F",  # 橄榄色
  "#1f77b4",  # 蓝色
  "#D68112",  # 浅蓝色
  "#C9A79E",  # 浅橙色
  "#98df8a",  # 浅绿色
  #  "#ff9896",  # 浅红色
  "#6C6B6F",  # 浅紫色
  "#CED2DD",  # 浅棕色
  #  "#f7b6d2",  # 浅粉色
  "#4E616A",  # 浅橄榄色
  #  "#9edae5",  # 浅青色
  "#ED6764",  # 深蓝色
  "#74A0A4"   # 深绿色
)
# 绘制百分比堆叠柱状图 
p2_argtype_percentage<-  
  ggplot(df_long, aes(x = variable, y = value, fill = type)) +
  geom_bar(stat = "identity", width = 1) +  # 使用实际值绘制柱状图
  scale_fill_manual(values = colorsss)+
  scale_x_discrete(limits=c('Cd11','Cr11','Cu11','phOH11','PNP11','PAP11','CK11',
                            'Cd21','Cr21','Cu21','phOH21','PNP21','PAP21','CK21',
                            'Seed'))+
  coord_flip() +  # 翻转坐标轴
  scale_y_continuous(limits=c(0,1),expand = c(0,0),labels = scales::percent) +  # 将y轴标签转换为百分比
  labs(x = "", y = "Percentage", fill = "ARG Types") +  # 设置轴标签和图例标题
  theme_minimal()+  # 使用简洁主题
  geom_vline(xintercept = c(7.5,14.5),linetype=2)+
  guides(fill=guide_legend(title = "",nrow =4))+
  theme(axis.text.x=element_text(angle = 0,vjust=0.5,hjust=0.5,color="black"),
        axis.text.y=element_text(color="black"),
        plot.background = element_rect(fill="white"), 
        panel.background = element_rect(fill="white"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.spacing.x = unit(0,"cm"),
        legend.position = 'top',
        plot.margin=unit(c(0.5,0.5,0.5,0.5),unit="cm"),
        axis.line.x=element_line(color="black"),
        axis.line.y.left = element_line(color="grey30"),
        axis.line.y.right = element_line(color="grey30"),
        axis.line.x.bottom = element_line(color="grey30"),
        axis.text.y.right = element_blank(),
        axis.ticks.y.right = element_blank())+
  guides(y = guide_axis(minor.ticks = TRUE))

p2_argtype_percentage
setwd("D:/wenjian/毕业设计/大论文/数据/meta/1reviewers")

ggsave("p2_argtype_percentage.pdf",p2_argtype_percentage,width = 7.2,height = 6.8,units = 'in')


