
#read in data，相关分析####
aroma <- read.csv('aromadeg_gene_out_table.csv',sep = '\t')
sarg <- read.csv('sarg_gene_out_table.csv',sep = '\t')
mrg <- read.csv('mrg_gene_out_table2.csv',sep = '\t')

## 合并数据框
merged_df <- merge(merge(aroma, sarg, by = "X", all = TRUE), mrg, by = "X", all = TRUE)

## NA  > 0
library(tidyverse)
merged_df <- merged_df %>% mutate_all(~replace(., is.na(.), 0))
rownames(merged_df) <- merged_df[,1]
merged_df <- merged_df[,-1]


library(corrplot)
library(ggplot2)
env.cor <- round(cor(merged_df, method = "spearman"),3) # round(),对输出结果取小数点前三位,计算pearson相关性系数
env.cor
library(ggcorrplot)
env.p <-round(cor_pmat(merged_df,method = "spearman"),3) # pearson系数
env.p




par(mfrow=c(1,1))
p1<- corrplot(corr =env.cor, p.mat = env.p,method = "square",
         type="upper",diag=FALSE,add=TRUE,
         tl.pos="n", tl.cex=0.5, tl.col="black",tl.srt = 45,tl.offset=0.5,
         insig="label_sig",sig.level = c(.001, .01, .05),
         pch.cex = 0.8,pch.col = "red") # tl.*参数依次设置变量标签的位置为左侧和右侧，字符大小为1，颜色为黑色，倾斜45度，距离热图边界0.5。pch*参数依次设置显著性标签大小为0.8，颜色为红色。
p2<-corrplot(corr = env.cor, method = "number",
         type="lower",
         tl.pos = "tp",tl.cex=0.5, tl.col="black",tl.srt = 45,tl.offset=0.5,
         cl.pos = "n",
         number.digits = 2,number.cex = 0.5,number.font = NULL
) # number*参数依次设置保留相关性系数3位小数点，字体大小为0.7，字体使用默认，只能设置par中的字体参数。后面会讲修改par字体参数。

#画这些ARG在各个样品中的分布####
setwd("D:/wenjian/毕业设计/大论文/数据/meta/ARG")
subtype_arg <- read.delim('normalized_cell.subtype.txt',row.names = 1)
colnames(subtype_arg) <- c('CG11','CG21','Cd11',
                        'Cd21','Cr11','Cr21',
                        'Cu11','Cu21','PAP11',
                        'PAP21','PNP11','PNP21',
                        'Seed','oh11','oh21')

argme <- subtype_arg[rownames(subtype_arg) %in% colnames(sarg)[-1],]


colnames(sarg)[-1]




#构建网络图####
aroma <- read.csv('aromadeg_gene_out2_table.csv',sep = '\t')
sarg <- read.csv('sarg_gene_out1_table.csv',sep = '\t')
mrg <- read.csv('mrg_gene_out_table.csv',sep = '\t')
## 合并数据框
merged_df <- merge(merge(aroma, sarg, by = "X", all = TRUE), mrg, by = "X", all = TRUE)

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
cooccur_allshu <- as.matrix(merged_df)

#此处由于点不多，不进行低丰度/低出现率属的过滤

#计算相关性系数；
sp.cor<- rcorr(cooccur_allshu,type="spearman") #这里计算的是所有物种之间的
#提取r、p值矩阵；
r.cor<-sp.cor$r
p.cor<-sp.cor$P
#使用Benjamini-Hochberg("FDR-BH")法进行多重检验校正；
p.adj <- p.adjust(p.cor, method="BH")
#确定物种间存在相互作用关系的阈值，将相关性R矩阵内不符合的数据转换为0
r.cor[p.cor>=0.05|abs(r.cor)<0.5] = 0
#对角线处的1不计
diag(r.cor) <- 0

r.cor[1:119,1:119] <- 0
r.cor[120:186,120:186] <- 0
r.cor[187:233,187:233] <- 0

r.cor[1:119,187:223] <- 0
r.cor[187:223,1:119] <- 0

#r.cor[1:119,120:186]<-0 #只留下mrg与arg的关系
#r.cor[120:186,1:119]<-0 #只留下mrg与arg的关系

r.cor[120:186,187:233]<-0 #只留下aroma与arg的关系
r.cor[187:233,120:186]<-0 #只留下aroma与arg的关系

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

V(net_allshu)$type <- c(rep('aroma',times=119),rep('arg',times=67),
                        rep('mrg',times=47)
                        )

#查看网络图的对象结构;
print(net_allshu)

#将网络图导出为"graphml"、"gml"格式，方便导入Gephi中使用；
write_graph(net_allshu, "argyuaroma.graphml", format="graphml")


#绘制mrg特殊ARGs的累积丰度####
#读取特殊ARGs
mrgarg <- read.csv('arg_mrg.csv')
mrgarg <- mrgarg[mrgarg$eigencentrality>0.5,]
mrgarg <- mrgarg[mrgarg$v_type == 'arg',]
#读取ARGs丰度表
subtype_arg <- read.delim('D:/wenjian/毕业设计/大论文/数据/meta/ARG/normalized_cell.subtype.txt',row.names = 1)
colnames(subtype_arg) <- c('CG11','CG21','Cd11',
                           'Cd21','Cr11','Cr21',
                           'Cu11','Cu21','PAP11',
                           'PAP21','PNP11','PNP21',
                           'Seed','oh11','oh21')

argme <- subtype_arg[rownames(subtype_arg) %in% c(mrgarg$v_name,
                                                  'multidrug__Pseudomonas aeruginosa CpxR',
                                                  'multidrug__YajC'
                                                  ),]
##转换成长表
library(reshape2)
argme$gene <- rownames(argme)
visdata <- melt(argme)

#######根据不同category绘制###
zone <- visdata
zone <- zone[]
zone$zo <- c(rep('no',times=12),rep('heavy-metal',times=36),rep('organic',times=24)
             ,rep('no',times=6),rep('organic',times=12))

library(ggplot2)
ggplot(data=zone,aes(variable,value,fill=zo))+
  ###geom_bar是绘制柱状图的函数
  geom_bar(stat="identity",
           position="stack", #如果是针对分组的柱形图，则position除了可以"identity"(不调整，组内前后重叠)、还包括“stack“（堆积，默认）；"fill"(按比例堆积)；“dodge“（分散开）
           # color="black", 
           size=0.25,
           width=0.5)+
  ###scale_fill_manual可以自己设置颜色
  #scale_fill_manual(values=c("red", "blue", "green","darkgreen","black"))+
  scale_y_continuous(expand = c(0,0))+ # 坐标轴延伸，确保图形元素覆盖至坐标
  labs(x = "",y = "realative abundance (genus)")+   #添加x，y轴名
  theme_classic()+ # 主题类型
  # scale_fill_discrete(limits=c("p:Proteobacteria","p:Verrucomicrobiota","p:Acidobacteriota",
  #                              "p:Actinobacteriota","p:Gemmatimonadota","p:Planctomycetota",
  #                              "p:Bacteroidota",)) #调整图例顺序
  scale_fill_manual(values = c("heavy-metal"="#98d09d","organic"="#d7e698",
                               "no"="#dadada","p:Actinobacteriota"="#fbf398",
                               "bbfasi"="#f7a895","p:Planctomycetota"="#e77381",
                               "p:Bacteroidota"="#9b8191","#8f888b"),
  )+
  ### 设置主题
  theme(panel.background=element_rect(fill='transparent'),   ##去掉底层阴影
        panel.grid=element_blank(),                          ##去掉网格线
        #  panel.border=element_blank(),                        ##去掉图的边界线
        panel.border=element_rect(fill='transparent',        ##设置图的边界线
                                  color='transparent'),
        axis.line=element_line(colour="black",size=0.9),     ##设置坐标轴的线条
        axis.title=element_text(face="bold",size = 13),      ##设置坐标轴文本字体
        # axis.title.x=element_blank(),                        ##删除x轴文本               
        axis.text=element_text(face = "bold",size = 10),     ##设置坐标轴标签字体
        axis.ticks=element_line(color='black'),              ##设置坐标轴刻度线
        # legend.position = 'none',
        # legend.position=c(0.15, 0.9),                        ##设置图例的位置
        # legend.direction = "horizontal",                      ##设置图列标签的排列方式
        #legend.title=element_text(face = "bold",size=12),    ##设置图例标题字体
        legend.title=element_blank(),                        ##去掉图例标题         
        # legend.text=element_text(face = "bold",size=12),     ##设置图例文本的字体
        legend.background=element_rect(linetype="solid",     ##设置图列的背景和线条颜色
                                       colour ="black"),
        axis.text.x = element_text(angle = 60,              ##设置横坐标文字的旋转角度
                                   vjust = 0.85,
                                   hjust = 0.75
        )
  )  


#绘制aroma特殊ARGs的累积丰度####
#读取特殊ARGs
aromaarg <- read.csv('arg_aroma.csv')
aromaarg <- aromaarg[aromaarg$eigencentrality>0.5,]
aromaarg <- aromaarg[aromaarg$v_type == 'arg',]
#读取ARGs丰度表
subtype_arg <- read.delim('D:/wenjian/毕业设计/大论文/数据/meta/ARG/normalized_cell.subtype.txt',row.names = 1)
colnames(subtype_arg) <- c('CG11','CG21','Cd11',
                           'Cd21','Cr11','Cr21',
                           'Cu11','Cu21','PAP11',
                           'PAP21','PNP11','PNP21',
                           'Seed','oh11','oh21')

argme <- subtype_arg[rownames(subtype_arg) %in% c(aromaarg$v_name,
                                                  'macrolide-lincosamide-streptogramin__macA',
                                                  'multidrug__Acinetobacter baumannii AmvA',
                                                  'quinolone__Acinetobacter baumannii AbaQ',
                                                  'beta_lactam__OXA-332',
                                                  'fosfomycin__Acinetobacter baumannii AbaF',
                                                  'beta_lactam__ADC-44'
                                                  ),]
##转换成长表
library(reshape2)
argme$gene <- rownames(argme)
visdata <- melt(argme)

#######根据不同category绘制###
zone <- visdata
zone <- zone[]
zone$zo <- c(rep('no',times=30),rep('heavy-metal',times=90),rep('organic',times=60)
             ,rep('no',times=15),rep('organic',times=30))

ggplot(data=zone,aes(variable,value,fill=zo))+
  ###geom_bar是绘制柱状图的函数
  geom_bar(stat="identity",
           position="stack", #如果是针对分组的柱形图，则position除了可以"identity"(不调整，组内前后重叠)、还包括“stack“（堆积，默认）；"fill"(按比例堆积)；“dodge“（分散开）
           # color="black", 
           size=0.25,
           width=0.5)+
  ###scale_fill_manual可以自己设置颜色
  #scale_fill_manual(values=c("red", "blue", "green","darkgreen","black"))+
  scale_y_continuous(expand = c(0,0))+ # 坐标轴延伸，确保图形元素覆盖至坐标
  labs(x = "",y = "realative abundance (genus)")+   #添加x，y轴名
  theme_classic()+ # 主题类型
  # scale_fill_discrete(limits=c("p:Proteobacteria","p:Verrucomicrobiota","p:Acidobacteriota",
  #                              "p:Actinobacteriota","p:Gemmatimonadota","p:Planctomycetota",
  #                              "p:Bacteroidota",)) #调整图例顺序
  scale_fill_manual(values = c("heavy-metal"="#98d09d","organic"="#d7e698",
                               "no"="#dadada","p:Actinobacteriota"="#fbf398",
                               "bbfasi"="#f7a895","p:Planctomycetota"="#e77381",
                               "p:Bacteroidota"="#9b8191","#8f888b"),
  )+
  ### 设置主题
  theme(panel.background=element_rect(fill='transparent'),   ##去掉底层阴影
        panel.grid=element_blank(),                          ##去掉网格线
        #  panel.border=element_blank(),                        ##去掉图的边界线
        panel.border=element_rect(fill='transparent',        ##设置图的边界线
                                  color='transparent'),
        axis.line=element_line(colour="black",size=0.9),     ##设置坐标轴的线条
        axis.title=element_text(face="bold",size = 13),      ##设置坐标轴文本字体
        # axis.title.x=element_blank(),                        ##删除x轴文本               
        axis.text=element_text(face = "bold",size = 10),     ##设置坐标轴标签字体
        axis.ticks=element_line(color='black'),              ##设置坐标轴刻度线
        # legend.position = 'none',
        # legend.position=c(0.15, 0.9),                        ##设置图例的位置
        # legend.direction = "horizontal",                      ##设置图列标签的排列方式
        #legend.title=element_text(face = "bold",size=12),    ##设置图例标题字体
        legend.title=element_blank(),                        ##去掉图例标题         
        # legend.text=element_text(face = "bold",size=12),     ##设置图例文本的字体
        legend.background=element_rect(linetype="solid",     ##设置图列的背景和线条颜色
                                       colour ="black"),
        axis.text.x = element_text(angle = 60,              ##设置横坐标文字的旋转角度
                                   vjust = 0.85,
                                   hjust = 0.75
        )
  )  





#mrg特殊ARGs的机理饼图####
#读取特殊ARGs
mrgarg <- read.csv('arg_mrg.csv')
mrgarg <- mrgarg[mrgarg$eigencentrality>0.5,]
mrgarg <- mrgarg[mrgarg$v_type == 'arg',]
#读取机理文件
jili <- read.delim('D:/wenjian/毕业设计/大论文/数据库/SARG/ARG_rank.txt')

mrgmejili <- jili[jili$Subtype %in% c(mrgarg$v_name,
                                'multidrug__Pseudomonas aeruginosa CpxR',
                                'multidrug__YajC'
                                ),] 
library(dplyr)
mrgmejili <- distinct(mrgmejili[,-1])
#全为外排泵
string_table <- table(mrgmejili$Mechanism.subgroup)

# 绘制饼图
pie(string_table, labels = paste(names(string_table), 
                                 "(", string_table, ")", sep = ""),
    main = "Pie Chart of String Vector")



#aroma特殊ARGs的机理饼图####
#读取特殊ARGs
aromaarg <- read.csv('arg_aroma.csv')
aromaarg <- aromaarg[aromaarg$eigencentrality>0.5,]
aromaarg <- aromaarg[aromaarg$v_type == 'arg',]

argmejili <- jili[jili$Subtype %in% c(aromaarg$v_name,
                                                  'macrolide-lincosamide-streptogramin__macA',
                                                  'multidrug__Acinetobacter baumannii AmvA',
                                                  'quinolone__Acinetobacter baumannii AbaQ',
                                                  'beta_lactam__OXA-332',
                                                  'fosfomycin__Acinetobacter baumannii AbaF',
                                                  'beta_lactam__ADC-44'
                                               ),]
library(dplyr)
argmejili <- distinct(argmejili[,-1],Subtype, .keep_all = TRUE)

# 统计每个字符串出现的次数
string_table <- table(argmejili$Mechanism.subgroup)

# 绘制饼图
pie(string_table, labels = paste(names(string_table), 
                                 "(", string_table, ")", sep = ""),
    main = "Pie Chart of String Vector")




