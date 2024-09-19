#16s

#read in data and treatment####
setwd('D:/wenjian/毕业设计/大论文/数据/16s')
data_folder='D:/wenjian/毕业设计/大论文/数据/16s'

#读取丰度信息
oturaw <- read.delim2('oturaw.txt',header = T,sep = "\t",row.names = 1)
library(vegan)
otu <- as.data.frame(t(rrarefy(t(oturaw),min(colSums(oturaw))))) #otu表抽平

#读取物种信息
tax <- read.delim("sintax.txt",header = F)
library(tidyr)
library(dplyr)
max_cols <- max(sapply(strsplit(tax$V4, ","), length))# 确定最大的拆分列数
col_names <- paste0("col", 1:max_cols)# 生成列名
taxonomy <- tax %>%
  separate(V4, into = col_names, sep = ",", fill = "right")# 拆分字符串到列

#合并物种与丰度
otu$asvid <- rownames(otu)
asv <- merge(otu,taxonomy,by.x = 'asvid',by.y = 'V1')[,c(-17,-18)]
##得到属级别的丰度表
shu <- select(asv,-col1,-col2,-col3,-col4,-col5,-col7)
rownames(shu) <- shu$asvid
shu <- select(shu,-asvid)
shu$col6[is.na(shu$col6)] <- paste("unassigned") #将NA改为unassigned
shu <- shu %>% 
  group_by(col6) %>% 
  summarise_all(sum) %>% 
  as.data.frame() #将重复属合并
rownames(shu) <- shu$col6
shu <- select(shu,-col6)

#计算各种α多样性指数
otu <- select(otu, -asvid)
richness <- colSums(otu > 0)
shannon <- diversity(t(otu),"shannon")
simpson <- t(otu) %>% diversity("simpson")
fisher <- t(otu) %>% diversity("invsimpson")
diversityindex <- data.frame(SampleID = rownames(t(otu)),
                             ShannonIndex = shannon,
                             simpson = simpson,
                             fisher = fisher,
                             richness = richness)

#读取metadata
metadata <- read.delim("metadata.txt",header = 1)

#替换字符串中的某个字符，可以用gsub(pattern = ',',replacement = '|',rowname)


#计算某个门的相对丰度占比并绘图####
##得到men级别的丰度表
men <- select(asv,-col1,-col6,-col3,-col4,-col5,-col7)
rownames(men) <- men$asvid
men <- select(men,-asvid)
men$col2[is.na(men$col2)] <- paste("unassigned") #将NA改为unassigned
men <- men %>% 
  group_by(col2) %>% 
  summarise_all(sum) %>% 
  as.data.frame() #将重复属合并
rownames(men) <- men$col2
men <- select(men,-col2)
##绘制物种堆积图,准备数据

# 根据value列的值对数据框的行重新排序
men <- men[order(men$Cd11,decreasing = T), ]
# 合并第三行及之后的所有行的值相加成一行
others <- colSums(men[9:nrow(men), ])
men <- rbind(men[1:8,],others)
rownames(men)[9] <- 'Others'

#计算百分比矩阵
# 计算每列数据相对该列之和的百分比并存入新的数据框
men_percent <- as.data.frame(lapply(men, function(x) x/sum(x)))

#转换
library(reshape2)
men$men <- rownames(men)
men_percent$men <- rownames(men)
df1 <- melt(men,)
df2 <- melt(men_percent,)
df <- cbind(df1,df2)
colnames(df)[4:6]<- c('men2','sample','percent')
#绘图
library(ggplot2)
library(wesanderson)
df %>%
  ggplot(aes(x = percent, y = variable,fill = men))+
  geom_bar(stat = "identity", position = "stack")+
  geom_text(aes(label = if_else(percent > 0.05, paste0(scales::percent(percent,accuracy = 0.1)), NA)), 
            position = position_fill(vjust = 0.5),size = 3, color = "#000000")+
  labs(x=NULL,y=NULL)+
  scale_x_continuous(expand = c(0,0))+
 # scale_fill_manual(values=wes_palette(9,"GrandBudapest2",type="continuous"))+
  theme(#legend.position = "none",
        plot.margin = unit(c(0.2,0.5,0.2,0.2), "cm"),
        plot.background = element_rect(color = NA, fill = "#F2F2F2"))

#计算某个特定门中某个纲的相对丰度占比####
##先把变形菌门挑出来
rownames(asv) <- asv$asvid
asv_prote <- asv[asv$col2=='p:Proteobacteria',]
asv_prote <- select(asv_prote,-col1,-col2,-col4,-col5,-col6,-col7)
asv_prote <- select(asv_prote,-asvid)
asv_prote$col3[is.na(asv_prote$col3)] <- paste("unassigned") #将NA改为unassigned
asv_prote <- asv_prote %>% 
  group_by(col3) %>% 
  summarise_all(sum) %>% 
  as.data.frame() #将重复属合并
rownames(asv_prote) <- asv_prote$col3
asv_prote <- select(asv_prote,-col3)

#将第三行的NA改成0
asv_prote[is.na(asv_prote)] <- 0

##绘制物种堆积图,准备数据

# 根据value列的值对数据框的行重新排序
asv_prote <- asv_prote[order(asv_prote$Cd11,decreasing = T), ]
# 合并第三行及之后的所有行的值相加成一行
#others <- colSums(men[3:nrow(men), ])
#men <- rbind(men[1:2,],others)
#rownames(men)[3] <- 'others'

#计算百分比矩阵
# 计算每列数据相对该列之和的百分比并存入新的数据框
asvprote_percent <- as.data.frame(lapply(asv_prote, function(x) x/sum(x)))

#转换
library(reshape2)
asv_prote$asvpro <- rownames(asv_prote)
asvprote_percent$asvpro <- rownames(asv_prote)
df1 <- melt(asv_prote,)
df2 <- melt(asvprote_percent,)
df <- cbind(df1,df2)
colnames(df)[4:6]<- c('gang2','sample','percent')
#绘图
library(wesanderson)
df %>%
  ggplot(aes(x = percent, y = variable,fill = asvpro))+
  geom_bar(stat = "identity", position = "stack")+
  geom_text(aes(label = if_else(percent > 0.05, paste0(scales::percent(percent,accuracy = 0.1)), NA)), 
            position = position_fill(vjust = 0.5),size = 3, color = "#000000")+
  labs(x=NULL,y=NULL)+
  scale_x_continuous(expand = c(0,0))+
  scale_fill_manual(values=wes_palette("Moonrise3",type="discrete"))+ #GrandBudapest2颜色板也还行
  theme(#legend.position = "none",
        plot.margin = unit(c(0.2,0.5,0.2,0.2), "cm"),
        plot.background = element_rect(color = NA, fill = "#F2F2F2"))


#绘制某个纲在各个样品中的丰度####
##挑出伽马变形菌纲
asv_gama <- asv[asv$col3=='c:Gammaproteobacteria' & asv$col3!='NA',]
asv_gama <- select(asv_gama,-col1,-col2,-col4,-col5,-col6,-col7)
asv_gama <- select(asv_gama,-asvid,-col3)
# 去掉数据框中所有包含NA值的行
asv_gama <- asv_gama[complete.cases(asv_gama), ]
# 计算数据框中每列的和作为一个新的行
sum_row <- as.data.frame(colSums(asv_gama, na.rm = TRUE))
sum_row$sample2 <- rownames(sum_row)
colnames(sum_row) <- c('gamapro','sample2')
#伽马变形菌的丰度合并到metadata
metadata <- bind_cols(metadata,sum_row)
#调整因子顺序，以便绘图时是想要的顺序
metadata$rongyu <- factor(metadata$rongyu,levels = c("y","n"))
metadata$category <- factor(metadata$category,levels = c("heavy_metal","organic","no"))
#计算显著性差异
library(rstatix) #add_xy_position
metadata_p_val1 <- metadata %>% group_by(category)%>%
  wilcox_test(gamapro ~rongyu) %>%
  adjust_pvalue(p.col="p",method="bonferroni") %>%
  add_significance(p.col="p.adj") %>% 
  add_xy_position(x="rongyu",dodge=0.8)

#绘图
library(ggplot2)
library(ggpubr) #stat_pvalue_manual
library(ggh4x) #facet_nested_wrap
metadata %>% ggplot(aes(rongyu,gamapro))+
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
  stat_pvalue_manual(metadata_p_val1,label="p.adj.signif",hide.ns=T,
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



#绘制α多样性的柱状图+箱线图####
index <- merge(diversityindex,metadata,by.x = "SampleID",by.y = 'sample')
index30d <- index[c(1,3,5,7,9,11,13,15),]
index60d <- index[c(2,4,6,8,10,12,14,15),]
library(ggplot2)
#绘制柱状图（30d）####
ggplot(index30d,aes(x=SampleID,y=simpson,fill=category))+
  geom_bar(stat="identity",
           position="dodge", #如果是针对分组的柱形图，则position除了可以"identity"(不调整，组内前后重叠)、还包括“stack“（堆积，默认）；"fill"(按比例堆积)；“dodge“（分散开）
           # color="black", 
           size=0.25,
           width=0.5)+
  ###scale_fill_manual可以自己设置颜色
  labs(x = "",y = "simpson")+   #添加x，y轴名
  theme_classic()+ # 主题类型
  # scale_fill_discrete(limits=c("p:Proteobacteria","p:Verrucomicrobiota","p:Acidobacteriota",
  #                              "p:Actinobacteriota","p:Gemmatimonadota","p:Planctomycetota",
  #                              "p:Bacteroidota",)) #调整图例顺序
  scale_fill_manual(values = c("heavy_metal"="cyan3","organic"="indianred1",
                               "p:Acidobacteriota"="#dadada","p:Actinobacteriota"="#fbf398",
                               "bbfasi"="orchid","p:Planctomycetota"="#e77381",
                               "p:Bacteroidota"="#9b8191","#8f888b"),)+
  scale_x_discrete(limits=c('Cd11','Cr11','Cu11','oh11','PNP11','PAP11','CG11','Seed'))+
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

#绘制柱状图（60d）####
ggplot(index60d,aes(x=SampleID,y=simpson,fill=category))+
  geom_bar(stat="identity",
           position="dodge", #如果是针对分组的柱形图，则position除了可以"identity"(不调整，组内前后重叠)、还包括“stack“（堆积，默认）；"fill"(按比例堆积)；“dodge“（分散开）
           # color="black", 
           size=0.25,
           width=0.5)+
  ###scale_fill_manual可以自己设置颜色
  labs(x = "",y = "simpson")+   #添加x，y轴名
  theme_classic()+ # 主题类型
  # scale_fill_discrete(limits=c("p:Proteobacteria","p:Verrucomicrobiota","p:Acidobacteriota",
  #                              "p:Actinobacteriota","p:Gemmatimonadota","p:Planctomycetota",
  #                              "p:Bacteroidota",)) #调整图例顺序
  scale_fill_manual(values = c("heavy_metal"="cyan3","organic"="indianred1",
                               "p:Acidobacteriota"="#dadada","p:Actinobacteriota"="#fbf398",
                               "bbfasi"="orchid","p:Planctomycetota"="#e77381",
                               "p:Bacteroidota"="#9b8191","#8f888b"),)+
  scale_x_discrete(limits=c('Cd21','Cr21','Cu21','oh21','PNP21','PAP21','CG21','Seed'))+
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

  


#绘制丰度前50的属的热图####
library(pheatmap)
row_max <- apply(shu, MARGIN = 1, max) #使用apply函数获取每一行的最大值,MARGIN = 1指定了函数是沿着行应用的
shu$MaxValue <- row_max #将最大值作为新的列添加到df中
shu_sorted <- shu[order(-shu$MaxValue), ]# 根据MaxValue列的值对df进行降序排序
top50_shu <- head(shu_sorted, 50)# 选取MaxValue值最大的前5行
top50_shu <- log10(top50_shu[,-16] + 1)

#删除CG11 CG21,对列排序
head(top50_shu)
top50_shu <- top50_shu[,c(-3,-4)]
top50_shu <- top50_shu[,order(colnames(top50_shu))]
top50_shu <- select(top50_shu,Cd11,Cd21,Cr11,Cr21,Cu11,Cu21,oh11,oh21,PAP11,PAP21,PNP11,PNP21,CG11,CG21,Seed)
#准备行注释文件
md <- read.delim('D:/wenjian/毕业设计/大论文/数据/meta/metadata.txt')#读取args数据
md <- md[order(md$sample),]
#hang <- md[c(-3,-4),c(1,3)]
hang <- md[,c(1,3)]
rownames(hang) <- hang$sample
library(dplyr)
hang <- select(hang,category)

#注释颜色
categorycolor <- c("cyan3",'indianred1','#8f888b')
names(categorycolor) <- c("heavy_metal","organic","no")
ann_colors <- list(category=categorycolor) #颜色设置
#绘图
pheatmap(top50_shu,
        border_color = NA,
        scale = 'none',
        # border_color='grey60',
        #cutree_cols=4,
        #cutree_rows=3,#把行聚为3大类
        #annotation_col=annot_col,
        #annotation_row = anno_row,
        show_rownames=T,
        cluster_cols = FALSE,#取消对列的聚类
        annotation_col = hang,  # 仅包括行注释列
        annotation_colors = ann_colors,
        angle_col = 45
        )

#otu的nmds分析####
#calculating stress values for NMDS,判断是否能用nmds分析。stress<=0.2,stressplot图中点基本都在线附近
dat <- t(otu)
# write.table(dat,file.path(data_folder,"intermediate.csv"),row.names=TRUE,na="", sep=",") # to write an intermediate file
# dat<- read.delim(file.path(data_folder,"intermediate.csv"),sep=",",row.names=1) # to read in an intermediate file (subset data)
# file.remove(file.path(data_folder,"intermediate.csv"))
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
col<-as.vector(metadata$categorycol) #color are respond with location
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
adonis <-  adonis2(dat ~ rongyu, metadata, permutations = 999, method = "bray")
adonis$R2
###基于bray-curtis距离进行anosim分析
anosim = anosim(dat, metadata$rongyu, permutations = 999, distance = "bray")
R <- anosim$statistic
p <- anosim$signif
plot(anosim)







#lefse####
##修改行名物种注释的格式
library(stringr)
rowname <- tax$V4
changename <- gsub(pattern = ',',replacement = '|',rowname)
changename <- gsub(pattern = ':',replacement = '_',changename)
changename <- sub(pattern = 'd',replacement = 'k',changename)
tax$V5 <- changename

asv <- merge(tax,otu,by.x = 'V1',by.y = 'asvid')[,-c(1:4)]
#去重复，否则会在imagegp中报错。
library(dplyr)
asv <- asv %>% 
  group_by(V5) %>% 
  summarise_all(sum) %>% 
  as.data.frame()

##将meta分组数据放到otu表中，方便lefse分析(imagegp只能有一行分组，linux里可以有多行分组)
asvma <- as.matrix(asv)
rownames(asvma) <- asvma[,1]
asvma <- asvma[,-1]

category <- metadata$category
# Pollution <- meta[,17]
# Sample <- meta[,1]
asvma <- rbind(category,asvma)
asvma30d <- asvma[,c(1,5,7,9,11,13)]
asvma60d <- asvma[,c(2,6,8,10,12,14)]
write.table(asvma,"asv_lefse_category.txt",sep = "\t",quote = F,row.names = T,col.names = F)
write.table(asvma30d,"asv_lefse_30d.txt",sep = "\t",quote = F,row.names = T,col.names = F)
write.table(asvma60d,"asv_lefse_60d.txt",sep = "\t",quote = F,row.names = T,col.names = F)

# asvnew <- read.table('asv_lefse.txt',sep = "\t")
# changename <- asvnew[,1]
# changename <- gsub(pattern = '_',replacement = '__',changename)
# asvnew[,1] <- changename
# asvnew <- asvnew[c(-1,-3),-53:-57]
# 
# pollutionname <- asvnew[1,]
# pollutionname <- gsub(pattern = 'no',replacement = 'no_PPCPs',pollutionname)
# pollutionname <- gsub(pattern = 'antibiotics',replacement = 'with_PPCPs',pollutionname)
# 
# asvnew[1,] <- pollutionname
# 
# write.table(asvnew,'asv_imagegp_lefse_NEW_PPCPs.txt',sep = "\t",quote = F,row.names = F,col.names = F)



#co-occurnet####
#用所有样本构建，划分module，然后看每个module能不能代表一种处理
#载入所需R包；
library(igraph)
library(Hmisc) #提供rcorr函数
library(ggplot2)
library(ggpubr)
library(ggsignif)
#将shu转成矩阵
cooccur_allshu <- as.matrix(shu)

#此处由于点不多，不进行低丰度/低出现率属的过滤

#计算相关性系数；
sp.cor<- rcorr(t(cooccur_allshu),type="spearman") #这里计算的是所有物种之间的
#提取r、p值矩阵；
r.cor<-sp.cor$r
p.cor<-sp.cor$P
#使用Benjamini-Hochberg("FDR-BH")法进行多重检验校正；
p.adj <- p.adjust(p.cor, method="BH")
#确定物种间存在相互作用关系的阈值，将相关性R矩阵内不符合的数据转换为0
r.cor[p.cor>=0.05|r.cor<0.6] = 0
#对角线处的1不计
diag(r.cor) <- 0

#使用邻接矩阵（即相关系数矩阵）创建网络；
net_allshu<-graph.adjacency(r.cor,weight=T,mode="undirected")
#去掉冗余的边（multiple edges、loop edges）；
net_allshu<-simplify(net_allshu)

#加入otu的门级别信息
matched_indices <- match(rownames(shu), taxonomy$col6)
matched_values <- taxonomy$col2[matched_indices]
V(net_allshu)$phy <- as.character(matched_values)
#加入不同处理组的物种是否存在的信息
metal_30d <- shu[,c(1,5,7)]
metal_30d$sum <- rowSums(metal_30d)
metal_30d[metal_30d>1] <- 1
V(net_allshu)$metal30d <- metal_30d$sum #第三十天，三种重金属

organic_30d <- shu[,c(9,11,13)]
organic_30d$sum <- rowSums(organic_30d)
organic_30d[organic_30d>1] <- 1
V(net_allshu)$organic30d <- organic_30d$sum #第三十天，三种有机物

metal_60d <- shu[,c(2,6,8)]
metal_60d$sum <- rowSums(metal_60d)
metal_60d[metal_60d>1] <- 1
V(net_allshu)$metal60d <- metal_60d$sum #第60天，三种重金属

organic_60d <- shu[,c(10,12,14)]
organic_60d$sum <- rowSums(organic_60d)
organic_60d[organic_60d>1] <- 1
V(net_allshu)$organic60d <- organic_60d$sum #第三十天，三种有机物

#提取权重
df_weight = E(net_allshu)$weight
# 设定边的宽度，这里我们将相关系数与边宽进行关联
E(net_allshu)$width = abs(df_weight)*5
#生成网络图的结点标签（OTU id）和degree属性；degree表示一个点的中心性
V(net_allshu)$label <- V(net_allshu)$name
V(net_allshu)$degree <- degree(net_allshu)

#查看网络图的对象结构;
print(net_allshu)

#将网络图导出为"graphml"、"gml"格式，方便导入Gephi中使用；
write_graph(net_allshu, "all_otushu_positive.graphml", format="graphml")
write_graph(net_allshu, "all_otushu_positive.gml", format="gml")


#随机森林，30d的数据、60d的数据、微生物冗余性分开做####
#30d的数据
# 1.1 导入数据
## 后续合并数据表时，物种名必须一样，所以导入的数据表尽量不要有特殊符号。
#dir()
#file.show("otu.csv")
random_otu <- as.data.frame(t(otu[,index30d$SampleID]))
random_otu$category <- index30d$category
dim(random_otu)
head(random_otu)

random_otu <- random_otu[-8,]

# 1.2 计算相对丰度
spe = random_otu
spe <- sweep(spe[1:ncol(spe)-1],1,rowSums(spe[1:ncol(spe)-1]),'/')*100
spe$category <- random_otu$category
head(spe)
# 2. random forest随机森林分析
library(randomForest)
set.seed(12345) #保证结果可重复
RF.best = randomForest(x = spe[-c(1132)], # 此处需要替换为自己的数据
                       y = factor(spe[,1132]),# 此处需要替换为自己的数据
                       importance=TRUE,# 输出变量重要性排序时必须设置
                      # maxnodes = 5,# 此处设置数值是之前调参的结果。
                      # mtry = 80,# 此处设置数值是之前调参的结果
                      # ntree = 3,# 此处设置数值是之前调参的结果
                       Perm = 999)
#RF.best
rfPermute::confusionMatrix(RF.best) 
summary(RF.best) # 会同时绘制分类准确率折线图。

#3.变量重要性图
# 3.1 提取每个变量对样本分类的重要性
##RF1$importance #包含分类数+2列数据，
##每个自变量对每个分类的平均正确性降低值(mean descrease in accuracy),
##后两列分别为变量对所有分类的MeanDecreaseAccuracy和MeanDecreaseGini(节点不纯度减少值)。
##两个值越大，变量的重要性越大。
#RF.best$importanceSD # 变量重要值的置换检验的标准误，最后一列为MeanDecreaseAccuracy置换检验的p值。
imp = data.frame(importance(RF.best),MDA.p = RF.best$importanceSD[4])
head(imp) # 提取变量重要性值及置换检验p值.

# 3.2 将变量按MeanDecreaseGini重要性降序排列
library(dplyr)
imp = arrange(imp,desc(MeanDecreaseGini)) 
head(imp) 

## 输出重要性排序结果到本地
write.csv(imp,"importance.csv",quote = FALSE)

# 3.3 提取重要性top10变量绘制条形图

## top10
imp = imp[1:10,]
imp

## 随机森林变量重要性条形图
##变量按MeanDecreaseGini降序排列：reorder(rownames(imp),MeanDecreaseGini)
library(ggplot2)
p1= ggplot(imp,aes(x=MeanDecreaseGini,y=reorder(rownames(imp),MeanDecreaseGini)))+
  geom_bar(position = position_dodge(),
           width = 0.5,
           stat = "identity",
           fill="steelblue")+ # 柱子的宽度与位置要保持一致，拼图时左/右才能时柱子对齐。
  theme_minimal() +
  xlab("Mean Decrease in Gini Index")+
  scale_y_discrete(expand=c(0,0))+
  scale_x_continuous(expand=c(0,0))+
  theme(axis.text.y = element_text(size = 16,colour = "black"),
        axis.text.x = element_text(size=14,color="black"),
        axis.title.x.bottom = element_text(size=16,color="black"),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        #panel.grid = element_blank() # 去除网格线
  )
p1

#4.单因素非参数差异检验图，对随机森林分析变量top10进行kruskal-Wallis检验并提取显著性标记。
#这里只标记单因素kruskal-Wallis差异检验结果，后面再写如何标记多重比较的结果。
#可以参考[R统计绘图-多变量单因素非参数差异检验及添加显著性标记图](https://mp.weixin.qq.com/s?__biz=MzIzMTgwMjk1NQ==&amp;mid=2247485484&amp;idx=1&amp;sn=f20464b6d5c4a8838808585a35e16d7c&amp;chksm=e89fd908dfe8501e3af6248b1ab73732827e4f45e247e59d863f1e7bbf881713fb19be691486&token=1100054514&lang=zh_CN#rd)，自行修改代码。

# 4.1 单因素非参数差异检验
library(rstatix)
library(dplyr)
## 提取top10变量丰度数据
BC = spe[c('category',rownames(imp))]
BC

## 批量对10个变量进行kruskal-Wallis检验及提取显著性标记
sig = matrix(data=NA,nrow = 10,ncol = 2,dimnames = list(names(BC[2:11]),c("species","sig")))
for (i in 2:11){
  assign(paste(names(BC[i]),"kw",sep="."),kruskal.test(BC[,i],as.factor(BC[,1])));
  tmp = get(paste(names(BC[i]),"kw",sep="."))
  sig[i-1,] = c(names(BC[i]),ifelse(tmp$p.value >0.05,"ns","*"));
} # assign()将kruskal.test()的运行结果赋予指定变量名
kw = mget(paste(names(BC[2:11]),"kw",sep=".")) # mget()提取变量值到kw对象
kw$OTU510.kw # 列表对象
##将结果输出到本地
capture.output(kw,file = "kw.list.txt",append = FALSE)

## 提取显著性标记
sig = as.data.frame(sig)
head(sig)
# 4.2 计算数据均值与标准误
BC1 = data.frame(name = rownames(BC),BC[c(1:11)])
BC1 

## 数据形式由"宽"变“长”
BC1 = gather(BC1,key="species","abundance",-name,-category)
BC1

## 计算均值和标准误
BC2 = BC1 %>% group_by(category,species) %>% 
  get_summary_stats(abundance,type = "mean_se") # 标准差用"mean_sd"计算
BC2

# 4.3 合并均值数据、显著性标记和imp值
da = merge(BC2,sig,by.x = c("species"),by.y = c("species"))
da$group = factor(da$category,levels = c("A","B","C"))
da[da$category %in% c("A","B"),"sig"] = NA # 组件差异检验显著性标记，保留再一个处理中。
da2 = merge(da,imp,by.x = c("species"),by.y = "row.names")
da2
# 4.4 绘制物种丰度条形图
p2 = ggplot(da2,
            aes(x=mean,y=reorder(species,MeanDecreaseGini),fill=category))+
  geom_bar(position = position_dodge(),
           width = 0.5,
           stat = "identity")+
  geom_errorbar(aes(xmin=mean-se,xmax=mean+se),
                width=0.3,
                color="black",size=0.5,
                position = position_dodge(0.5))+
  geom_text(aes(x=(mean+se),label=sig,group=category),
            hjust =0,nudge_x=0.5, # 显著性标签水平调整
            nudge_y=0.075, # 显著性标签垂直调整，hjust/vjust:0(right/bottom)和1(top/left) 。
            angle = -90, # 逆时针旋转90°
            size=4.75,fontface="bold")+
  theme_minimal() +
  xlab("Abundance")+
  scale_y_discrete(expand=c(0,0))+
  scale_x_continuous(expand=c(0,0))+
  theme(axis.text.y = element_blank(), # 物种顺序一致，为了后续拼图，隐藏物种名。
        #axis.text.y = element_text(size=16,color="black"),
        axis.text.x = element_text(size=14,color="black"),
        axis.title.x.bottom = element_text(size=16,color="black"),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        #panel.grid = element_blank() # 去除网格线
        legend.text = element_text(size=14,color="black"),
        legend.title = element_text(size=16,color="black")
  )
p2 # 柱子高低不一致，显著性标记输出到本地后，需要手动调整一下。

# 拼图
library(patchwork)
pdf("随机森林丰度组合图.pdf", width = 12,height = 8,family = "Times")
p1+p2+plot_layout(nrow = 1, widths = c(2, 2),
                  guides = "collect")
dev.off()

##交叉验证帮助选择特定数量的 OTUs;暂时停滞***
#5 次重复十折交叉验证
set.seed(123)
otu_train.cv <- replicate(5, rfcv(spe[-ncol(spe)], spe$category, cv.fold = 5,step = 1.5), simplify = FALSE)
otu_train.cv
#提取验证结果绘图
otu_train.cv <- data.frame(sapply(otu_train.cv, '[[', 'error.cv'))
otu_train.cv$otus <- rownames(otu_train.cv)
otu_train.cv <- reshape2::melt(otu_train.cv, id = 'otus')
otu_train.cv$otus <- as.numeric(as.character(otu_train.cv$otus))
#拟合线图
library(ggplot2)
library(splines)  #用于在 geom_smooth() 中添加拟合线，或者使用 geom_line() 替代 geom_smooth() 绘制普通折线

p <- ggplot(otu_train.cv, aes(otus, value)) +
  geom_smooth(se = FALSE,  method = 'glm', formula = y~ns(x, 6)) +
  theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent')) +
  labs(title = '',x = 'Number of OTUs', y = 'Cross-validation error')

p

#NCM,这里不太适用？####
#查看样品序列条数
summary(colSums(otu))
#将OTU表转置为行名为样品名称，列名为OTU名称的表
otu_ncm<-t(otu)
head(rowSums(otu_ncm))
#取重金属30d样品
ncm_metal30d<-otu_ncm[c(1,5,7),] #重金属-30d
# ncm_metal30d<-otu_ncm[c(9,11,13),] #有机物处理-30d
# ncm_metal30d<-otu_ncm[c(2,6,8),] #重金属-60d
# ncm_metal30d<-otu_ncm[c(10,12,14),] #有机物处理-60d
# ncm_metal30d<-otu_ncm[c(1,5,7,2,6,8),] #重金属-all-6个
# ncm_metal30d<-otu_ncm[c(9,11,13,10,12,14),] #有机物-all-6个
# ncm_metal30d<-otu_ncm[c(3,4,15),] #不加污染物-3个

#Get immigration probility m and R2
#using Non-linear least squares (NLS) to calculate R2:
#spp: A community table with taxa as rows and samples as columns
spp<-ncm_metal30d
#remove other variable to save memory
#rm(otu2);rm(otu_niche);rm(otu_biom2);rm(otu);
#get the mean of abundance of each sample
N <- mean(apply(spp, 1, sum))
#get the mean of species relative abundance in metacommmunity
p.m <- apply(spp, 2, mean)
p.m <- p.m[p.m != 0]
#get the percentage of each species in metacommmunity
p <- p.m/N
#get the binary data of community abundance matrix
spp.bi <- 1*(spp>0)
#get the frequncy of species occurrence in metacommunity
freq <- apply(spp.bi, 2, mean)
freq <- freq[freq != 0]
#get a table record species percentage and occurrence frequency in metacommunity
C <- merge(p, freq, by=0)
#sort the table according to occurence frquency of each species
C <- C[order(C[,2]),]
C <- as.data.frame(C)
#delete rows containning zero
C.0 <- C[!(apply(C, 1, function(y) any(y == 0))),]
p <- C.0[,2]
freq <- C.0[,3]
names(p) <- C.0[,1]
names(freq) <- C.0[,1]
d = 1/N
##Fit model parameter m (or Nm) using Non-linear least squares (NLS)
library(minpack.lm)
m.fit <- nlsLM(freq ~ pbeta(d, N*m*p, N*m*(1-p), lower.tail=FALSE),start=list(m=0.1))
m.fit #get the m value

m.ci <- confint(m.fit, 'm', level=0.95)

freq.pred <- pbeta(d, N*coef(m.fit)*p, N*coef(m.fit)*(1-p), lower.tail=FALSE)
library(Hmisc)
pred.ci <- binconf(freq.pred*nrow(spp), nrow(spp), alpha=0.05, method="wilson", return.df=TRUE)
# get the R2 value
Rsqr <- 1 - (sum((freq - freq.pred)^2))/(sum((freq - mean(freq))^2))
Rsqr

#Data visulization
bacnlsALL <-data.frame(p,freq,freq.pred,pred.ci[,2:3])
inter.col<-rep('black',nrow(bacnlsALL))
inter.col[bacnlsALL$freq <= bacnlsALL$Lower]<-'#A52A2A'#define the color of below points
inter.col[bacnlsALL$freq >= bacnlsALL$Upper]<-'#29A6A6'#define the color of up points
library(grid)
grid.newpage()
pushViewport(viewport(h=0.6,w=0.6))
pushViewport(dataViewport(xData=range(log10(bacnlsALL$p)), yData=c(0,1.02),extension=c(0.02,0)))
grid.rect()
grid.points(log10(bacnlsALL$p), bacnlsALL$freq,pch=20,gp=gpar(col=inter.col,cex=0.7))
grid.yaxis()
grid.xaxis()
grid.lines(log10(bacnlsALL$p),bacnlsALL$freq.pred,gp=gpar(col='blue',lwd=2),default='native')

grid.lines(log10(bacnlsALL$p),bacnlsALL$Lower ,gp=gpar(col='blue',lwd=2,lty=2),default='native') 
grid.lines(log10(bacnlsALL$p),bacnlsALL$Upper,gp=gpar(col='blue',lwd=2,lty=2),default='native')  
grid.text(y=unit(0,'npc')-unit(2.5,'lines'),label='Mean Relative Abundance (log10)', gp=gpar(fontface=2)) 
grid.text(x=unit(0,'npc')-unit(3,'lines'),label='Frequency of Occurance',gp=gpar(fontface=2),rot=90) 
#grid.text(x=unit(0,'npc')-unit(-1,'lines'), y=unit(0,'npc')-unit(-15,'lines'),label='Mean Relative Abundance (log)', gp=gpar(fontface=2)) 
#grid.text(round(coef(m.fit)*N),x=unit(0,'npc')-unit(-5,'lines'), y=unit(0,'npc')-unit(-15,'lines'),gp=gpar(fontface=2)) 
#grid.text(label = "Nm=",x=unit(0,'npc')-unit(-3,'lines'), y=unit(0,'npc')-unit(-15,'lines'),gp=gpar(fontface=2))
#grid.text(round(Rsqr,2),x=unit(0,'npc')-unit(-5,'lines'), y=unit(0,'npc')-unit(-16,'lines'),gp=gpar(fontface=2))
#grid.text(label = "Rsqr=",x=unit(0,'npc')-unit(-3,'lines'), y=unit(0,'npc')-unit(-16,'lines'),gp=gpar(fontface=2))
draw.text <- function(just, i, j) {
  grid.text(paste("Rsqr=",round(Rsqr,3),"\n","Nm=",round(coef(m.fit)*N)), x=x[j], y=y[i], just=just)
  #grid.text(deparse(substitute(just)), x=x[j], y=y[i] + unit(2, "lines"),
  #          gp=gpar(col="grey", fontsize=8))
}
x <- unit(1:4/5, "npc")
y <- unit(1:4/5, "npc")
draw.text(c("centre", "bottom"), 4, 1)

#C-score####
#链接到 github 安装 EcoSimR 包
# library(devtools)
# install_github('GotelliLab/EcoSimR')

#加载 EcoSimR 包
library(EcoSimR)

#读取示例数据
comm <- otu

#根据 EcoSimR 包的方法，计算 C-score，需转为 0-1 丰度矩阵
comm[comm>0] <- 1

#按分组拆分为 ENV1 和 ENV2 两组数据
comm_ENV1 <- comm[ ,c(1,5,7)]
comm_ENV2 <- comm[ ,c(9,11,13)]
comm_ENV3 <- comm[ ,c(3,15)]

comm_ENV1 <- comm[ ,c(2,6,8)]
comm_ENV2 <- comm[ ,c(10,12,14)]
comm_ENV3 <- comm[ ,c(3,15)]




#计算基于 C-score 的零模型，详情 ?cooc_null_model
#本示例基于 1000 次随机置换模拟获得 C-score 零分布（实际使用时可能要设置高一些，如 Mo 等（2021）的 30000 次随机，不过可能会比较慢）
set.seed(123)
Cscore_ENV1 <- cooc_null_model(comm_ENV1, algo = 'sim9', metric = 'c_score', nReps = 1000, 
                               saveSeed = FALSE, burn_in = 500, algoOpts = list(), metricOpts = list(), suppressProg = FALSE)

summary(Cscore_ENV1)  #ENV1 群落的 C-score 和 SES
plot(Cscore_ENV1, type = 'hist')  #C-score 观测值（经验 C-score）和 C-score 模拟值（C-score 零分布）的比较

set.seed(123)
Cscore_ENV2 <- cooc_null_model(comm_ENV2, algo = 'sim9', metric = 'c_score', nReps = 1000, 
                               saveSeed = FALSE, burn_in = 500, algoOpts = list(), metricOpts = list(), suppressProg = FALSE)

summary(Cscore_ENV2)  #ENV2 群落的 C-score 和 SES
plot(Cscore_ENV2, type = 'hist')  #C-score 观测值（经验 C-score）和 C-score 模拟值（C-score 零分布）的比较

set.seed(123)
Cscore_ENV3 <- cooc_null_model(comm_ENV3, algo = 'sim9', metric = 'c_score', nReps = 1000, 
                               saveSeed = FALSE, burn_in = 500, algoOpts = list(), metricOpts = list(), suppressProg = FALSE)

summary(Cscore_ENV3)  #ENV2 群落的 C-score 和 SES
plot(Cscore_ENV3, type = 'hist')  #C-score 观测值（经验 C-score）和 C-score 模拟值（C-score 零分布）的比较




#根据上述统计好的 C-score 观测值（经验 C-score）、C-score 模拟值（C-score 零分布）以及 SES 等，作个手动记录输出
result <- data.frame(group = c('metal30', 'organic30'), 
                     Cscore_obs = c(0.36216, 0.34826),  #C-score 观测值（经验 C-score）
                     Cscore_sim = c(0.36316 , 0.35032),  #C-score 模拟值（C-score 零分布）的均值
                     SES = c(-0.75336 , -6.2226)  #标准化效应大小（SES）
)
result <- data.frame(group = c('metal60', 'organic60'), 
                     Cscore_obs = c(0.32267, 0.43735),  #C-score 观测值（经验 C-score）
                     Cscore_sim = c(0.32223, 0.44268),  #C-score 模拟值（C-score 零分布）的均值
                     SES = c(0.27701, -4.7453)  #标准化效应大小（SES）
)
write.table(result, 'result.txt', row.names = FALSE, sep = '\t', quote = FALSE)

#仿照 Mo 等（2021）的可视化方法，绘制三 Y 坐标图展示 C-score 观测值（经验 C-score）、C-score 模拟值（C-score 零分布）以及 SES
library(ggplot2)

p.Cscore_obs <- ggplot(data = result, aes(x = group, y = Cscore_obs)) +
  geom_col(fill = 'gray30', width = 0.6) + 
  theme(panel.grid = element_blank(), panel.background = element_rect(fill = NA, color = 'gray30')) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, max(result$Cscore_obs, result$Cscore_sim)*1.2)) +
  labs(x = '', y = 'C-score obs')

p.Cscore_obs

p.Cscore_sim <- ggplot(data = result, aes(x = group, y = Cscore_sim)) +
  geom_col(fill = 'blue', width = 0.4) + 
  theme(panel.grid = element_blank(), panel.background = element_rect(fill = NA, color = NA),  
        axis.text.y = element_text(color = 'blue'), axis.ticks.y = element_line(color = 'blue'), 
        axis.title.y = element_text(color = 'blue'), axis.line.y = element_line(color = 'blue')) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, max(result$Cscore_obs, result$Cscore_sim)*1.2)) +
  labs(x = '', y = 'C-score sim\n')

p.Cscore_sim

p.SES <- ggplot(data = result, aes(x = group, y = SES)) +
  geom_point(color = 'red', shape = 15, size = 5) + 
  theme(panel.grid = element_blank(), panel.background = element_rect(fill = NA, color = NA),  
        axis.text.y = element_text(color = 'red'), axis.ticks.y = element_line(color = 'red'), 
        axis.title.y = element_text(color = 'red'), axis.line.y = element_line(color = 'red')) +
  labs(x = '', y = 'Standardized Effect Size (SES)')

p.SES

#备注：自定义函数 y3_plot() 的来源：https://mp.weixin.qq.com/s/Wl01G8_6-e0GgBLnbrK74A
library(ggplot2)
library(gtable)
library(grid)
y3_plot <- function(gp1, gp2, gp3) {
  p1 <- ggplotGrob(gp1)
  p2 <- ggplotGrob(gp2)
  p3 <- ggplotGrob(gp3)
  
  # Get the location of the plot panel in p1.
  # These are used later when transformed elements of p2 are put back into p1
  pp <- c(subset(p1$layout, name == 'panel', se = t:r))
  
  # Overlap panel for second plot on that of the first plot
  p1 <- gtable_add_grob(p1, p2$grobs[[which(p2$layout$name == 'panel')]], pp$t, pp$l, pp$b, pp$l)
  p1 <- gtable_add_grob(p1, p3$grobs[[which(p3$layout$name == 'panel')]], pp$t, pp$l, pp$b*0.0001, pp$l)
  
  # Then proceed as before:
  
  # ggplot contains many labels that are themselves complex grob; 
  # usually a text grob surrounded by margins.
  # When moving the grobs from, say, the left to the right of a plot,
  # Make sure the margins and the justifications are swapped around.
  # The function below does the swapping.
  # Taken from the cowplot package:
  # https://github.com/wilkelab/cowplot/blob/master/R/switch_axis.R 
  
  hinvert_title_grob <- function(grob){
    
    # Swap the widths
    widths <- grob$widths
    grob$widths[1] <- widths[3]
    grob$widths[3] <- widths[1]
    grob$vp[[1]]$layout$widths[1] <- widths[3]
    grob$vp[[1]]$layout$widths[3] <- widths[1]
    
    # Fix the justification
    grob$children[[1]]$hjust <- 1 - grob$children[[1]]$hjust 
    grob$children[[1]]$vjust <- 1 - grob$children[[1]]$vjust 
    grob$children[[1]]$x <- unit(1, 'npc') - grob$children[[1]]$x
    grob
  }
  
  #### add p3
  # Get the y axis title from p3
  index <- which(p3$layout$name == 'ylab-l') # Which grob contains the y axis title?
  ylab <- p3$grobs[[index]]                # Extract that grob
  ylab <- hinvert_title_grob(ylab)         # Swap margins and fix justifications
  
  # Put the transformed label on the right side of p1
  p1 <- gtable_add_cols(p1, p3$widths[p3$layout[index, ]$l], pp$r)
  p1 <- gtable_add_grob(p1, ylab, pp$t, pp$r + 1, pp$b, pp$r + 1, clip = 'off', name = 'ylab-r')
  
  # Get the y axis from p3 (axis line, tick marks, and tick mark labels)
  index <- which(p3$layout$name == 'axis-l')  # Which grob
  yaxis <- p3$grobs[[index]]                  # Extract the grob
  
  # yaxis is a complex of grobs containing the axis line, the tick marks, and the tick mark labels.
  # The relevant grobs are contained in axis$children:
  #   axis$children[[1]] contains the axis line;
  #   axis$children[[2]] contains the tick marks and tick mark labels.
  
  # First, move the axis line to the left
  yaxis$children[[1]]$x <- unit.c(unit(0, 'npc'), unit(0, 'npc'))
  
  # Second, swap tick marks and tick mark labels
  ticks <- yaxis$children[[2]]
  ticks$widths <- rev(ticks$widths)
  ticks$grobs <- rev(ticks$grobs)
  
  # Third, move the tick marks
  ticks$grobs[[1]]$x <- ticks$grobs[[1]]$x - unit(1, 'npc') + unit(3, 'pt')
  
  # Fourth, swap margins and fix justifications for the tick mark labels
  ticks$grobs[[2]] <- hinvert_title_grob(ticks$grobs[[2]])
  
  # Fifth, put ticks back into yaxis
  yaxis$children[[2]] <- ticks
  
  # Put the transformed yaxis on the right side of p1
  p1 <- gtable_add_cols(p1, p3$widths[p3$layout[index, ]$l], pp$r)
  p1 <- gtable_add_grob(p1, yaxis, pp$t, pp$r + 1, pp$b, pp$r + 1, clip = 'off', name = 'axis-r')
  
  #### add p2
  # Get the y axis title from p2
  index <- which(p2$layout$name == 'ylab-l') # Which grob contains the y axis title?
  ylab <- p2$grobs[[index]]                # Extract that grob
  ylab <- hinvert_title_grob(ylab)         # Swap margins and fix justifications
  
  # Put the transformed label on the right side of p1
  p1 <- gtable_add_cols(p1, p2$widths[p2$layout[index, ]$l], pp$r)
  p1 <- gtable_add_grob(p1, ylab, pp$t, pp$r + 1, pp$b, pp$r + 1, clip = 'off', name = 'ylab-r')
  
  # Get the y axis from p2 (axis line, tick marks, and tick mark labels)
  index <- which(p2$layout$name == 'axis-l')  # Which grob
  yaxis <- p2$grobs[[index]]                  # Extract the grob
  
  # yaxis is a complex of grobs containing the axis line, the tick marks, and the tick mark labels.
  # The relevant grobs are contained in axis$children:
  #   axis$children[[1]] contains the axis line;
  #   axis$children[[2]] contains the tick marks and tick mark labels.
  
  # First, move the axis line to the left
  yaxis$children[[1]]$x <- unit.c(unit(0, 'npc'), unit(0, 'npc'))
  
  # Second, swap tick marks and tick mark labels
  ticks <- yaxis$children[[2]]
  ticks$widths <- rev(ticks$widths)
  ticks$grobs <- rev(ticks$grobs)
  
  # Third, move the tick marks
  ticks$grobs[[1]]$x <- ticks$grobs[[1]]$x - unit(1, 'npc') + unit(3, 'pt')
  
  # Fourth, swap margins and fix justifications for the tick mark labels
  ticks$grobs[[2]] <- hinvert_title_grob(ticks$grobs[[2]])
  
  # Fifth, put ticks back into yaxis
  yaxis$children[[2]] <- ticks
  
  # Put the transformed yaxis on the right side of p1
  p1 <- gtable_add_cols(p1, p2$widths[p2$layout[index, ]$l], pp$r)
  p1 <- gtable_add_grob(p1, yaxis, pp$t, pp$r + 1, pp$b, pp$r + 1, clip = 'off', name = 'axis-r')
  grid.newpage()
  grid.draw(p1)
}
y3_plot <- function(gp1, gp2, gp3) {
  p1 <- ggplotGrob(gp1)
  p2 <- ggplotGrob(gp2)
  p3 <- ggplotGrob(gp3)
  
  # Get the location of the plot panel in p1.
  pp <- c(subset(p1$layout, name == 'panel', se = t:r))
  
  # Overlap panel for second plot on that of the first plot
  p1 <- gtable_add_grob(p1, p2$grobs[[which(p2$layout$name == 'panel')]], pp$t, pp$l, pp$b, pp$l)
  p1 <- gtable_add_grob(p1, p3$grobs[[which(p3$layout$name == 'panel')]], pp$t, pp$l, pp$b*0.0001, pp$l)
  
  # Then proceed as before:
  # ...
  
  # Finally, draw the plot
  grid.newpage()
  grid.draw(p1)
} 
y3_plot(gp1 = p.Cscore_obs, gp2 = p.Cscore_sim, gp3 = p.SES)
#自定义函数，用于组合 ggplot2 绘图结果构建三个坐标轴
#改编自先前绘制双坐标轴的一个自定义函数，原双坐标轴函数可见：https://mp.weixin.qq.com/s/DnotJG4ic8fqFwikkf5QKQ

#Specificity-Occupancy，SPEC-OCCU分析####
#读取 OTU 丰度表和样本分组
##得到属级别的丰度表
men <- select(asv,-col1,-col6,-col3,-col4,-col5,-col7)
colnames(men)[17] <- 'Taxonomy'
colnames(men)[1] <- 'OTU'
otu <- men
rownames(otu) <- men$OTU
otu <- otu[,-1]
otu <- na.omit(otu)

group <- metadata
rownames(group) <- metadata$sample

#30d和60d分开做
#30d
#otu <- otu[,c(1,5,7,9,11,13,16)]
#group <- group[c(1,5,7,9,11,13),]
#60d
#otu <- otu[,c(2,6,8,10,12,14,16)]
#group <- group[c(2,6,8,10,12,14),]

##计算各组样本中，OTU 的特异性（specificity）和占有率（occupancy）

Nindividuals_S <- rep(0, nrow(otu))
for (i in unique(group$category)) {
  otu_group_i <- otu[ ,rownames(subset(group, category == i))]
  Nindividuals_S <- Nindividuals_S + rowMeans(otu_group_i)  #计算 Nindividuals S
}

spec_occu <- NULL
for (i in unique(group$category)) {
  otu_group_i <- otu[ ,rownames(subset(group, category == i))]
  Nindividuals_SH <- apply(otu_group_i, 1, mean)  #计算 Nindividuals SH
  Specificity <- Nindividuals_SH / Nindividuals_S  #计算 Specificity
  Nsites_H <- ncol(otu_group_i)  #计算 Nsites H
  Nsites_SH <- apply(otu_group_i, 1, function(x) sum(x>0))  #计算 Nsites SH
  Occupancy <- Nsites_SH / Nsites_H  #计算 Occupancy
  spec_occu_group_i <- data.frame(category = i, OTU = rownames(otu_group_i), Specificity = Specificity, Occupancy = Occupancy, Abundance_mean = rowMeans(otu_group_i), Taxonomy = otu$Taxonomy)
  spec_occu <- rbind(spec_occu, spec_occu_group_i)  #合并各组统计
}
head(spec_occu)  #该数据框包含各组中各 OTU 的名称、特异性、占有率、平均丰度、类群信息等

#输出统计表格
#write.table(spec_occu, 'spec_occu.txt', row.names = FALSE, sep = '\t', quote = FALSE)

#绘制 SPEC-OCCU 图
library(ggplot2)

spec_occu <- na.omit(spec_occu)


p <- ggplot(spec_occu, aes(Occupancy, Specificity)) +
  geom_point(aes(size = log10(Abundance_mean), color = Taxonomy)) +
  geom_jitter(aes(size = log10(Abundance_mean), color = Taxonomy)) +  #如果觉得普通散点图中的点重叠严重不好看，可以仿照 Gweon et al (2021) 使用抖动点图来展示
  scale_size(breaks = c(-1, -2, -3, -4), labels = c(expression(10^{-1}), expression(10^{-2}), expression(10^{-3}), expression(10^{-4})), range = c(0, 4)) +
  scale_color_manual(values = c('#F59D1F', '#768CC5', '#9BC648', '#794779', '#A19E9E'), limits = c('p:Proteobacteria', 'p:Nitrospirota', 'p:Gemmatimonadota', 'p:Chloroflexi', 'Others')) +
  theme(panel.grid = element_blank(), panel.background = element_rect(fill = 'white', color = 'gray30'), legend.key = element_blank()) +
  facet_wrap(~category, ncol = 3) +
  scale_x_continuous(breaks = c(0, 0.5, 1), expand = c(0, 0), limit = c(0, 1)) +
  scale_y_continuous(breaks = c(0, 0.5, 1), expand = c(0, 0), limit = c(0, 1)) +
  labs(x = 'Occupancy', y = 'Specificity', size = 'Average relative abundance', color = 'Taxonomy') +
  coord_cartesian(clip = 'off')

p
##识别各组样本中的特化种（specialist species）

#在上述统计表格中直接根据特异性和占有率 ≥0.7 做筛选，保留下的 OTU 即为特化种
spec_occu_specialist <- subset(spec_occu, Specificity >= 0.7 & Occupancy >= 0.7)
head(spec_occu_specialist)

#输出统计表格
write.table(spec_occu_specialist, 'spec_occu_specialist_all.txt', row.names = FALSE, sep = '\t', quote = FALSE)

#在上述 SPEC-OCCU 图中添加阈值线
p + 
  geom_segment(aes(x = 0.7, xend = 1, y = 0.7, yend = 0.7), linetype = 2)+
  geom_segment(aes(x = 0.7, xend = 0.7, y = 0.7, yend = 1), linetype = 2)
#RDA或者procurtis/VPA####

#SVM####


