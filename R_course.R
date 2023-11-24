####R语言基础####
# 加载包
library(tidyverse)

#R镜像 快速选择
#chooseBioCmirror() #这样下载会更快

# 设置目录
setwd("Day1/Day1.1/Day1.1.1")

# 读取txt文件
c <- read.table("x.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)

# 复制一份文件并重命名
a <- c

# 输出txt文件
write.table(a,"a.txt",sep = "\t",row.names = T,col.names = NA,quote = F)

# 输出csv文件
write.csv(a, file="a.csv")

# 读取csv文件
d <- read.csv("a.csv", row.names = 1)

# 数据框行列反转
x <- t(x)

# 变成数据框
x <- as.data.frame(x)

# class判断数据类型
class(x)

# $提取数据框中的一列
class(x$`TCGA-DD-AACK-01A`)

# substr 提取
substr("stellalxy", 1, 4)

# c()创建集合
jihe <- c("a","b","c")
jihe

# 提取/删除数据框行列（切片）
q = x[, 1:3]
p = x[1:3,]
r = x[1:3, 5:10]
s = x[-1,] # 去掉第一行
t = x[,-1]
u = x[,-(2:4)]
v <- x[,-c(1,3)] #去掉第一列和第三列

# 运用传导符%>% shift+ctrl+M
x <- x %>% t() %>% as.data.frame() # 传导符让很多指令串在一起执行

# duplicated函数判断是否重复
a <- c("a","b","a","b","c")
duplicated(a)
# 将结果反过来
!duplicated(a)
# 将判断结果为TRUE的字符提取出来
a <- a[!duplicated(a)]

#inner_join
#tribble创建简易数据框
class1 <- tribble(
  ~'名次',~'姓名',
  '第一名','王某人',
  '第二名','张某人',
  '第三名','李某人'
) # ～代表列名
class2 <- tribble(
  ~'名次',~'姓名',
  '第一名','胡某人',
  '第二名','刘某人',
  '第四名','于某人'
)
# 合并数据框（合并共有的列）
inner_join(class1, class2, by="名次")
# 合并左边的数据框的列
left_join(class1, class2, by="名次")
# 合并右边的数据框的列
right_join(class1, class2, by="名次")


####TCGA-LUAD数据下载####
setwd("../../../TCGA-LUAD")
setwd("TCGAdata")
library(tidyverse)

# 安装TCGAbiolinks包
#install.packages("BiocManager")
library(BiocManager)
#BiocManager::install("TCGAbiolinksGUI.data")
#BiocManager::install("remotes",force = TRUE)
#BiocManager::install("ExperimentHub",force=TRUE)
#BiocManager::install("TCGAbiolinks")
library(TCGAbiolinks)
cancer_type = "TCGA-LUAD" # 癌症类型
# TCGA肿瘤缩写：https://www.jianshu.com/p/3c0f74e85825
# 下载肺癌数据
expquery <- GDCquery(project = cancer_type,
                     data.category = "Transcriptome Profiling",
                     data.type = "Gene Expression Quantification",
                     workflow.type = "STAR - Counts")
GDCdownload(expquery, directory = "GDCdata")
expquery2 <- GDCprepare(expquery, directory = "GDCdata",summarizedExperiment = T)
save(expquery2, file="luad.gdc_2022.rda")

setwd("TCGA-LUAD")
setwd("TCGAdata")
library(tidyverse)
# 加载rda文件
load("luad.gdc_2022.rda")
load("gene_annotation_2022.rda")
# table分组计数
table(gene_annotation_2022$type)
# TCGA数据库表达谱形式有三种：counts，FPKM，TPMS
# 我们用counts和TPMS
# TCGA中counts只用来做差异分析，真正的基因表达谱是TPMS
# 提取counts
counts <- expquery2@assays@data@listData[["unstranded"]]
# counts为基因表达谱，里面都是整数，每一列为一个样本，每一行都是一个基因
# 为counts数据框添加行名与列名
colnames(counts) <- expquery2@colData@rownames
rownames(counts) <- expquery2@rowRanges@ranges@NAMES
# 拆解以上代码
counts1 <- expquery2@assays@data@listData[["unstranded"]]
colnames(counts1) <- expquery2@colData@rownames
rownames(counts1) <- expquery2@rowRanges@ranges@NAMES
# 将counts1转换为数据框
counts1 <- as.data.frame(counts1)
# 将counts1中的行名变成数据框中的一个列
counts1 <- rownames_to_column(counts1, var="ENSEMBL")
# 以上这么做的目的是为了实施inner_join
counts1 <- inner_join(counts1, gene_annotation_2022, "ENSEMBL")
# 对counts1的symbol列去重复
counts1 <- counts1[!duplicated(counts1$symbol),]
# 对counts实施相同的操作
counts <- counts %>% as.data.frame() %>% rownames_to_column("ENSEMBL") %>% 
  inner_join(gene_annotation_2022, "ENSEMBL") %>%
  .[!duplicated(.$symbol),]

# 大致判断counts和counts1是否相同
a <- c("a", "b", "a", "b", "c")
b <- c("a", "b", "b", "a", "c")
identical(a, b)
identical(colnames(counts), colnames(counts1))
# counts与counts1列名相同
identical(rownames(counts), rownames(counts1))
# counts与counts1行名相同
# R语言中数据框的行名不能重复
# 去除counts1中的行名
rownames(counts1) <- NULL
# 将symbol列变为counts1的行名
counts1 <- column_to_rownames(counts1, var = "symbol")
# counts一样操作
rownames(counts) <- NULL
counts <- column_to_rownames(counts, var="symbol")

# table查看基因类型数量
table(counts$type)
# 将counts中type为protein coding的行都提取出来
counts <- counts[counts$type == "protein_coding",]
# 将counts中的第一列与最后一列去除
counts <- counts[, -c(1, ncol(counts))]
# ncol计算数据框中有多少列
# 把TCGA barcode切割为16位字符，并去除重复样本
colnames(counts) <- substr(colnames(counts),1,16)
counts <- counts[,!duplicated(colnames(counts))]
table(substr(colnames(counts),14,16))
# 01A是肿瘤样本，11A是正常样本
# 提取01A样本
counts01A <- counts[,substring(colnames(counts),14,16)=="01A"]
# 提取11A样本
counts11A <- counts[,substring(colnames(counts),14,16)=="11A"]


####tpms####
# 提取tpms数据
tpms <- expquery2@assays@data@listData[["tpm_unstrand"]]
colnames(tpms) <- expquery2@colData@rownames
rownames(tpms) <- expquery2@rowRanges@ranges@NAMES
tpms <- tpms %>% 
  as.data.frame() %>% 
  rownames_to_column("ENSEMBL") %>% 
  inner_join(gene_annotation_2022,"ENSEMBL") %>% 
  .[!duplicated(.$symbol),]
rownames(tpms) <- NULL
tpms <- column_to_rownames(tpms,"symbol")
tpms <- tpms[tpms$type=="protein_coding",]
tpms <- tpms[,-c(1, ncol(tpms))]
colnames(tpms) <- substring(colnames(tpms),1,16)
tpms <- tpms[,!duplicated(colnames(tpms))]
tpms01A <- tpms[,substr(colnames(tpms),14,16)=="01A"]
tpms11A <- tpms[,substr(colnames(tpms),14,16)=="11A"]

# 判断counts和tpms的行列名是否一致
identical(rownames(counts01A), rownames(counts11A))
identical(rownames(tpms01A), rownames(tpms11A))
identical(rownames(counts01A), rownames(tpms01A))
identical(colnames(counts01A), colnames(tpms01A))
# 保存counts和tpms数据
write.table(counts01A,"counts01A.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
write.table(counts11A,"counts11A.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
write.table(tpms01A,"tpms01A.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
write.table(tpms11A,"tpms11A.txt",sep = "\t",row.names = T,col.names = NA,quote = F)

# cbind为按列合并，rbind为按行合并
counts <- cbind(counts01A, counts11A)
tpms <- cbind(tpms01A, tpms11A)
write.table(counts,"counts.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
write.table(tpms,"tpms.txt",sep = "\t",row.names = T,col.names = NA,quote = F)


####tpms_log2####
range(tpms) #查看数据范围
range(tpms01A)
range(tpms11A)
# 数据范围太大，需要对数据进行等比例缩小，使用log2
tpms_log2 <- log2(tpms+1)
range(tpms_log2)
tpms01A_log2 <- log2(tpms01A+1)
range(tpms01A_log2)
tpms11A_log2 <- log2(tpms11A+1)
range(tpms11A_log2)
#保存
write.table(tpms_log2,"tpms_log2.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
write.table(tpms01A_log2,"tpms01A_log2.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
write.table(tpms11A_log2,"tpms11A_log2.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
# 表达谱整理完毕


####ESTIMATE####
# 计算患者免疫评分与肿瘤纯度
# 肿瘤组织中含有基质组分，免疫组分以及肿瘤组分
# 基质组分和免疫组分占比越大，肿瘤组分占比越小
# 使用estimate计算肿瘤纯度
setwd("../")
setwd("ESTIMATE")
#安装包
library(utils)
#install.packages("estimate", repos="http://R-Forge.R-project.org")
library(tidyverse)
library(estimate)

# 读取肿瘤患者01A的表达谱
exp <- read.table("tpms01A_log2.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
# 计算免疫评分
# 筛选共有基因
filterCommonGenes(input.f="tpms01A_log2.txt", output.f = "tpms01A_log2.gct", id="GeneSymbol")
estimateScore("tpms01A_log2.gct", "tpms01A_log2_estimate_score.txt", platform = "affymetrix")

# 提取结果并整理
ESTIMATE_result = read.table("tpms01A_log2_estimate_score.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
# 去除第一列
ESTIMATE_result <- ESTIMATE_result[,-1]
# 将第一行变成列名（第一行为样本名称）
colnames(ESTIMATE_result) <- ESTIMATE_result[1,]
ESTIMATE_result <- as.data.frame(t(ESTIMATE_result[-1,]))
rownames(ESTIMATE_result) <- colnames(exp)
# 保存
write.table(ESTIMATE_result, file = "ESTIMATE_result.txt",sep = "\t",row.names = T,col.names = NA,quote = F)


####生存信息整理####
#xena官网：https://xenabrowser.net/datapages/
# 下载生存信息
setwd("../")
setwd("Survival_data")
library(tidyverse)
# 手动导入OS.txt取名survival
# 提取OS和OS.time
survival <- survival[,2:3]
survival <- survival %>% rownames_to_column("sample")
survival$name <- paste0(survival$sample, "A")
table(substring(survival$name,14,16))
rownames(survival) <- survival$name
survival <- survival[,2:3]
# OS表示生存状态，0表示生存，1表示死亡
# 合并生存信息与表达谱
tpms01A_log2 <- read.table("tpms01A_log2.txt", sep = "\t",row.names = 1,check.names = F,header = T)
# 取表达谱和生存信息中样本名称的交集
a <- intersect(colnames(tpms01A_log2), rownames(survival))
table(substr(a,14,16))
exp_01A <- tpms01A_log2[,a]
surv_01A <- survival[a,]
# 将exp01A进行行列转换
exp_01A <- exp_01A %>% t() %>% as.data.frame()
identical(rownames(exp_01A), rownames(surv_01A))
# 将表达谱与生存信息合并为一个数据框
exp_surv_01A <- cbind(surv_01A, exp_01A)
# 保存
write.table(exp_surv_01A,"exp_surv_01A.txt",sep = "\t",row.names = T,col.names = NA,quote = F)

#合并生存信息与ESTIMATE
ESTIMATE_result <- read.table("ESTIMATE_result.txt", sep = "\t",row.names = 1,check.names = F,header = T)
identical(rownames(ESTIMATE_result), rownames(surv_01A))
ESTIMATE_result_surv_01A <- cbind(surv_01A, ESTIMATE_result)
# 保存
write.table(ESTIMATE_result_surv_01A,"ESTIMATE_result_surv_01A.txt",sep = "\t",row.names = T,col.names = NA,quote = F)


####根据ESTIMATE_result高低组做生存分析####
setwd("../")
setwd("survival")
surv <- read.table("ESTIMATE_result_surv_01A.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
# 将OS.time的单位从天变成年
surv$OS.time <- surv$OS.time/365

# median 中位数
# 找到immune score的中位数
# 根据immune score的中位数把样本分成high和low两组
surv$group <- ifelse(surv$ImmuneScore > median(surv$ImmuneScore), "High", "Low")
# 将High和Low从字符型转变为因子型数据
surv$group <- factor(surv$group, levels=c("Low", "High"))
class(surv$group)
table(surv$group)

# 生存分析
#install.packages("survival")
library(survival)
fitd <- survdiff(Surv(OS.time, OS) ~ group,
                 data = surv,
                 na.action = na.exclude)
pValue <- 1 - pchisq(fitd$chisq, length(fitd$n) - 1)

# 拟合生存曲线
fit <- survfit(Surv(OS.time, OS) ~ group, data=surv)
summary(fit)
p.lab <- paste0("P", ifelse(pValue < 0.001, "< 0.001", paste0(" = ", round(pValue, 3))))

# 画图
#install.packages("survminer")
library(survminer)
ggsurvplot(fit,
           data = surv,
           pval = p.lab,
           conf.int = TRUE, # 显示置信区间
           risk.table = TRUE, # 显示风险表
           risk.table.col = "strata",
           palette = "jco", # 配色采用jco
           legend.labs = c("Low", "High"), # 图例
           size = 1,
           xlim = c(0,20), # x轴长度
           break.time.by = 5, # x轴步长为5
           legend.title = "ImmuneScore",
           surv.median.line = "hv", # 限制垂直和水平的中位生存
           ylab = "Survival probability (%)", # 修改y轴标签
           xlab = "Time (Years)", # 修改x轴标签
           ncensor.plot = TRUE, # 显示删失图块
           ncensor.plot.height = 0.25,
           risk.table.y.text = FALSE)
dev.off()

# 找到stromal score的中位数
surv$group <- ifelse(surv$StromalScore > median(surv$StromalScore), "High", "Low")
surv$group <- factor(surv$group, levels = c("Low", "High"))
class(surv$group)
table(surv$group)
#生存分析
fitd <- survdiff(Surv(OS.time, OS) ~ group,
                 data = surv,
                 na.action = na.exclude)
pValue <- 1 - pchisq(fitd$chisq, length(fitd$n) - 1)
# 拟合生存曲线
fit <- survfit(Surv(OS.time, OS) ~ group, data=surv)
summary(fit)
p.lab <- paste0("P", ifelse(pValue < 0.001, "< 0.001", paste0(" = ", round(pValue, 3))))
# 画图
ggsurvplot(fit,
           data = surv,
           pval = p.lab,
           conf.int = TRUE, # 显示置信区间
           risk.table = TRUE, # 显示风险表
           risk.table.col = "strata",
           palette = "jco", # 配色采用jco
           legend.labs = c("Low", "High"), # 图例
           size = 1,
           xlim = c(0,20), # x轴长度
           break.time.by = 5, # x轴步长为5
           legend.title = "StromalScore",
           surv.median.line = "hv", # 限制垂直和水平的中位生存
           ylab = "Survival probability (%)", # 修改y轴标签
           xlab = "Time (Years)", # 修改x轴标签
           ncensor.plot = TRUE, # 显示删失图块
           ncensor.plot.height = 0.25,
           risk.table.y.text = FALSE)
dev.off()

# ESTIMATE score生存分析
surv$group <- ifelse(surv$ESTIMATEScore > median(surv$ESTIMATEScore), "High", "Low")
surv$group <- factor(surv$group, levels = c("Low", "High"))
class(surv$group)
table(surv$group)
#生存分析
fitd <- survdiff(Surv(OS.time, OS) ~ group,
                 data = surv,
                 na.action = na.exclude)
pValue <- 1 - pchisq(fitd$chisq, length(fitd$n) - 1)
# 拟合生存曲线
fit <- survfit(Surv(OS.time, OS) ~ group, data=surv)
summary(fit)
p.lab <- paste0("P", ifelse(pValue < 0.001, "< 0.001", paste0(" = ", round(pValue, 3))))
# 画图
ggsurvplot(fit,
           data = surv,
           pval = p.lab,
           conf.int = TRUE, # 显示置信区间
           risk.table = TRUE, # 显示风险表
           risk.table.col = "strata",
           palette = "jco", # 配色采用jco
           legend.labs = c("Low", "High"), # 图例
           size = 1,
           xlim = c(0,20), # x轴长度
           break.time.by = 5, # x轴步长为5
           legend.title = "ESTIMATEScore",
           surv.median.line = "hv", # 限制垂直和水平的中位生存
           ylab = "Survival probability (%)", # 修改y轴标签
           xlab = "Time (Years)", # 修改x轴标签
           ncensor.plot = TRUE, # 显示删失图块
           ncensor.plot.height = 0.25,
           risk.table.y.text = FALSE)
dev.off()


####整理TCGA临床信息####
setwd("../")
setwd("clinical")
library(tidyverse)
load("luad.gdc_2022.rda")
# 提取临床信息
clinical <- as.data.frame(expquery2@colData) %>% .[!duplicated(.$sample),]
# 提取需要的临床信息数据
clinical <- clinical[, c("gender", "age_at_index", "ajcc_pathologic_stage", "ajcc_pathologic_t","ajcc_pathologic_n","ajcc_pathologic_m")]

class(clinical$gender)
class(clinical$age_at_index)
class(clinical$ajcc_pathologic_stage)
class(clinical$ajcc_pathologic_t)
class(clinical$ajcc_pathologic_n)
class(clinical$ajcc_pathologic_m)

table(clinical$gender)
table(clinical$age_at_index)
table(clinical$ajcc_pathologic_stage)
table(clinical$ajcc_pathologic_t)
table(clinical$ajcc_pathologic_n)
table(clinical$ajcc_pathologic_m)

# 把stage列里的A和B都去掉
clinical$ajcc_pathologic_stage <- gsub("A", "", clinical$ajcc_pathologic_stage)
clinical$ajcc_pathologic_stage <- gsub("B", "", clinical$ajcc_pathologic_stage)

# 把T分期列里的a和b都去掉
clinical$ajcc_pathologic_t <- gsub("a", "", clinical$ajcc_pathologic_t)
clinical$ajcc_pathologic_t <- gsub("b", "", clinical$ajcc_pathologic_t)

# 把m分期里的a和b去掉
clinical$ajcc_pathologic_m <- gsub("a", "", clinical$ajcc_pathologic_m)
clinical$ajcc_pathologic_m <- gsub("b", "", clinical$ajcc_pathologic_m)

# 提取01A临床数据
rownames(clinical) <- substring(rownames(clinical),1,16)

# 将基因表达谱和临床数据合并并保存
exp01A <- read.table("tpms01A_log2.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
# 在clinical中提取表达谱中的513个样本数据
clinical01A <- clinical[colnames(exp01A),]
exp01A <- exp01A %>% t() %>% as.data.frame()
identical(rownames(exp01A), rownames(clinical01A))
clinical.expr01A <- cbind(clinical01A, exp01A)
write.table(clinical.expr01A,"clinical.expr01A.txt",sep = "\t",row.names = T,col.names = NA,quote = F)

# 将ESTIMATE_result和临床信息合并
ESTIMATE_result <- read.table("ESTIMATE_result.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
identical(rownames(ESTIMATE_result), rownames(clinical01A))
clinical.ESTIMATE_result01A <- cbind(clinical01A, ESTIMATE_result)
write.csv(clinical.ESTIMATE_result01A,file = "clinical.ESTIMATE_result01A.csv")
#解螺旋：https://www.helixlife.cn/class


####差异分析####
# 免疫评分
setwd("../")
setwd("Immune_DEG")
library(BiocManager)
# 下载做差异分析的包
#BiocManager::install("DESeq2")
library(DESeq2)
library(tidyverse)
# TCGA差异分析用counts来做，因为把01A患者分组做差异分析所以读取01A
counts_01A <- read.table("counts01A.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
# 因为是用免疫评分分组所以读取ESTIMATE_result
estimate <- read.table("ESTIMATE_result.txt", sep = "\t",row.names = 1,check.names = F,header = T)
# 整理分组信息
x <- "ImmuneScore"
# 求immune score列的中位数
med <- as.numeric(median(estimate[,x]))
estimate <- as.data.frame(t(estimate))
identical(colnames(estimate), colnames(counts_01A))

# 高于免疫评分中位数的为high，低于为low
conditions = data.frame(sample=colnames(counts_01A),
           group=factor(ifelse(estimate[x,]>med, "high", "low"),levels=c("low","high")))
conditions <- column_to_rownames(conditions, "sample")

# 差异分析准备工作
dds <- DESeqDataSetFromMatrix(countData = counts_01A,
                              colData = conditions,
                              design = ~ group)
# 开始差异分析
dds <- DESeq(dds)
resultsNames(dds)
# 提取结果
res <- results(dds)
# 保存结果
save(res, file="DEG_ImmuneScore.Rda")


####热图绘制####
DEG <- as.data.frame(res)
# DEG中数值为正的基因是高表达的，负的是低表达的
# 高表达为在high组里高表达，在low组里低表达
# 读取表达谱
exp <- read.table("tpms01A_log2.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
# 添加上下调信息
logFC_cutoff <- 1
# 两个筛选条件
type1 = (DEG$padj < 0.05)&(DEG$log2FoldChange < -logFC_cutoff)
type2 = (DEG$padj < 0.05)&(DEG$log2FoldChange > logFC_cutoff)
# 根据上调或下调倍数对基因进行分类
DEG$change <- ifelse(type1,"DOWN",ifelse(type2, "UP", "NOT"))
table(DEG$change)

# 下载pheatmap包
#install.packages('pheatmap')
library(dplyr)
library(pheatmap)
# 提取差异基因表达谱
a <- filter(DEG, change == "UP")
b <- filter(DEG, change == "DOWN")
c <- rbind(a,b)
d <- rownames(c)
exp_diff <- exp[d,]
# 根据exp_diff绘制热图
# 设置分组信息
annotation_col <- conditions
# 对exp_diff列的顺序进行处理
a <- filter(annotation_col, group == "high")
b <- filter(annotation_col, group == "low")
exp_diff_high <- exp_diff[,rownames(a)]
exp_diff_low <- exp_diff[,rownames(b)]
exp_diff <- cbind(exp_diff_high, exp_diff_low)
# 画图
pheatmap(exp_diff,
         annotation_col=annotation_col,
         scale = "row",
         show_rownames = F,
         show_colnames =F,
         color = colorRampPalette(c("navy", "white", "red"))(50),
         cluster_cols =F,
         cluster_rows = T,
         fontsize = 10,
         fontsize_row=3,
         fontsize_col=3)
dev.off()

####基质评分####
setwd("../")
setwd("Stromal_DEG")
library(BiocManager)
library(DESeq2)
library(tidyverse)
#TCGA差异分析用counts来做 因为是把01A患者分组做差异分析所以读取01A
counts_01A <- read.table("counts01A.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
#因为是用免疫评分分组所以读取ESTIMATE_result
estimate <- read.table("ESTIMATE_result.txt", sep = "\t",row.names = 1,check.names = F,header = T)
#整理分组信息
x <- "StromalScore"
med <- as.numeric(median(estimate[,x]))
estimate <- as.data.frame(t(estimate))
identical(colnames(counts_01A),colnames(estimate))

conditions=data.frame(sample=colnames(counts_01A),
                      group=factor(ifelse(estimate[x,]>med,"high","low"),levels = c("low","high"))) %>% 
  column_to_rownames("sample")
#差异分析准备工作
dds <- DESeqDataSetFromMatrix(
  countData = counts_01A,
  colData = conditions,
  design = ~ group)

#开始差异分析
dds <- DESeq(dds)
#这句很重要
resultsNames(dds)
#提取结果
res <- results(dds)
save(res,file="DEG_StromalScore.Rda")

####热图绘制####
DEG <- as.data.frame(res)
#读取表达谱
exp <- read.table("tpms01A_log2.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
#添加上下调信息
logFC_cutoff <- 1
type1 = (DEG$padj < 0.05)&(DEG$log2FoldChange < -logFC_cutoff)
type2 = (DEG$padj < 0.05)&(DEG$log2FoldChange > logFC_cutoff)
DEG$change = ifelse(type1,"DOWN",ifelse(type2,"UP","NOT"))
table(DEG$change)

library(pheatmap)
#提取差异基因表达谱
a <- filter(DEG,change == 'UP')
b <- filter(DEG,change == 'DOWN')
c <- rbind(a,b)
d <- rownames(c)
exp_diff <- exp[d,]
#设置分组信息
annotation_col <- conditions
#对exp_diff 列的顺序进行处理
a <- filter(annotation_col,group == 'high')
b <- filter(annotation_col,group == 'low')
exp_diff_high <- exp_diff[,rownames(a)]
exp_diff_low <- exp_diff[,rownames(b)]
exp_diff <- cbind(exp_diff_high,exp_diff_low)
#开始画图
pheatmap(exp_diff,
         annotation_col=annotation_col,
         scale = "row",
         show_rownames = F,
         show_colnames =F,
         color = colorRampPalette(c("navy", "white", "red"))(50),
         cluster_cols =F,
         cluster_rows = T,
         fontsize = 10,
         fontsize_row=3,
         fontsize_col=3)
#保存图片 调整大小
dev.off()#关闭画板


####将两次差异分析的差异基因取交集####
setwd("../")
setwd("Immune_Stromal_DEG")
load("DEG_ImmuneScore.Rda")
library(tidyverse)
DEG <- as.data.frame(res)
# 添加上下调信息
logFC_cutoff = 1
type1 <- (DEG$padj < 0.05)&(DEG$log2FoldChange < -logFC_cutoff)
type2 <- (DEG$padj < 0.05)&(DEG$log2FoldChange > logFC_cutoff)
DEG$change <- ifelse(type1, "DOWN", ifelse(type2, "UP", "NOT"))
table(DEG$change)
# 提取上下调基因
a <- filter(DEG, change == "UP")
b <- filter(DEG, change == "DOWN")
write.csv(a, file="Immune_up.csv")
write.csv(b, file="Immune_down.csv")

load("DEG_StromalScore.Rda")
DEG <- as.data.frame(res)
# 添加上下调信息
logFC_cutoff = 1
type1 <- (DEG$padj < 0.05)&(DEG$log2FoldChange < -logFC_cutoff)
type2 <- (DEG$padj < 0.05)&(DEG$log2FoldChange > logFC_cutoff)
DEG$change <- ifelse(type1, "DOWN", ifelse(type2, "UP", "NOT"))
table(DEG$change)
# 提取上下调基因
a <- filter(DEG, change == "UP")
b <- filter(DEG, change == "DOWN")
write.csv(a, file="Stromal_up.csv")
write.csv(b, file="Stromal_down.csv")



####用交集差异基因做富集分析####
setwd("../")
setwd("FUJI_Immune_Stromal_DEG")
library(tidyverse)
library(BiocManager)
# 安装加载包
#BiocManager::install("clusterProfiler")
#BiocManager::install("org.Hs.eg.db")
library(org.Hs.eg.db)
#org.Hs.eg.db包主要注释人类基因:用于不同数据库ID间的转化
library(clusterProfiler)
# 导入DEG_final.txt
# 导入immune或stromal差异分析结果
load("DEG_ImmuneScore.Rda")
DEG <- as.data.frame(res)
DEG <- DEG[DEG_final$SYMBOL,]
DEG <- rownames_to_column(DEG, "SYMBOL")
# 基因ID转换
genelist <- bitr(DEG$SYMBOL, fromType="SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
DEG <- inner_join(DEG, genelist, by = "SYMBOL")
# 要用ENTREZID去做富集分析

####GO####
ego <- enrichGO(gene = DEG$ENTREZID,
                OrgDb = org.Hs.eg.db,
                ont = "all",
                pAdjustMethod = "BH",
                minGSSize = 1,
                pvalueCutoff =0.05, 
                qvalueCutoff =0.05,
                readable = TRUE)
ego_res <- ego@result
save(ego, ego_res, file="GO_DEG_final.Rda")
# GO富集分析从三个方面：BF，CC以及MF去看这些基因富集在哪些功能里面


####KEGG####
kk <- enrichKEGG(gene = DEG$ENTREZID,
                 organism = "hsa",
                 pvalueCutoff = 0.1,
                 qvalueCutoff =0.1)
kk_res <- kk@result
save(kk, kk_res, file="KEGG_DEG_final.Rda")
# KEGG富集分析总结这些基因都富集在哪个通路里面

# 网络图
#install.packages("ggnewscale")
library(ggnewscale)
List <- DEG$log2FoldChange
names(List) <- DEG$ENTREZID
head(List)
List <- sort(List, decreasing = T)
# 画图
# GO
cnetplot(ego, foldchange=List, circular = TRUE, colorEdge = TRUE)
# KEGG
cnetplot(kk, foldChange = List, circular = TRUE, colorEdge = TRUE)
dev.off()



####PPI####
#String：https://cn.string-db.org/cgi/input?sessionId=bXmYsv7CnUrH&input_page_active_form=multiple_identifiers


####COX####
# COX回归分析找出基因表达与肿瘤患者生存时间和生存状态的关系
setwd("COX")
#install.packages("survival")
#install.packages("forestplot")
library(survival)
library(forestplot)
library(tidyverse)
exp_surv_01A <- read.table("exp_surv_01A.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
# 我们只看差异基因和生存时间的关系
# 提取生存时间和差异基因
surv.expr <- cbind(exp_surv_01A[,1:2], exp_surv_01A[,DEG_final$SYMBOL])

# COX分析
Coxoutput <- NULL
for (i in 3:ncol(surv.expr)){
  g <- colnames(surv.expr)[i]
  cox <- coxph(Surv(OS.time, OS) ~ surv.expr[,i], data = surv.expr)
  coxSummary = summary(cox)
  Coxoutput <- rbind.data.frame(Coxoutput,
                                data.frame(gene = g,
                                           HR = as.numeric(coxSummary$coefficients[,"exp(coef)"])[1],
                                           z = as.numeric(coxSummary$coefficients[,"z"])[1],
                                           pvalue = as.numeric(coxSummary$coefficients[,"Pr(>|z|)"])[1],
                                           lower = as.numeric(coxSummary$conf.int[,3][1]),
                                           upper = as.numeric(coxSummary$conf.int[,4][1]),
                                           stringsAsFactors = F),
                                stringsAsFactors = F)
}
Coxoutput <- arrange(Coxoutput,pvalue)
# 筛选p value小于0.005的基因
gene_sig <- Coxoutput[Coxoutput$pvalue < 0.005,]
write.csv(gene_sig, file="gene_sig.csv")
topgene <- gene_sig

# 画森林图
tabletext <- cbind(c("Gene",topgene$gene),
                   c("HR",format(round(as.numeric(topgene$HR),3),nsmall = 3)),
                   c("lower 95%CI",format(round(as.numeric(topgene$lower),3),nsmall = 3)),
                   c("upper 95%CI",format(round(as.numeric(topgene$upper),3),nsmall = 3)),
                   c("pvalue",format(round(as.numeric(topgene$p),3),nsmall = 3)))
forestplot(labeltext=tabletext,
           mean=c(NA,as.numeric(topgene$HR)),
           lower=c(NA,as.numeric(topgene$lower)), 
           upper=c(NA,as.numeric(topgene$upper)),
           graph.pos=5,# 图在表中的列位置
           graphwidth = unit(.25,"npc"),# 图在表中的宽度比
           fn.ci_norm="fpDrawDiamondCI",# box类型选择钻石
           col=fpColors(box="#00A896", lines="#02C39A", zero = "black"),# box颜色
           
           boxsize=0.4,# box大小固定
           lwd.ci=1,
           ci.vertices.height = 0.1,ci.vertices=T,# 显示区间
           zero=1,# zero线横坐标
           lwd.zero=1.5,# zero线宽
           xticks = c(0.5,1,1.5),# 横坐标刻度根据需要可随意设置
           lwd.xaxis=2,
           xlab="Hazard ratios",
           txt_gp=fpTxtGp(label=gpar(cex=1.2),# 各种字体大小设置
                          ticks=gpar(cex=0.85),
                          xlab=gpar(cex=1),
                          title=gpar(cex=1.5)),
           hrzl_lines=list("1" = gpar(lwd=2, col="black"), # 在第一行上面画黑色实线
                           "2" = gpar(lwd=1.5, col="black"), # 在第一行标题行下画黑色实线
                           "53" = gpar(lwd=2, col="black")), # 在最后一行上画黑色实线
           lineheight = unit(.75,"cm"),# 固定行高
           colgap = unit(0.3,"cm"),
           mar=unit(rep(1.5, times = 4), "cm"),
           new_page = F
)
dev.off()


####BTK在肿瘤样本与正常样本中的表达差异####
# 柱状图
setwd("../")
setwd("BTK")
library(tidyverse)
tpms01A_log2 <- read.table("tpms01A_log2.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
tpms11A_log2 <- read.table("tpms11A_log2.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
gene <- "BTK"
a <- tpms01A_log2[gene,]
b <- tpms11A_log2[gene,]
# t转换
a <- a %>% t() %>% as.data.frame()
b <- b %>% t() %>% as.data.frame()
write.csv(a, file="BTK_01A.csv")
write.csv(b, file="BTK_11A.csv")

# 配对图绘制
tpms01A_log2 <- tpms01A_log2 %>% t() %>% as.data.frame()
tpms11A_log2 <- tpms11A_log2 %>% t() %>% as.data.frame()
rownames(tpms01A_log2) <- substring(rownames(tpms01A_log2),1,12)
rownames(tpms11A_log2) <- substring(rownames(tpms11A_log2),1,12)
a <- intersect(rownames(tpms01A_log2), rownames(tpms11A_log2))
tpms01A_log2 <- tpms01A_log2[a,]
tpms11A_log2 <- tpms11A_log2[a,]
peidui <- cbind(tpms11A_log2[,gene], tpms01A_log2[,gene])
peidui <- as.data.frame(peidui)
write.csv(peidui, file="peidui.csv")


####根据BTK高低组做生存分析####
setwd("../")
setwd("survival")
surv <- read.table("exp_surv_01A.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
# 将OS.time的单位从天变成年
surv$OS.time <- surv$OS.time/365

# median 中位数
# 找到BTK的中位数
# 根据BTK的中位数把样本分成high和low两组
surv$group <- ifelse(surv$BTK > median(surv$BTK), "High", "Low")
# 将High和Low从字符型转变为因子型数据
surv$group <- factor(surv$group, levels=c("Low", "High"))
class(surv$group)
table(surv$group)

# 生存分析
#install.packages("survival")
library(survival)
fitd <- survdiff(Surv(OS.time, OS) ~ group,
                 data = surv,
                 na.action = na.exclude)
pValue <- 1 - pchisq(fitd$chisq, length(fitd$n) - 1)

# 拟合生存曲线
fit <- survfit(Surv(OS.time, OS) ~ group, data=surv)
summary(fit)
p.lab <- paste0("P", ifelse(pValue < 0.001, "< 0.001", paste0(" = ", round(pValue, 3))))

# 画图
#install.packages("survminer")
library(survminer)
ggsurvplot(fit,
           data = surv,
           pval = p.lab,
           conf.int = TRUE, # 显示置信区间
           risk.table = TRUE, # 显示风险表
           risk.table.col = "strata",
           palette = "jco", # 配色采用jco
           legend.labs = c("Low", "High"), # 图例
           size = 1,
           xlim = c(0,20), # x轴长度
           break.time.by = 5, # x轴步长为5
           legend.title = "BTK",
           surv.median.line = "hv", # 限制垂直和水平的中位生存
           ylab = "Survival probability (%)", # 修改y轴标签
           xlab = "Time (Years)", # 修改x轴标签
           ncensor.plot = TRUE, # 显示删失图块
           ncensor.plot.height = 0.25,
           risk.table.y.text = FALSE)
dev.off()


####不同分期BTK的表达####
setwd("../")
setwd("BTK")
library(tidyverse)
clinical.expr01A <- read.table("clinical.expr01A.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
gene <- "BTK"
clinical_BTK <- cbind(clinical.expr01A[,1:6],clinical.expr01A[,gene])
write.csv(clinical_BTK, file="clinical_BTK.csv")


####BTK差异分析####
setwd("../")
setwd("BTK_DEG")
library(DESeq2)
library(tidyverse)
counts_01A <- read.table("counts01A.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
exp <- read.table("tpms01A_log2.txt", sep = "\t",row.names = 1,check.names = F,header = T)
identical(colnames(counts_01A), colnames(exp))
gene <- "BTK"
med <- median(as.numeric(exp[gene,]))

conditions=data.frame(sample=colnames(exp),
                      group=factor(ifelse(exp[gene,]>med,"high","low"),levels = c("low","high"))) %>% 
  column_to_rownames("sample")

dds <- DESeqDataSetFromMatrix(countData = counts_01A, colData = conditions, design= ~ group)
dds <- DESeq(dds)
resultsNames(dds)
res <- results(dds)
save(res, file="DEG_BTK.Rda")


####GSEA富集分析####
library(org.Hs.eg.db)
library(clusterProfiler)
DEG <- as.data.frame(res) %>% arrange(padj)
DEG <- DEG %>% rownames_to_column("Gene")
geneList = DEG[,3]
names(geneList) <- as.character(DEG[,"Gene"])
head(geneList)
geneList <- sort(geneList, decreasing = T)
head(geneList)

# GSEA基因集：https://zhuanlan.zhihu.com/p/504101161

msigdb_GMTs <- "msigdb_v7.0_GMTs"
msigdb <-"h.all.v7.0.symbols.gmt"
# 读取上面指定的gmt文件
kegmt <- read.gmt(file.path(msigdb_GMTs, msigdb))
# 设置种子
set.seed(1)
# GSEA分析
gsea <- GSEA(geneList, TERM2GENE = kegmt)
# 转换成数据框
gsea_result_df <- as.data.frame(gsea)
save(gsea, gsea_result_df, file="GSEA_BTK_h.all.rda")
# 绘图
# 安装enrichplot
library(enrichplot)
# 单个结果绘制
gseaplot2(gsea, 1, color="red", pvalue_table = T)
# 多个结果绘制
# A
gseaplot2(gsea, geneSetID = c(1,2,3,4,5,6,8,10), subplots = 1:3)
# B
gseaplot2(gsea, geneSetID = c(7,9,11,13,14,16,17), subplots = 1:3)
# 画前十条
gseaplot2(gsea, geneSetID = 1:10, subplots=1:3)
dev.off()

####换C7跑####
msigdb_GMTs <- "msigdb_v7.0_GMTs"
msigdb <-"c7.all.v7.0.symbols.gmt"
kegmt <- read.gmt(file.path(msigdb_GMTs, msigdb))
# 设置种子
set.seed(1)
# GSEA分析
gsea <- GSEA(geneList, TERM2GENE = kegmt)
# 转换成数据框
gsea_result_df <- as.data.frame(gsea)
save(gsea, gsea_result_df, file="GSEA_BTK_c7.rda")
# 绘图
library(enrichplot)
# 单个结果绘制
gseaplot2(gsea, 1, color="red", pvalue_table = T)
# 多个结果绘制
# C
gseaplot2(gsea, geneSetID = 1:7, subplots = 1:3)
# D
gseaplot2(gsea, 293, color="red",pvalue_table = T)
# 画前十条
gseaplot2(gsea, geneSetID = 1:10, subplots=1:3)
dev.off()


####cibersort####
setwd("../")
setwd("cibersort")
#install.packages("e1071")
#install.packages("parallel")
#BiocManager::install("preprocessCore", version="3.18")
library(e1071)
library(parallel)
library(preprocessCore)
library(tidyverse)
# cibersort计算每个样本里面22种免疫细胞所占比例
source("CIBERSORT.R")
sig_matrix <- "LM22.txt"
mixture_file <- "tpms01A_log2.txt"
res_cibersort <- CIBERSORT(sig_matrix, mixture_file, perm=100, QN=TRUE)
res_cibersort <- res_cibersort[,1:22]
ciber.res <- res_cibersort[, colSums(res_cibersort)>0]
ciber.res <- as.data.frame(ciber.res)
write.table(ciber.res, "ciber.res.txt",sep = "\t",row.names = T,col.names = NA,quote = F)


####cibersort彩虹图####
mycol <- ggplot2::alpha(rainbow(ncol(ciber.res)), 0.7)
par(bty="o", mgp = c(2.5,0.3,0), mar = c(2.1,4.1,2.1,10.1),tcl=-.25,las = 1,xpd = F)
barplot(as.matrix(t(ciber.res)),
        border = NA, # 柱子无边框
        names.arg = rep("",nrow(ciber.res)), # 无横坐标样本名
        yaxt = "n", # 先不绘制y轴
        ylab = "Relative percentage", # 修改y轴名称
        col = mycol) # 采用彩虹色板
axis(side = 2, at = c(0,0.2,0.4,0.6,0.8,1), # 补齐y轴添加百分号
     labels = c("0%","20%","40%","60%","80%","100%"))
legend(par("usr")[2]-20, # 
       par("usr")[4], 
       legend = colnames(ciber.res), 
       xpd = T,
       fill = mycol,
       cex = 0.6, 
       border = NA, 
       y.intersp = 1,
       x.intersp = 0.2,
       bty = "n")
dev.off()



# 分组比较图
a <- ciber.res
exp <- read.table("tpms01A_log2.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
med <- median(as.numeric(exp["BTK",]))
exp <- exp %>% t() %>% as.data.frame()
# mutate新增一列数据框
exp <- exp %>% mutate(group=factor(ifelse(exp$BTK>med, "high", "low"), levels=c("low", "high")))
class(exp$group)
identical(rownames(a), rownames(exp))
a$group <- exp$group
a <- rownames_to_column(a, "sample")
library(ggsci)
library(tidyr)
library(ggpubr)
b <- gather(a,key=CIBERSORT,value = Fraction,-c(group,sample))
ggboxplot(b, x = "CIBERSORT", y = "Fraction",
          fill = "group", palette = "lancet")+
  stat_compare_means(aes(group = group),
                     method = "wilcox.test",
                     label = "p.signif",
                     symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1),
                                      symbols = c("***", "**", "*", "ns")))+
  theme(text = element_text(size=10),
        axis.text.x = element_text(angle=45, hjust=1)) 
dev.off()


####相关性热图####
#install.packages("ggstatsplot")
#install.packages("ggcorrplot")
#install.packages("corrplot")
library(ggstatsplot)
library(ggcorrplot)
library(corrplot)
cor<-sapply(ciber.res,function(x,y) cor(x,y,method="spearman"),ciber.res)
rownames(cor) <- colnames(ciber.res)

ggcorrplot(cor, 
           hc.order = TRUE, #使用hc.order进行排序
           type = "upper", #图片位置
           outline.color = "white",#轮廓颜色
           lab = TRUE,#true为在图上添加相关系数
           ggtheme = ggplot2::theme_gray, #指ggplot2函数对象，默认值为thememinimal
           colors = c("#01468b", "white", "#ee0000"))
dev.off()


####基因与cibersort相关性散点图####
exp = read.table("tpms01A_log2.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
exp <- exp["BTK",]
ciber = read.table("ciber.res.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
ciber <- ciber %>% t() %>% as.data.frame()
rownames(ciber) <- gsub(" ",".",rownames(ciber))
identical(colnames(ciber), colnames(exp))
exp_ciber <- rbind(exp, ciber)
exp_ciber <- exp_ciber %>% t() %>% as.data.frame()
#install.packages("ggside")
library(ggstatsplot)
library(ggside)
ggscatterstats(data = exp_ciber, #要分析的数据
               y = BTK, #设置Y轴
               x = B.cells.naive,#设置X轴
               type = "nonparametric", 
               margins = "both",#是否显示 边缘，默认为true                                      
               xfill = "#01468b", #x轴边缘图形的颜色
               yfill = "#ee0000", #y轴边缘图形的颜色
               marginal.type = "densigram")#在图片坐标轴边缘添加图形类型
dev.off()
