###甲基化数据处理
##一、GSE153712 
library(meffil)
options(mc.cores=40) 
qc.objects1 <- meffil.qc(AIBL_samplesheet, cell.type.reference="blood gse167998", verbose=TRUE)
qc.objects2 <- meffil.qc(ADNI_samplesheet, cell.type.reference="blood gse167998", verbose=TRUE)
qc.objects3 <- meffil.qc(GSE208623_samplesheet, cell.type.reference="blood gse167998", verbose=TRUE)
qc.objects <- c(qc.objects1, qc.objects2,qc.objects3)
qc.summary <- meffil.qc.summary(qc.objects)
meffil.qc.report(qc.summary, output.file="qc/report.html")#qc报告
qc.objects <- meffil.remove.samples(qc.objects, qc.summary$bad.samples$sample.name)
norm.objects <- meffil.normalize.quantiles(qc.objects, number.pcs=10)
#normalization报告
norm.summary <- meffil.normalization.summary(norm.objects, pcs=10)
meffil.normalization.report(norm.summary, output.file="normalization/report.html")#
#保存β值
meffil.normalize.samples(
  norm.objects,
  just.beta=T,
  remove.poor.signal=T,
  cpglist.remove=qc.summary$bad.cpgs$name,
  gds.filename="total_beta.gds",
  verbose=T)
####过滤探针
load("~/zhutianshu/combined/filtered_sites.Rdata")
#常染色体探针
autosomal.sites <- meffil::meffil.get.autosomal.sites("epic")
#文献报道低质量探针
zhou <- read.csv("~/zhutianshu/combined/zhou_2016.csv")
std_names <- zhou$ID
filtered_sites <- setdiff(autosomal.sites, std_names)
#质控去掉的低质量探针
load("~/zhutianshu/combined/qc.summary.rda")
bad.sites <- qc.summary$bad.cpgs$name
filtered_sites <- setdiff(filtered_sites, bad.sites)
#读取保留探针数据
beta <- meffil.gds.methylation("~/zhutianshu/combined/total_beta.gds",sites=filtered_sites)
#按0.95过滤缺失值
#先过滤缺失多的样本
threshold <- 0.95
max_na <- floor((1 - threshold) * nrow(beta))
final_betas <- beta[, colSums(is.na(beta)) < max_na]
#过滤缺失多的探针
final_betas <- t(beta)
max_na <- floor((1 - threshold) * nrow(final_betas))
final_betas <- final_betas[, colSums(is.na(final_betas)) < max_na]
####缺失值填充
library(softImpute)
fita=softImpute(final_betas)
final_betas <- complete(final_betas, fita)
######找出AIBL中的数据
library(dplyr)
library(tibble)
AIBL_samplesheet <- read.csv("~/zhutianshu/GSE153712/yuanshishuju/AIBL_samplesheet.csv")
AIBL_beta<- final_betas[rownames(final_betas) %in% AIBL_samplesheet$Sample_Name, ]
AIBL_beta=as.data.frame(AIBL_beta)%>%
  rownames_to_column(var="Sample_Name")%>%
  merge(AIBL_samplesheet[,c(1,6)],.,by="Sample_Name")%>%
  column_to_rownames(var="Sample_Name")
save(AIBL_beta,file="~/zhutianshu/combined/BioM2/AIBL_beta.RData")
#找出对应位点注释文件
my_featureAnno <- read.csv("~/zhutianshu/combined/zhou_featureAnno.csv")
select_col <- intersect(colnames(AIBL_beta),my_featureAnno$ID)
AIBL_fea <- my_featureAnno[my_featureAnno$ID%in%select_col,]
save(AIBL_fea,file="~/zhutianshu/combined/BioM2/AIBL_fea.RData")


#二. GSE156984 
library(data.table)
GSE156984_STG <- fread("~/zhutianshu/GSE156984/GSE156984_STG_Matrix_processed.txt.gz",sep="\t")
GSE156984_STG <-GSE156984_STG%>%
  select(!matches("_Detection_Pval"))
GSE156984_STG1 <- GSE156984_STG%>%
  column_to_rownames(var="V1")%>%
  t()%>%
  as.data.frame()
#探针过滤
zhou <- read.csv("~/zhutianshu/combined/zhou_2016.csv")
std_names <- zhou$ID
select_col <- setdiff(colnames(GSE156984_STG1),std_names)
GSE156984_STG <- GSE156984_STG1[,select_col]
save(GSE156984_STG,file="~/zhutianshu/GSE156984/GSE156984_STG.RData")
#加label
GSE156984_anno=read.csv("~/zhutianshu/GSE156984/GSE156984_anno.csv")
GSE156984_STG <- GSE156984_STG%>%
  rownames_to_column(var='sample')%>%
  merge(GSE156984_anno[,c(4,6)],.,by='sample')%>%
  column_to_rownames(var='sample')
save(GSE156984_STG,file="~/zhutianshu/combined/BioM2/GSE156984_STG.RData")
##########选择适合的注释文件
my_featureAnno <- read.csv("~/zhutianshu/combined/zhou_featureAnno.csv")
std_names <- my_featureAnno$ID
select_col <- intersect(colnames(GSE156984_STG),std_names)
GSE156984_STG_fea <- my_featureAnno[my_featureAnno$ID%in%select_col,]
save(GSE156984_STG_fea,file="~/zhutianshu/combined/BioM2/GSE156984_STG_fea.RData")


#三. GSE134379/GSE66531
##samplesheet
GSE66531_anno <- read.csv("GSE66351_series_matrix.csv",header = F)
GSE66531_anno <- t(GSE66531_anno)%>%
  as.data.frame()
GSE66531_anno1 <- GSE66531_anno[2:nrow(GSE66531_anno),]
colnames(GSE66531_anno1) <- c("Sample_Name","region","cell type","label","AGE","Sex","Slide","Array")
GSE66531_anno2 <- GSE66531_anno1%>%
  dplyr::mutate(across(3,~sapply(strsplit(.x, " "), `[`, 3)))%>%
  dplyr::mutate(across(c(4:8),~sapply(strsplit(.x, " "), `[`, 2)))%>%
  dplyr::mutate(label=ifelse(label=="AD",1,0))%>%
  filter(!(region=="Temporal Cortex"&`cell type`=='bulk'))
write.csv(GSE66531_anno2,"~/zhutianshu/GSE66531/GSE66531_samplesheet.csv",row.names = F)
GSE66531_samplesheet <- meffil::meffil.read.samplesheet("~/zhutianshu/GSE66531/")
write.csv(GSE66531_samplesheet,"~/zhutianshu/GSE66531/GSE66531_samplesheet.csv")
GSE66531_samplesheet <- read.csv("~/zhutianshu/GSE66531/GSE66531_samplesheet.csv")
GSE134379_anno <- read.csv("GSE134379_series_matrix.csv",header = F)
GSE134379_anno <- t(GSE134379_anno)%>%
  as.data.frame()
library(stringr)
GSE134379_anno1 <- GSE134379_anno[2:nrow(GSE134379_anno),]
colnames(GSE134379_anno1) <- c("Sample_Name","region","Sex","AGE","label")
GSE134379_anno2 <- GSE134379_anno1%>%
  dplyr::mutate(across(2,~sapply(strsplit(.x, " "), `[`, 3)))%>%
  dplyr::mutate(across(c(3:5),~sapply(strsplit(.x, " "), `[`, 2)))%>%
  dplyr::mutate(label=ifelse(label=="AD",1,0))
shell-'ls > name.txt'
name_file <- read.table("~/zhutianshu/GSE134379/name.txt")
name_file <- head(name_file, -1)
name_file2 <- name_file%>%
  dplyr::mutate(Sample_Name=sapply(strsplit(V1, "_"), `[`, 1))%>%
  dplyr::mutate(Array=sapply(strsplit(V1, "_"), `[`, 3))%>%
  dplyr::mutate(Slide=sapply(strsplit(V1, "_"), `[`, 2))%>%
  distinct(Sample_Name,.keep_all = T)
GSE134379_anno3 <- merge(name_file2,GSE134379_anno2,by="Sample_Name")%>%
  dplyr::select(-2)
write.csv(GSE134379_anno3,"~/zhutianshu/GSE134379/GSE134379_samplesheet.csv",row.names = F)
GSE134379_samplesheet <- meffil::meffil.read.samplesheet("~/zhutianshu/GSE134379/")
write.csv(GSE134379_samplesheet,"~/zhutianshu/GSE134379/GSE134379_samplesheet.csv")
#meffil处理原始数据
library(meffil)
options(mc.cores=40)
GSE134379_samplesheet <-read.csv("~/zhutianshu/GSE134379/GSE134379_samplesheet.csv")
GSE66531_samplesheet <-read.csv("~/zhutianshu/GSE66531/GSE66531_samplesheet.csv")
qc.objects1 <- meffil.qc(GSE134379_samplesheet, cell.type.reference="guintivano dlpfc", verbose=TRUE)
qc.objects2 <- meffil.qc(GSE66531_samplesheet, cell.type.reference="guintivano dlpfc", verbose=TRUE)
qc.objects <- c(qc.objects1, qc.objects2)
save(qc.objects,file="brain_qc.objects.rda")
qc.summary <- meffil.qc.summary(qc.objects)
save(qc.summary,file="brain_qc.summary.rda")
#meffil.qc.report(qc.summary, output.file="qc/report.html")#qc报告
qc.objects <- meffil.remove.samples(qc.objects, qc.summary$bad.samples$sample.name)
norm.objects <- meffil.normalize.quantiles(qc.objects, number.pcs=10)
save(norm.objects,file="brain_norm.objects.rda")
#normalization报告
#norm.summary <- meffil.normalization.summary(norm.objects, pcs=10)
#meffil.normalization.report(norm.summary, output.file="normalization/report.html")#
#保存β值
meffil.normalize.samples(
  norm.objects,
  just.beta=T,
  remove.poor.signal=T,
  cpglist.remove=qc.summary$bad.cpgs$name,
  gds.filename="brain_total_beta.gds",
  verbose=T)
#读取gds数据
beta <- meffil::meffil.gds.methylation("~/zhutianshu/combined/brain_total_beta.gds", sites=filtered_sites)
#先过滤缺失多的样本
threshold <- 0.95
max_na <- floor((1 - threshold) * nrow(beta))
final_betas <- beta[, colSums(is.na(beta)) < max_na]
#过滤缺失多的探针
final_betas <- t(beta)
max_na <- floor((1 - threshold) * nrow(final_betas))
final_betas <- final_betas[, colSums(is.na(final_betas)) < max_na]
####缺失值填充
library(softImpute)
fita=softImpute(final_betas)
final_betas <- complete(final_betas, fita)
save(final_betas,file="~/zhutianshu/combined/BioM2/brain_beta.RData")
######找出GSE134379中的数据
library(dplyr)
library(tibble)
GSE134379_samplesheet <- read.csv("~/zhutianshu/GSE134379/GSE134379_samplesheet.csv")
###middle temporal gyrus (MTG)
name <- GSE134379_samplesheet[GSE134379_samplesheet$region=='MTG',]$Sample_Name
MTG_beta<- final_betas[rownames(final_betas) %in% name, ]
MTG_beta=as.data.frame(MTG_beta)%>%
  rownames_to_column(var="Sample_Name")%>%
  merge(GSE134379_samplesheet[,c(1,7)],.,by="Sample_Name")%>%
  column_to_rownames(var="Sample_Name")
save(MTG_beta,file="~/zhutianshu/combined/BioM2/MTG_beta.RData")
####找出GSE66531中的数据
GSE66531_samplesheet <- read.csv("~/zhutianshu/GSE66531/GSE66531_samplesheet.csv")
###froCortex
name <- GSE66531_samplesheet[GSE66531_samplesheet$region=='Frontal Cortex',]$Sample_Name
froCortex_beta<- final_betas[rownames(final_betas) %in% name, ]
froCortex_beta=as.data.frame(froCortex_beta)%>%
  rownames_to_column(var="Sample_Name")%>%
  merge(GSE66531_samplesheet[,c(1,4)],.,by="Sample_Name")%>%
  column_to_rownames(var="Sample_Name")
save(froCortex_beta,file="~/zhutianshu/combined/BioM2/froCortex_beta.RData")
###neuron
name <- GSE66531_samplesheet[GSE66531_samplesheet$cell.type=='Neuron',]$Sample_Name
Neuron_beta<- final_betas[rownames(final_betas) %in% name, ]
Neuron_beta=as.data.frame(Neuron_beta)%>%
  rownames_to_column(var="Sample_Name")%>%
  merge(GSE66531_samplesheet[,c(1,4)],.,by="Sample_Name")%>%
  column_to_rownames(var="Sample_Name")
save(Neuron_beta,file="~/zhutianshu/combined/BioM2/Neuron_beta.RData")
###glia
name <- GSE66531_samplesheet[GSE66531_samplesheet$cell.type=='Glia',]$Sample_Name
Glia_beta<- final_betas[rownames(final_betas) %in% name, ]
Glia_beta=as.data.frame(Glia_beta)%>%
  rownames_to_column(var="Sample_Name")%>%
  merge(GSE66531_samplesheet[,c(1,4)],.,by="Sample_Name")%>%
  column_to_rownames(var="Sample_Name")
save(Glia_beta,file="~/zhutianshu/combined/BioM2/Glia_beta.RData")


#二、选出BP通路
xx <- as.list(GOTERM)
go_terms <- data.frame(
  GO_ID = names(xx),
  Description = sapply(xx, function(term) Term(term)),
  Ontology= sapply(xx, function(Ontology) Ontology(Ontology)),
  stringsAsFactors = FALSE
)
go_terms <- go_terms%>%
  filter(Ontology=="BP")
#通路基因注释筛选
go_BP <- as.list(org.Hs.egGO2ALLEGS)
go_BP <- go_BP[names(go_BP) %in% go_terms$GO_ID]
save(go_BP,file="~/zhutianshu/combined/BioM2/go_BP.RData")


