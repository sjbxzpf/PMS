###甲基化数据处理
##一、GSE153712 
library(meffil)
options(mc.cores=40) 
qc.objects <- meffil.qc(AIBL_samplesheet, cell.type.reference="blood gse167998", verbose=TRUE)
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
save(AIBL_beta,file="~/zhutianshu/combined/BioM2/AIBL_beta.RData")
AIBL_beta<- final_betas[rownames(final_betas) %in% AIBL_samplesheet$Sample_Name, ]
AIBL_beta=as.data.frame(AIBL_beta)%>%
  rownames_to_column(var="Sample_Name")%>%
  merge(AIBL_samplesheet[,c(1,6)],.,by="Sample_Name")%>%
  column_to_rownames(var="Sample_Name")
#找出对应位点注释文件
my_featureAnno <- read.csv("~/zhutianshu/combined/zhou_featureAnno.csv")
select_col <- intersect(colnames(AIBL_beta),my_featureAnno$ID)
AIBL_fea <- my_featureAnno[my_featureAnno$ID%in%select_col,]
save(AIBL_fea,file="~/zhutianshu/combined/BioM2/AIBL_fea.RData")



#二. GSE66531
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
#meffil处理原始数据
library(meffil)
options(mc.cores=40)
GSE66531_samplesheet <-read.csv("~/zhutianshu/GSE66531/GSE66531_samplesheet.csv")
qc.objects <- meffil.qc(GSE66531_samplesheet, cell.type.reference="guintivano dlpfc", verbose=TRUE)
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
#三. ADNI
samplesheet <- read.csv("~/zhutianshu/ADNI/ADNI_meth/ADNI_Methylation_SampleAnnotation_20170530.csv")
library(dplyr)
library(ADNIMERGE)
samplesheet_bl <- samplesheet%>%
  group_by(RID)%>%
  arrange(Edate)%>%
  slice(1)%>%
  ungroup()%>%
  dplyr::select(-c(3,5,6))
#添加非纵向临床信息
clin <- adnimerge[,c(1,7,9:11,15)]
samplesheet_AGE <- clin%>%
  group_by(RID)%>%
  arrange(EXAMDATE)%>%
  slice(1)%>%
  ungroup()%>%
  merge(samplesheet_bl,.,by="RID")%>%
  mutate(TimeInterval = floor(as.numeric(difftime(Edate, EXAMDATE, units = "days"))/365)) %>%
  mutate(TimeInterval=ifelse(TimeInterval<0,0,TimeInterval))%>%
  mutate(AGE=AGE+TimeInterval)%>%
  dplyr::select(-c(6,11))
#添加诊断
samplesheet_clin <- adnimerge[,c(1,7,61)]%>%
  filter(!is.na(DX))%>%
  merge(samplesheet_AGE,.,by="RID")%>%
  mutate(tm=abs(as.numeric(difftime(Edate,EXAMDATE, units = "days"))))%>%
  group_by(RID)%>%
  arrange(tm)%>%
  slice(1)%>%
  dplyr::select(-c(10,12))%>%
  mutate(DX=as.numeric(factor(DX,
                              levels=c("CN","MCI","Dementia"),
                              labels=c(1,2,3))))%>%
  dplyr::rename(Sex=PTGENDER)

#修改为处理甲基化所需格式
samplesheet_merffil <- samplesheet_clin%>%
  mutate(Sex=ifelse(Sex=="Male","M","F"))%>%
  dplyr::rename(Sample_Name=barcodes)
samplesheet_path <- meffil.create.samplesheet("~/zhutianshu/ADNI/ADNI_meth/ADNI_iDAT_files/")
ADNI_samplesheet <- samplesheet_path[,c(1,6)]%>%
  merge(.,samplesheet_merffil,by="Sample_Name")%>%
  dplyr::select(-1)%>%
  dplyr::rename(Sample_Name=RID,label=DX)
write.csv(ADNI_samplesheet,"~/zhutianshu/ADNI/ADNI_meth/ADNI_iDAT_files/ADNI_samplesheet.csv",row.names = F)
#############################################################################################################
#处理甲基化原始数据
#/opt/R/4.3.1/lib/R/bin/R
library(meffil)
options(mc.cores=40)
ADNI_samplesheet <- read.csv("ADNI_samplesheet.csv")
qc.objects <- meffil.qc(ADNI_samplesheet, cell.type.reference="blood gse167998", verbose=TRUE)
qc.summary <- meffil.qc.summary(qc.objects)
meffil.qc.report(qc.summary, output.file="qc/report.html")#qc报告
qc.objects <- meffil.remove.samples(qc.objects, qc.summary$bad.samples$sample.name)
norm.objects <- meffil.normalize.quantiles(qc.objects, number.pcs=10)
meffil.normalize.samples(
  norm.objects,
  just.beta=T,
  remove.poor.signal=T,
  cpglist.remove=qc.summary$bad.cpgs$name,
  gds.filename="ADNI_beta.gds",
  verbose=T)
#过滤探针
autosomal.sites <- meffil::meffil.get.autosomal.sites("epic")
zhou <- read.csv("~/zhutianshu/combined/zhou_2016.csv")
std_names <- zhou$ID
filtered_sites <- setdiff(autosomal.sites, std_names)
bad.sites <- qc.summary$bad.cpgs$name
filtered_sites <- setdiff(filtered_sites, bad.sites)
beta <- meffil.gds.methylation("ADNI_beta.gds",sites=filtered_sites)
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
save(final_betas,file="final_betas.rds")
library(dplyr)
library(tibble)
MCI_beta<- final_betas[rownames(final_betas) %in% ADNI_samplesheet[ADNI_samplesheet$label==2,]$Sample_Name, ]
MCI_beta=as.data.frame(MCI_beta)%>%
  rownames_to_column(var="Sample_Name")%>%
  merge(ADNI_samplesheet[,c(2,10)],.,by="Sample_Name")%>%
  column_to_rownames(var="Sample_Name")
save(MCI_beta,file="MCI_beta.rds")
CNAD_beta<- final_betas[rownames(final_betas) %in% ADNI_samplesheet[ADNI_samplesheet$label!=2,]$Sample_Name, ]
CNAD_beta=as.data.frame(CNAD_beta)%>%
  rownames_to_column(var="Sample_Name")%>%
  merge(ADNI_samplesheet[,c(2,10)],.,by="Sample_Name")%>%
  column_to_rownames(var="Sample_Name")%>%
  mutate(label=ifelse(label==1,1,2))
save(CNAD_beta,file="CNAD_beta.rds")
#nohup /opt/R/4.3.1/lib/R/bin/Rscript ADNI_meffi.R > output.log 2>&1 &

#找出对应位点注释文件
my_featureAnno <- read.csv("~/zhutianshu/combined/zhou_featureAnno.csv")
select_col <- intersect(colnames(final_betas),my_featureAnno$ID)
ADNI_fea <- my_featureAnno[my_featureAnno$ID%in%select_col,]
save(ADNI_fea,file="~/zhutianshu/ADNI/ADNI_PMS/ADNI_fea.RData")
###################################################################################
#四、选出BP通路
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


