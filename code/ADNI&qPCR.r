load("CNAD_ROC.RData")
a=CNAD_ROC[["Prediction"]]
CNAD_PMS <- data.frame(sample = character(), prediction = numeric())
for (i in 1:length(a)) {
  df <- a[[i]]
  CNAD_PMS <- rbind(CNAD_PMS, df)  # 合并数据框
}
CNAD_PMS <- CNAD_PMS%>%
  dplyr::rename(Sample_Name=sample)%>%
  merge(ADNI_samplesheet,.,by="Sample_Name")%>%
  dplyr::select(-c(2,4,5))
#添加表型信息
PMS_clin <- adnimerge[,c(1,7,16,18,20:22,27,55,62)]
PMS_clin[, 5:7] <- apply(PMS_clin[, 5:7], 2, function(x) suppressWarnings(as.numeric(as.character(x))))
CNAD_PMS_clin_bl <- CNAD_PMS%>%
  dplyr::rename(RID=Sample_Name)%>%
  merge(.,PMS_clin,by="RID")%>%
  mutate(TimeInterval = abs(as.numeric(difftime(Edate, EXAMDATE, units = "days"))))%>%
  group_by(RID)%>%
  arrange(TimeInterval)%>%
  slice(1)%>%
  ungroup()%>%
  dplyr::select(-c(2,18))%>%
  dplyr::rename(Edate=EXAMDATE)
##########################################################################
##########################################################################################
#纵向研究的函数
library(arm)
library(lme4)
library(lmerTest)

lmer_data <- function(data, variable_list) {
  result_df <- data.frame()
  for (variable in variable_list) {
    model <- lmer(as.formula(paste0("scale(", variable, ") ~ PMS*scale(month) + AGE + Sex + PTEDUCAT + APOE4 + (1+month|RID)")), data = data)
    params <- summary(model)$coefficients[c(2, 3, 8), c("Estimate", "Std. Error", "Pr(>|t|)")]
    row_names_vector <- rownames(params)
    params <- cbind(type = row_names_vector, params)
    rownames(params) <- NULL
    params <- cbind(Variable = variable, params)
    result_df <- rbind(result_df, params)
  }
  return(result_df)
}  


#纵向ADSP改变
ADSP=read.csv("~/ADNI/ADNI_BIO/plasma_mix/cognition/ADSP_PHC_COGN_10_05_22_20Dec2023.csv")
ADSP=ADSP[, c(1,6,14,15,16,17)]
ADSP$EXAMDATE <- as.Date(ADSP$EXAMDATE)
ADSP_PMS_clin <- CNAD_PMS_clin_bl[,c(1:8)]%>%
  merge(.,ADSP,by="RID")%>%
  mutate(month = as.numeric(difftime(EXAMDATE,Edate , units = "days")/30))%>%
  filter(month>=0)%>%
  mutate(PMS=ifelse(prediction>0.5,1,0))%>%
  group_by(RID) %>%
  filter(n() > 1)
ADSP_result <- data.frame()
variables_to_replace <- c("PHC_MEM", "PHC_EXF", "PHC_LAN")
ADSP_result <- rbind(ADSP_result,lmer_data(CNAD_PMS_clin, variables_to_replace))

#纵向FDG改变
FDG_PMS_clin <- CNAD_PMS_clin_bl[,c(1:8)]%>%
  merge(.,PMS_clin[,c(1,2,3)],by="RID")%>%
  mutate(month = as.numeric(difftime(EXAMDATE,Edate , units = "days")/30))%>%
  filter(month>=0)%>%
  mutate(PMS=ifelse(prediction>0.5,1,0))%>%
  filter(!is.na(FDG))%>%
  group_by(RID) %>%
  filter(n() > 1)
ADSP_result <- rbind(ADSP_result,lmer_data(FDG_PMS_clin, "FDG"))
#纵向av45
ucberkeleyav45 <- data.frame(ucberkeleyav45)
av45_z <- ucberkeleyav45[, c(2, 4,19)]
av45_z <- av45_z%>%
  dplyr::rename(av45=SUMMARYSUVR_COMPOSITE_REFNORM)
av45_PMS_clin <- CNAD_PMS_clin_bl[,c(1:8)]%>%
  merge(.,av45_z,by="RID")%>%
  mutate(month = as.numeric(difftime(EXAMDATE,Edate , units = "days")/30))%>%
  filter(month>=0)%>%
  mutate(PMS=ifelse(prediction>0.5,1,0))%>%
  filter(!is.na(av45))%>%
  group_by(RID) %>%
  filter(n() > 1)
ADSP_result <- rbind(ADSP_result,lmer_data(av45_PMS_clin, 'av45'))
#纵向PTAU
PTAU_PMS_clin <- CNAD_PMS_clin_bl[,c(1:8)]%>%
  merge(.,PMS_clin[,c(1,2,7)],by="RID")%>%
  mutate(month = as.numeric(difftime(EXAMDATE,Edate , units = "days")/30))%>%
  filter(month>=0)%>%
  mutate(PMS=ifelse(prediction>0.5,1,0))%>%
  filter(!is.na(PTAU))%>%
  group_by(RID) %>%
  filter(n() > 1)
ADSP_result <- rbind(ADSP_result,lmer_data(PTAU_PMS_clin, "PTAU"))

#纵向MMSE
MMSE_PMS_clin <- CNAD_PMS_clin_bl[,c(1:8)]%>%
  merge(.,PMS_clin[,c(1,2,8)],by="RID")%>%
  mutate(month = as.numeric(difftime(EXAMDATE,Edate , units = "days")/30))%>%
  filter(month>=0)%>%
  mutate(PMS=ifelse(prediction>0.5,1,0))%>%
  filter(!is.na(MMSE))%>%
  group_by(RID) %>%
  filter(n() > 1)
ADSP_result <- rbind(ADSP_result,lmer_data(MMSE_PMS_clin, "MMSE"))
#纵向的mPACCdigit
mPACCdigit_PMS_clin <- CNAD_PMS_clin_bl[,c(1:8)]%>%
  merge(.,PMS_clin[,c(1,2,10)],by="RID")%>%
  mutate(month = as.numeric(difftime(EXAMDATE,Edate , units = "days")/30))%>%
  filter(month>=0)%>%
  mutate(PMS=ifelse(prediction>0.5,1,0))%>%
  filter(!is.na(mPACCdigit))%>%
  group_by(RID) %>%
  filter(n() > 1)
ADSP_result <- rbind(ADSP_result,lmer_data(mPACCdigit_PMS_clin, "mPACCdigit"))
write.csv(ADSP_result,"PMS_lmer.csv",row.names = F)
#斜率绘图
library(sloper)
library(ggplot2)
library(colorspace)
myplot_slopes <- function(data = NULL,
                          r_id = 'r_id', # column identifying raters
                          response = "response", # column of ratings or scores
                          contingency = "contingency", # column of rating contingencies
                          xname = NULL,
                          yname = NULL,
                          groupfactor = NULL,
                          gflevels = NULL,
                          linear = T
) {
  # do this super annoying thing to get the actual variables
  df <- data
  r_id <- df[[r_id]]
  response <- df[[response]]
  contingency <- df[[contingency]]
  groupfactor <- factor(df[[groupfactor]])
  if(!is.null(gflevels)){
    invisible(capture.output(levels(groupfactor) <- dput(gflevels)))
  }
  p <- ggplot(df, aes(x = contingency, y = response)) +
    geom_smooth(aes(group = r_id, color = groupfactor,
                    linetype = groupfactor),
                method = "lm",
                se = F, linewidth = .2,alpha = 0.5) +
    scale_linetype_manual(values = c("0" = "dashed", "1" = "solid")) +
    scale_color_manual(values = c("0" = "steelblue", "1" = "darkorange")) +
    theme_classic() +
    ylab(yname) +
    xlab(xname) +
    theme(legend.position = "bottom", legend.title=element_blank())
  return(p)
}
myplot_slopes(PTAU_PMS_clin, r_id = "RID",response = "PTAU", contingency = "month",
              yname = "CSF PTAU", xname = "month",gflevels=c(0,1),
              groupfactor = "PMS")
#基线统计
PMS_bl_ADSP <- CNAD_PMS_clin_bl[,-10]%>%
  merge(.,ADSP,by="RID")%>%
  mutate(month = as.numeric(difftime(EXAMDATE,Edate , units = "days")/30))%>%
  filter(month>=0)%>%
  mutate(PMS=ifelse(prediction>0.5,1,0))%>%
  group_by(RID) %>%
  filter(n() > 1)%>%
  arrange(month)%>%
  slice(1)
PMS_bl_av45 <- PMS_bl[,-16]%>%
  merge(.,av45_z,by="RID")%>%
  mutate(month = as.numeric(difftime(EXAMDATE,Edate , units = "days")/30))%>%
  filter(month>=0)%>%
  group_by(RID) %>%
  filter(n() > 1)%>%
  arrange(month)%>%
  slice(1)
PMS_bl <-  PMS_bl_ADSP%>%
  merge(.,PMS_bl_av45[,c(1,24)],by="RID",all=T)


PMS_s <- CNAD_PMS_clin_bl%>%
  merge(.,ADSP,by="RID")%>%
  mutate(month = as.numeric(difftime(EXAMDATE,Edate , units = "days")/30))%>%
  filter(month>=0)%>%
  mutate(PMS=ifelse(prediction>0.5,1,0))%>%
  group_by(RID) %>%
  filter(n() > 1)%>%
  arrange(desc(month))%>%
  slice(1)%>%
  dplyr::rename(month_s=month)%>%
  dplyr::select(c(1,22))
PMS_bl_st <- merge(PMS_bl,PMS_s,by='RID')%>%
  dplyr::select(-c(-1,8,10:11,14,16,20,21))%>%
  mutate(label=as.factor(label),PMS=as.factor(PMS),APOE4=as.factor(APOE4))
Table_S4 <-  merge(PMS_bl,PMS_s,by='RID')%>%
  dplyr::select(-c(8,10:11,14,16,20,21))
write.csv(Table_S4,"Table_S4.csv",row.names = F)
library(gtsummary)
library(gt)
#分组统计
table1 <-
  tbl_summary(
    PMS_bl_st,
    include = everything(),
    by = label, # split table by group
    statistic = list(
      all_continuous() ~ "{mean} ({sd})",
      all_categorical() ~ "{n} / {N} ({p}%)"
    ),
    missing = "no" # don't list missing data separately
  ) %>%
  add_n() %>% # add column with total number of non-missing observations
  add_p(all_continuous()~'t.test') %>% # test for a difference between groups
  modify_header(label = "**Variable**") %>% # update the column header
  bold_labels()#%>%
#as_gt() %>%
#gt::tab_source_note(gt::md("*This data is simulated*"))
#保存
as_gt(table1)
as_gt(table1) %>%
  gt::gtsave(filename = "table2.docx")
#整体统计
table2 <-
  tbl_summary(
    PMS_bl_st,
    include = everything(),
    statistic = list(
      all_continuous() ~ "{mean} ({sd})",
      all_categorical() ~ "{n} / {N} ({p}%)"
    ),
    missing = "no" # don't list missing data separately
  ) %>%
  add_n()
#保存
as_gt(table2)%>%
  gt::gtsave(filename = "table2.docx")       

###############################################################
###################qPCR计算
####函数定义
#随机抽三个数据构成新的表
library(qPCRtools)
all_combinations <- expand.grid(lapply(test, function(x) combn(x, 3, simplify = FALSE)))
#短数据--长数据--两个table
create_tables <- function(df) {
  long_data <- df %>%
    pivot_longer(
      cols = -BioRep,  # 排除“样本名”列
      names_to = c("Group", "Gene"),  # 将列名拆分为“group”和“基因”
      names_sep = "\\.",  # 以“.”为分隔符
      values_to = "Cq"  # 值列命名为“Cq”
    ) %>%
    select(Group, BioRep, Gene, Cq) %>%
    mutate(Position = paste0("A", row_number()))
  cq.table <- long_data %>%
    select(Position, Cq)
  design.table <- long_data %>%
    select(Position, Group, BioRep, Gene)
  return(list(cq.table = cq.table, design.table = design.table))
}
#计算显著性
run_CalExp2ddCt <- function(tables) {
  cq.table <- tables$cq.table
  design.table <- tables$design.table
  res <- CalExp2ddCt(
    cq.table,
    design.table,
    ref.gene = "G", ## 内参
    ref.group = "Ctr", ## 对照
    stat.method = "t.test", ## 统计方法
    remove.outlier = FALSE,
    fig.type = "bar",
    fig.ncol = NULL
  )
  return(res)
}