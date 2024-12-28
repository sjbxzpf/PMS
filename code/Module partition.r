#########
#定义相似公式  
jaccard_similarity <- function(set1, set2) {
  intersection <- length(intersect(set1, set2))
  union <- length(union(set1, set2))
  return(intersection / union)
}
#min_modules函数
min_modules <- function(data,num) {
  Cor_list <- list()
  for (x in 1:length(num)){
    set.seed(1)
    #数据处理
    a <- as.data.frame(t(data$PathwaysMatrix))
    a2 <- a[-1,]
    #模块划分
    result <- predictModules(a2, n.comp = num[x])
    S.vp = partition(result$S, method = "ann") 
    S.vp.split = splitMatrix(S.vp)
    vp.counts = apply(S.vp.split, 2, sum)
    S.vp.split <- S.vp.split[,vp.counts>5]
    colnames(S.vp.split) <- seq_len(ncol(S.vp.split))
    ModuleResult <- data.frame(ID = character(), cluster = character(), stringsAsFactors = FALSE)
    for (i in 1:nrow(S.vp.split)) {
      for (j in 1:ncol(S.vp.split)) {
        if (S.vp.split[i, j] == 1) {
          ModuleResult <- rbind(ModuleResult, data.frame(ID = rownames(S.vp.split)[i], 
                                                         cluster = colnames(S.vp.split)[j]))
        }
      }
    }
    #计算模块之间的相关性并筛选
    clusters <- ModuleResult %>%
      group_by(cluster) %>%
      summarize(IDs = list(ID))
    #设定阈值
    threshold <- 0.15
    # 过滤相似性过高的集合
    to_remove <- vector()
    for (i in 1:(nrow(clusters) - 1)) {
      for (j in (i + 1):nrow(clusters)) {
        if (j %in% to_remove) next
        sim <- jaccard_similarity(clusters$IDs[[i]], clusters$IDs[[j]])
        if (sim > threshold) {
          to_remove <- c(to_remove, j)
        }
      }
    }
    filtered_clusters <- clusters[-to_remove,]
    if (nrow(filtered_clusters) <= 8) {
      message("Filtered clusters are less than or equal to 8, moving to the next iteration.")
      Cor_list[[x]] <- data.frame(module = NA, cor = 0) 
      next
    }
    #过滤模块
    ModuleResult <- ModuleResult[ModuleResult$cluster%in%filtered_clusters$cluster,]
    #计算模块水平
    datax=data$PathwaysMatrix[,ModuleResult$ID]
    datax=moduleEigengenes(datax,ModuleResult$cluster)$eigengenes
    #计算模块相关性  
    datax$label=data$PathwaysMatrix[,'label']
    datax=datax[,c(ncol(datax),(1:(ncol(datax)-1)))]
    Cor=data.frame(module=rownames(as.data.frame(stats::cor(datax)))[-1],cor=as.data.frame(stats::cor(datax))[-1,1])
    Cor_list[[x]] <- Cor}
    Cor_list <- Filter(Negate(is.null), Cor_list)
    max_Cor_df_index <- which.max(sapply(Cor_list, function(m) max(m$cor, na.rm = TRUE)))
    num.comp <- num[max_Cor_df_index]
    return(list(Cor_list = Cor_list, num.comp = num.comp))
}
#######brain_Modules函数########
#找到需要得数据框并再次处理
brain_Modules <- function(data,num.comp) {
  set.seed(1)
  a <- as.data.frame(t(data$PathwaysMatrix))
  a2 <- a[-1,]
  #模块划分
  result <- predictModules(a2, n.comp = num.comp)
  S.vp = partition(result$S, method = "ann") 
  S.vp.split = splitMatrix(S.vp)
  vp.counts = apply(S.vp.split, 2, sum)
  S.vp.split <- S.vp.split[,vp.counts>5]
  colnames(S.vp.split) <- seq_len(ncol(S.vp.split))
  ModuleResult <- data.frame(ID = character(), cluster = character(), stringsAsFactors = FALSE)
  for (i in 1:nrow(S.vp.split)) {
    for (j in 1:ncol(S.vp.split)) {
      if (S.vp.split[i, j] == 1) {
        ModuleResult <- rbind(ModuleResult, data.frame(ID = rownames(S.vp.split)[i], 
                                                       cluster = colnames(S.vp.split)[j]))
      }
    }
  }
  #计算模块之间的相关性并筛选
  clusters <- ModuleResult %>%
    group_by(cluster) %>%
    summarize(IDs = list(ID))
  #设定阈值
  threshold <- 0.15
  # 过滤相似性过高的集合
  to_remove <- vector()
  for (i in 1:(nrow(clusters) - 1)) {
    for (j in (i + 1):nrow(clusters)) {
      if (j %in% to_remove) next
      sim <- jaccard_similarity(clusters$IDs[[i]], clusters$IDs[[j]])
      if (sim > threshold) {
        to_remove <- c(to_remove, j)
      }
    }
  }
  filtered_clusters <- clusters[-to_remove,]
  if (nrow(filtered_clusters) <= 8) {
    message("Filtered clusters are less than or equal to 8, moving to the next iteration.")
    Cor_list[[x]] <- data.frame(module = NA, cor = 0) 
    next
  }
  #过滤模块
  ModuleResult <- ModuleResult[ModuleResult$cluster%in%filtered_clusters$cluster,]
  #计算模块水平
  datax=data$PathwaysMatrix[,ModuleResult$ID]
  datax=moduleEigengenes(datax,ModuleResult$cluster)$eigengenes
  #计算模块相关性  
  datax$label=data$PathwaysMatrix[,'label']
  datax=datax[,c(ncol(datax),(1:(ncol(datax)-1)))]
  Cor=data.frame(module=rownames(as.data.frame(stats::cor(datax)))[-1],cor=as.data.frame(stats::cor(datax))[-1,1])
  #通过MEtrait构造DE_PathwaysModule
  vp.counts.df <- data.frame(
    module = paste0("ME",colnames(S.vp.split)),
    Num_pathways = as.vector(apply(S.vp.split, 2, sum)),
    stringsAsFactors = FALSE
  )
  p.values <- sapply(2:ncol(datax), function(i) {
    cor.test(datax[, 1], datax[, i])$p.value
  })
  Cor_P <- data.frame(module = colnames(datax)[-1], pvalue = p.values)
  DE_PathwaysModule <- Cor%>%
    mutate(Fraction=80)%>%
    merge(.,Cor_P,by="module")%>%
    mutate(adjust_pvalue=p.adjust(pvalue, method = "fdr"))%>%
    merge(.,vp.counts.df,by="module")%>%
    dplyr::select(1,6,3,4,5,2)%>%
    arrange(desc(cor))%>%
    slice_head(n=12)
  #构造Matrix
  Matrix=data$PathwaysMatrix
  #组合成MTG_Modules
  brainModules <- list(
    ModuleResult = ModuleResult,
    DE_PathwaysModule = DE_PathwaysModule,
    Matrix = Matrix
  )
  return(brainModules)
}
##汇总函数
construct_Modules <- function(PATH, num){
  test <- min_modules(PATH,num)
  num.comp <- test[["num.comp"]]
  Modules <- brain_Modules(PATH,num.comp)
  return(Modules)
}
#运行示例
library(DEXICA)
library(dplyr)
library(WGCNA)
num=seq(10,400,20)
aa <- construct_Modules(MTG_PATH,num)
MTG_Modules <- aa


