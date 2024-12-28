#绘图
#load("Neuron_Modules.RData")
Modules=froCortex_Modules
froCortex_Modules$DE_PathwaysModule$module[1:4]
modules=c(2,14,3,8)
#小提琴,pdf尺寸5，7
img=list()
for(i in 1:length(modules)){
  pic=VisMultiModule(PathwaysModule_obj=Modules,volin=TRUE,control_label=0,module=modules[i])
  img[[i]]=pic
}
d=0.1
img[[1]]+theme(plot.margin = unit(c(d,d,d,d), "cm"))+
  img[[2]]+theme(plot.margin = unit(c(d,d,d,d), "cm"))+
  img[[3]]+theme(plot.margin = unit(c(d,d,d,d), "cm"))+
  img[[4]]+theme(plot.margin = unit(c(d,d,d,d), "cm"))
#词云图
output=paste0('Module',modules,'_WordCloud.png')
for(i in 1:length(modules)){
  ModulesInner = ShowModule(Modules,modules[i],exact=F)
  my_graph=VisMultiModule(ShowModule_obj=ModulesInner,exact =F)
  saveWidget(my_graph,'tmp.html',selfcontained = F)
  webshot('tmp.html',output[i])
}
ff=list.files(pattern = 'png')
gg=do.call(cbind,lapply(1:4,function(x) image_read(ff[x])))
ff2=gsub('_WordCloud.png','',ff)
ff2=gsub('ME','Module',ff2)
ba='#FFFFFF'
d=2
image_ggplot(gg[[1]])+labs(title = ff2[1])+
  theme(text = element_text(family = "serif", size = 10, color = "black",face='bold'),
        panel.border = element_rect(color = "black", fill = NA, size = 2),
        plot.margin = unit(c(d, d, d, d), "mm"),
        panel.background = element_rect(fill = ba, color = ba),
        plot.background = element_rect(fill = ba, color = ba))+
  image_ggplot(gg[[2]])+labs(title = ff2[2])+
  theme(text = element_text(family = "serif", size = 10, color = "black",face='bold'),
        panel.border = element_rect(color = "black", fill = NA, size = 2),
        plot.margin = unit(c(d, d, d, d), "mm"),
        panel.background = element_rect(fill = ba, color = ba),
        plot.background = element_rect(fill = ba, color = ba))+
  image_ggplot(gg[[3]])+labs(title = ff2[3])+
  theme(text = element_text(family = "serif", size = 10, color = "black",face='bold'),
        panel.border = element_rect(color = "black", fill = NA, size = 2),
        plot.margin = unit(c(d, d, d, d), "mm"),
        panel.background = element_rect(fill = ba, color = ba),
        plot.background = element_rect(fill = ba, color = ba))+
  image_ggplot(gg[[4]])+labs(title = ff2[4])+
  theme(text = element_text(family = "serif", size = 10, color = "black",face='bold'),
        panel.border = element_rect(color = "black", fill = NA, size = 2),
        plot.margin = unit(c(d, d, d, d), "mm"),
        panel.background = element_rect(fill = ba, color = ba),
        plot.background = element_rect(fill = ba, color = ba))
######ROC曲线
library(pROC)
load("~/zhutianshu/combined/BioM2/Glia_ROC.RData")
a=Glia_ROC[["Prediction"]]
combined_df <- data.frame(sample = character(), prediction = numeric())
for (i in 1:length(a)) {
  df <- a[[i]]
  combined_df <- rbind(combined_df, df)  # 合并数据框
}
combined_df <- combined_df%>%
  dplyr::rename(Sample_Name=sample)%>%
  merge(GSE66531_samplesheet[,c(2,5)],.,by="Sample_Name")

roc1 <- roc(combined_df$label, combined_df$prediction, levels=c("0","1"))
#先绘制1条ROC曲线
#pdf("roc_curve_output.pdf", width = 6, height = 6)
#par(pty="s",las=1)
plot(roc1,
     print.auc=F,print.auc.cex=0.85,
     print.auc.x=0.7,print.auc.y=0.4,
     auc.polygon=TRUE,auc.polygon.col="#fff7f7",
     grid=c(0.2,0.2),grid.col=c("skyblue","skyblue"),
     print.thres=F,print.thres.cex=0.85,
     main="ROC curves", col='tomato',
     legacy.axes=TRUE)
ci_obj <- ci.auc(roc1)
text(0.5, 0.2, paste("AUC:", round(roc1$auc, 2), " [", round(ci_obj[1], 2), "; ", round(ci_obj[3], 2), "]", sep=""),
     cex=0.85, col="tomato")