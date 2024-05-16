#####Focused array plot####
#####1.read Loess normalized data####
setwd('D:/LC autoantibodies/LC_diagnosis autoantibodies/LUADlowdensity/fit data for plot and plot')
LD.IgG<-read.table('./normalization data of IgG after impute add group.csv',sep = ',',head=T)
LD.IgM<-read.table('./normalization data of IgM after impute add group.csv',sep = ',',head=T)
filt<-read.table('./select IgG and IgM.csv',sep=',',header = T)

library(scatterplot3d)
library(ggsci)
library(plotly)
library(RColorBrewer)
library(dplyr)

#####Scatterplot3d used origin software with manual operation#####


##### Candidates boxviolin plot####
library(tidyr)
library(dplyr)
library(ggplot2)
library(ggsci)
library(ggpubr)
library(purrr)
library(stringr)
setwd('D:/LC autoantibodies/LC_diagnosis autoantibodies/LUADlowdensity/fit data for plot and plot')
LD.IgG<-read.table('./normalization data of IgG after impute add group.csv',sep = ',',head=T)
LD.IgM<-read.table('./normalization data of IgM after impute add group.csv',sep = ',',head=T)
filt<-read.table('./select IgG and IgM.csv',sep=',',header = T)

filtID<-unique(filt$ID)
filtName<-unique(filt$featureName)

#IgG
datIgG<-LD.IgG%>%select(c(1:3,colnames(LD.IgG)[colnames(LD.IgG)%in%filtID]))
colnames(datIgG)[-c(1:3)]<-filtName[match(colnames(datIgG[-c(1:3)]),filtID)]
datIgGL<-gather(datIgG,key='AAbs',value = "Fluorescene", -c(1:3))#长转宽
datIgGL$Group1<-factor(datIgGL$Group1,levels = c('EarlyLUAD','Control'))
datIgGL$Group2<-factor(datIgGL$Group2,levels = c('EarlyLUAD','BLD','NHC'))
#IgM
datIgM<-LD.IgM%>%select(c(1:3,colnames(LD.IgM)[colnames(LD.IgM)%in%filtID]))
colnames(datIgM)[-c(1:3)]<-filtName[match(colnames(datIgM[-c(1:3)]),filtID)]
datIgML<-gather(datIgM,key='AAbs',value = "Fluorescene", -c(1:3))#长转宽
datIgML$Group1<-factor(datIgML$Group1,levels = c('EarlyLUAD','Control'))
datIgML$Group2<-factor(datIgML$Group2,levels = c('EarlyLUAD','BLD','NHC'))


cp=c('EarlyLUAD vs BLD','EarlyLUAD vs NHC','EarlyLUAD vs Control')
cc=list(c('EarlyLUAD','BLD'),c('EarlyLUAD','NHC'),c('EarlyLUAD','Control'))
#IgG
selectAAbs<-map(cp,function(x)filt[filt$Comparison==x&filt$Type=='IgG',2])
selectAbdat<-map2(selectAAbs,cc,function(x,y)if(all(y==cc[[3]])) {
  datIgGL[datIgGL$AAbs%in%x&datIgGL$Group1%in%y,c(2,4,5)]
} else {
  datIgGL[datIgGL$AAbs%in%x&datIgGL$Group2%in%y,3:5]
})

#IgM
selectAAbs2<-map(cp,function(x)filt[filt$Comparison==x&filt$Type=='IgM',2])
selectAbdat2<-map2(selectAAbs2,cc,function(x,y)if(all(y==cc[[3]])) {
  datIgML[datIgML$AAbs%in%x&datIgML$Group1%in%y,c(2,4,5)]
} else {
  datIgML[datIgML$AAbs%in%x&datIgML$Group2%in%y,3:5]
})


#  initial submission version  --------------------------------------------------------------
viobox=function(dat){
  if(names(dat)[1]=='Group2'){
    ggplot(dat, aes(AAbs, Fluorescene,color=Group2))+ 
      geom_jitter(alpha=0.32,
                  position=position_jitterdodge(jitter.width = 0.2, 
                                                jitter.height = 0, 
                                                dodge.width = 0.8))+
      geom_boxplot(alpha=0.33,width=0.3,
                   position=position_dodge(width=0.8),
                   size=0.6,outlier.colour = NA)+
      geom_violin(alpha=0.4,width=1.2,
                  position=position_dodge(width=0.8),
                  size=0.75)+
      scale_color_manual(values = unique(pal_lancet("lanonc")(9)[c(2,1)]))+
      theme_classic()+
      theme(panel.background=element_rect(fill="white",colour="black",size=0.25),
            axis.line=element_line(colour="black",size=0.25),
            axis.title=element_text(size=13,face="plain",color="black"),
            axis.text = element_text(size=10,face="plain",color="black"),
            legend.position='right', 
            legend.title = element_text(size = 9), 
            legend.text = element_text(size = 9),
            legend.key.width = unit(0.8, "cm"),
            legend.key.height = unit(0.8, "cm"),
            axis.text.x= element_text(angle = 45,hjust = 0.9))+
      stat_compare_means(method = 't.test', aes(label = paste0(..p.signif..)),show.legend = F)
  }else{ ggplot(dat, aes(AAbs, Fluorescene,color=Group1))+ 
      geom_jitter(alpha=0.32,
                  position=position_jitterdodge(jitter.width = 0.2, 
                                                jitter.height = 0, 
                                                dodge.width = 0.8))+
      geom_boxplot(alpha=0.33,width=0.3,
                   position=position_dodge(width=0.8),
                   size=0.6,outlier.colour = NA)+
      geom_violin(alpha=0.4,width=1.2,
                  position=position_dodge(width=0.8),
                  size=0.75)+
      scale_color_manual(values = unique(pal_lancet("lanonc")(9)[c(2,1)]))+
      theme_classic()+
      theme(panel.background=element_rect(fill="white",colour="black",size=0.25),
            axis.line=element_line(colour="black",size=0.25),
            axis.title=element_text(size=13,face="plain",color="black"),
            axis.text = element_text(size=10,face="plain",color="black"),
            legend.position='right', 
            legend.title = element_text(size = 9), 
            legend.text = element_text(size = 9),
            legend.key.width = unit(0.8, "cm"),
            legend.key.height = unit(0.8, "cm"),
            axis.text.x= element_text(angle = 45,hjust = 0.9))+
      stat_compare_means(method = 't.test', aes(label = paste0(..p.signif..)),show.legend = F)
  }
  }


files = str_c(str_c('./vioboxplot/vioboxplot of selected AAbs in FA for',cp,'.IgG'), ".pdf")
walk2(files,map(selectAbdat,viobox),~ggsave(.x,.y,width=16,height = 6,units = 'in',dpi = 300))


files2 = str_c(str_c('./vioboxplot/vioboxplot of selected AAbs in FA for',cp,'.IgM'), ".pdf")
walk2(files2,map(selectAbdat2,viobox),~ggsave(.x,.y,width=16,height = 6,units = 'in',dpi = 300))


# revision version  --------------------------------------------------------------
library(ggbeeswarm)
library(ggpubr)
library(ggnewscale)
names(selectAbdat)<-cp
names(selectAbdat2)<-cp
selectAbdat<-map2(selectAbdat,list(c('EarlyLUAD','BLD'),c('EarlyLUAD','NHC'),c('EarlyLUAD','Control')),
                 function(x,y){colnames(x)[1]<-'Group'
                             x<-x%>%mutate(factor(Group,levels=y))
                             return(x)})
selectAbdat2<-map2(selectAbdat2,list(c('EarlyLUAD','BLD'),c('EarlyLUAD','NHC'),c('EarlyLUAD','Control')),
                  function(x,y){colnames(x)[1]<-'Group'
                  x<-x%>%mutate(factor(Group,levels=y))
                  return(x)})
#mean sd
selectAbdat_stat <- purrr::map2(selectAbdat,list(c('EarlyLUAD','BLD'),c('EarlyLUAD','NHC'),c('EarlyLUAD','Control')),
                                 function(x,y){
                                   dfstat<- x %>% 
                                     group_by(AAbs,Group) %>% 
                                     mutate(upper =  quantile(Fluorescene, 0.75),
                                            lower = quantile(Fluorescene, 0.25),
                                            mean = mean(Fluorescene),
                                            median = median(Fluorescene),
                                            sd = sd(Fluorescene))%>%mutate(Group=factor(Group,levels=y))
                                   return(dfstat)
                                 })

selectAbdat_stat2 <- purrr::map2(selectAbdat2,list(c('EarlyLUAD','BLD'),c('EarlyLUAD','NHC'),c('EarlyLUAD','Control')),
                                function(x,y){
                                  dfstat<- x %>% 
                                    group_by(AAbs,Group) %>% 
                                    mutate(upper =  quantile(Fluorescene, 0.75),
                                           lower = quantile(Fluorescene, 0.25),
                                           mean = mean(Fluorescene),
                                           median = median(Fluorescene),
                                           sd = sd(Fluorescene))%>%mutate(Group=factor(Group,levels=y))
                                  return(dfstat)
                                })

pvalue<-map(selectAbdat_stat,function(x){
  filt$P.Value[match(unique(x$AAbs),filt$featureName)]
})
pvalue2<-map(selectAbdat_stat2,function(x){
  filt$P.Value[match(unique(x$AAbs),filt$featureName)]
})


vioscatter <- function(df, dfstatics, x, y, group, ylab, pval) {
  p1 <- ggplot(df, aes(x = .data[[x]], y = .data[[y]])) +
    geom_violin(
      aes(fill = .data[[group]]),
      position = position_dodge(width = 0.8),
      color = "grey",
      trim = T,
      adjust = 3
    ) +
    scale_fill_manual(values = c('#EDEDED', '#EDEDED')) +
    geom_quasirandom(
      aes(color = .data[[group]]),
      width = 0.15,
      dodge.width = 0.8,
      size = 0.4,
      alpha = 0.7
    ) +
    scale_color_manual(values = c('#CD5B45', '#4682B4'),
                       labels = unique(df[, group])) +
    theme_classic() +
    labs(x = " ", y = ylab) +
    theme(
      axis.title.y = element_text(colour = 'black', size = 16),
      axis.text = element_text(colour = 'black', size = 14),
      axis.text.x = element_text(
        angle = 45,
        vjust = 0.7,
        hjust = 0.8
      ),
      axis.line = element_line(size = 1),
      legend.title = element_text(size = 14),
      legend.text = element_text(size = 14)
    ) +
    guides(color = guide_legend(override.aes = list(size = 4)))
  
  p2 <- p1 + new_scale_color() +
    geom_errorbar(
      data = dfstatics,
      aes(
        x = .data[[x]],
        color = .data[[group]],
        ymin = mean - sd,
        ymax = mean + sd
      ),
      position = position_dodge(width = 0.8),
      width = 0.3,
      size = 0.7
    ) +
    scale_color_manual(values = c('black', 'black')) +
    new_scale_fill() +
    stat_summary(
      fun = "mean",
      geom = "crossbar",
      linewidth = 0.28,
      mapping = aes(
        ymin = ..y..,
        ymax = ..y..,
        group = .data[[group]]
      ),
      position = position_dodge(width = 0.8),
      size = 0.20
    ) +
    stat_summary(
      aes(fill = .data[[group]]),
      geom = "point",
      position = position_dodge(width = 0.8),
      fun = mean,
      shape = 21,
      size = 2.5,
      stroke = 0.75
    ) +
    scale_fill_manual(values = c('#E65454', '#36649B')) +
    geom_signif(
      xmin = seq(0.8, length(unique(df[, x])) - 0.2),
      xmax = seq(1.2, length(unique(df[, x])) + 0.2),
      annotations = sapply(pval, function(p) {
        ifelse(p < 0.001, '***', ifelse(p < 0.01, '**', '*'))
      }),
      y_position = rep(max(df[, y]), length(unique(df[, x]))),
      textsize = 5,
      tip_length = c(0, 0),
      size = 0.5
    ) +
    guides(fill = guide_legend(title = 'Mean'))
  p2
}

files = str_c(str_c('./vioboxplot/vioScatter of selected AAbs in FA for',cp,'.IgG'), ".pdf")
walk2(files,pmap(list(selectAbdat, selectAbdat_stat, pvalue), function(x,y,z)vioscatter(df=x,dfstatics = y,x='AAbs',
                                                                                        y='Fluorescene',group='Group',ylab='Fluorescene Intensity',pval=z)),
            ~ggsave(.x,.y,width=14,height = 8,units = 'in',dpi = 300))


files2 = str_c(str_c('./vioboxplot/vioScatter of selected AAbs in FA for',cp,'.IgM'), ".pdf")
walk2(files2,pmap(list(selectAbdat, selectAbdat_stat, pvalue), function(x,y,z)vioscatter(df=x,dfstatics = y,x='AAbs',
                                                                                        y='Fluorescene',group='Group',ylab='Fluorescene Intensity',pval=z)),
      ~ggsave(.x,.y,width=14,height = 8,units = 'in',dpi = 300))


