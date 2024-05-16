#changelog
#2024.3.25
#delete the bottom plot of Heatmap for concise display according to the reviewer's question

####
library(dplyr)
library(readr)
####1.IgG,IgM data after normalization####
setwd('D:/LC autoantibodies/LC_diagnosis autoantibodies/112samples')
#filter proteins based on FC, P and sensitivity
Loessfilt<-read.table('./filt data for plot/Loess select.csv',sep=',',header=T)
#read nomalized IgG data
Lf.IgG<-Loessfilt%>%filter(Type=='IgG')
norm.IgG<-read.table('./Loess/originalData/highdensity/normalization.IgG/normdata2del.csv',sep=',',header=T)
Lfdat.IgG<-Lf.IgG[1:5]%>%left_join(norm.IgG[,-c(1:3,6)],by='ID')  
#setNames for proteins with different ID but same name 
Lfdat.IgG$Name.x[duplicated(Lfdat.IgG$Name.x)]<-c("IRAK4.1","TIPIN.1","TIPIN.2","BIN1.1","ZGRF1.1","BIN1.2","PHF10.1")

#read nomalized IgM data
Lf.IgM<-Loessfilt%>%filter(Type=='IgM')
norm.IgM<-read.table('./Loess/originalData/highdensity/normalization.IgM/normdata2del.IgM.csv',sep=',',header=T)
Lfdat.IgM<-Lf.IgM[1:5]%>%left_join(norm.IgM[,-c(1:3,6)],by='ID')  
#setNames for proteins with different ID but same name
Lfdat.IgM$Name.x[duplicated(Lfdat.IgM$Name.x)]<-c("KIZ.1","TAGLN2.1","EZR.1")



##read clinical information 
clininfo.IgG<-read_csv("filt data for plot/Clinic info for HuProt_IgG.csv", 
                       locale = locale(encoding = "GB2312"))
Smoke.IgG<-ifelse(clininfo.IgG$Smoke=='Yes','Y',ifelse(clininfo.IgG$Smoke=='No','N',''))
Smoke.IgG2<-replace(Smoke.IgG, is.na(Smoke.IgG), "")
Alcohol.IgG<-ifelse(clininfo.IgG$Alcohol=='Yes','Y',ifelse(clininfo.IgG$Alcohol=='No','N',''))
Alcohol.IgG2<-replace(Alcohol.IgG, is.na(Alcohol.IgG), "")

clininfo.IgM<-read_csv('./filt data for plot/Clinic info for HuProt_IgM.csv',locale = locale(encoding = "GB2312"))
Smoke.IgM<-ifelse(clininfo.IgM$Smoke=='Yes','Y',ifelse(clininfo.IgM$Smoke=='No','N',''))
Smoke.IgM2<-replace(Smoke.IgM, is.na(Smoke.IgM), "")
Alcohol.IgM<-ifelse(clininfo.IgM$Alcohol=='Yes','Y',ifelse(clininfo.IgM$Alcohol=='No','N',''))
Alcohol.IgM2<-replace(Alcohol.IgM, is.na(Alcohol.IgM), "")

####2.extract the IgG,IgM data of proteins filtered by supplementary methods####
#proteins filtered by supplementary methods
supfilt<-read_csv('./filt data for plot/supplemtal methods select.csv',locale = locale(encoding = "GB2312"))
#Zscore normalized data
Zsnorm.IgG<-read.table('./Zscore/Zscorenormalization/IgG/allSamples_mergedheadfix.norm.txt',sep='\t',header=T,fileEncoding = 'GB2312')
Zsnorm.IgM<-read.table('./Zscore/Zscorenormalization/IgM/allSamples_mergedheadfix.norm.IgM.csv',sep=',',header=T,fileEncoding = 'GB2312')

sf.IgG<-supfilt%>%filter(ID%in%c('JHU04931.B3C28','JHU04987.B3C7'))#EarlyLUAD_NHC,GIMAP4,WWP2
sfdat.IgG<-sf.IgG[1:5]%>%left_join(Zsnorm.IgG[-c(1:3)],by='ID')

sf.IgM<-supfilt%>%filter(ID%in%c('JHU02685.B4C28','JHU06755.B7C6'))#EarlyLUAD_NHC,DCTPP1,GDA
sfdat.IgM<-sf.IgM[1:5]%>%left_join(Zsnorm.IgM[-c(1:3)],by='ID')


######integrate data of Loess normalized and supplementary method#####
r1<-81 
dat.IgG<-rbind(Lfdat.IgG[1:r1,],sfdat.IgG,Lfdat.IgG[-(1:r1),])

r2<-70 
dat.IgM<-rbind(Lfdat.IgM[1:r2,],sfdat.IgM,Lfdat.IgM[-(1:r2),])




#####3.plot complex Heatmap####
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(ggplot2)

mat.IgG<-as.matrix(dat.IgG[,7:ncol(dat.IgG)])
rownames(mat.IgG)<-dat.IgG$Name.x
mat.IgGs<-apply(mat.IgG, 1, scale)##row scale
rownames(mat.IgGs) = colnames(mat.IgG)
mat.IgGs = t(mat.IgGs)

mat.IgM<-as.matrix(dat.IgM[,7:ncol(dat.IgM)])
rownames(mat.IgM)<-dat.IgM$Name.x
mat.IgMs<-apply(mat.IgM, 1, scale)##row scale
rownames(mat.IgMs) = colnames(mat.IgM)
mat.IgMs = t(mat.IgMs)

#####3.1IgG data Heatmap####
f = colorRamp2(seq(min(mat.IgGs), max(mat.IgGs), length = 3), c("blue", "#EEEEEE", "red"), space = "RGB")
ht_opt$TITLE_PADDING = unit(c(1.5, 1.5), "points")## 默认情况下标题顶部没有空间，现在我们增加4pt的空间

colGender=brewer.pal(3,'Set2')[c(1,2)]
names(colGender)=c('Female','Male')
colStage=brewer.pal(7,'Set1')[c(2,5,1)]
names(colStage)=c('Stage 0','Stage IA','Stage IB')

######3.11top annototation####
hatop=HeatmapAnnotation(
        df=data.frame(Gender=clininfo.IgG$Gender,
                      Stage=clininfo.IgG$Stage),
        col = list(Gender = colGender,
                   Stage = colStage),
        #Age=anno_points(clininfo$Age,pch = 16,gp = gpar(col= 'red'),size = unit(1.0, "mm")),
        Age=anno_barplot(clininfo.IgG$Age,gp=gpar(fill='#CDB7B5'),height = unit(0.4, "cm")),
        Smoke=anno_text(Smoke.IgG2,location = 0.5,just = "center",rot =360,
                        gp = gpar(fill=ifelse(Smoke.IgG2=='Y','#F08080',ifelse(Smoke.IgG2=='N','#B4CDCD','lightgrey')),
                                  fontsize =5), 
                        width = max_text_width(Smoke.IgG2) * 1.5,
                        height = unit(0.28, "cm")),
        Alcohol=anno_text(Alcohol.IgG2,location = 0.5,just = "center",rot =360,
                          gp = gpar(fill=ifelse(Alcohol.IgG2=='Y','#F08080',ifelse(Alcohol.IgG2=='N','#B4CDCD','lightgrey')),
                                    fontsize =5), 
                          width = max_text_width(Alcohol.IgG) * 1.5,
                          height = unit(0.28, "cm")),
        simple_anno_size = unit(0.15, "cm"), 
        #gp = gpar(col = "black"),
        height = unit(1.3, "cm"),
        border = c(Gender= F,Stage=F))
draw(hatop)

#######3.12Bottom annototation####
# fP = colorRampPalette(c('green','orange','red'))(88)
# fI = colorRampPalette(c('green','orange','red'))(88)
# fL = colorRampPalette(c('blue','yellow','red'))(88)
# 
# habottom=columnAnnotation(PMD = anno_barplot(clininfo$PMD,gp=gpar(fill=fP[rank(clininfo$PMD[1:88])])),
#                           IMD=anno_barplot(clininfo$IMD,gp=gpar(fill=fI[rank(clininfo$IMD[1:88])])),
#                           LN=anno_barplot(clininfo$LN,gp=gpar(fill=fL[rank(clininfo$LN[1:88])])),
#                           height = unit(2.5, "cm"))

########3.13Left annototation####
r1<-81 
Se.IgG<-c(Loessfilt$Se[1:r1],sf.IgG$Se,Loessfilt$Se[(r1+1):227])
FC.IgG<-c(Loessfilt$logFC[1:r1],sf.IgG$logFC,Loessfilt$logFC[(r1+1):227])
fSe = colorRamp2(seq(min(Se.IgG), max(Se.IgG), length = 3), c('#F5FFFA','#EEA2AD','#F08080'), space = "RGB")
haright=rowAnnotation(#logFC=anno_points(FC.IgG,pch = 20,gp = gpar(col= rainbow(12))),
                     Se=Se.IgG,col=list(Se=fSe),border=T,  
                     simple_anno_size = unit(0.3, "cm"),
                     width = unit(0.85, "cm"))
draw(haright)

#########3.14 IgG final Heatmap####
dev.new()
ht<-Heatmap(mat.IgGs,col = f,border_gp = gpar(lty = 2, col = "red"),
        rect_gp = gpar(col = "white", lty = 0.5, lwd = 1),name = "Fluorescene",
        column_title = c('EarlyLUAD','BLD','NHC'),
        row_title= c('EarlyLUAD vs BLD','EarlyLUAD vs NHC','EarlyLUAD vs Control'),
        column_title_gp = gpar(fill = c("#CD5555", "#3A5FCD", "#3CB371"), font = 1),
        row_title_gp = gpar(font=0.5, fill = c("#CD5555", "#3A5FCD", "#3CB371")),
        row_names_gp = gpar(fontsize =2.5),
        column_names_gp = gpar(fontsize =3.5,col = c(rep("#CD5555", 58), rep("#3A5FCD", 30),rep('#3CB371',24))),
        cluster_rows = F, cluster_columns = F,
        row_names_centered = TRUE, column_names_centered = TRUE,column_names_rot = 45,
        row_split = factor(c(rep("EarlyLUAD vs BLD", 81),rep("EarlyLUAD vs NHC", 105),rep("EarlyLUAD vs Control", 43)),
                           levels = c("EarlyLUAD vs BLD","EarlyLUAD vs NHC","EarlyLUAD vs Control")),
        column_split = factor(c(rep("EarlyLUAD", 58),rep("BLD", 30),rep("NHC", 24)),
                              levels = c('EarlyLUAD','BLD','NHC')),
        top_annotation =hatop,
        #bottom_annotation = habottom,
        right_annotation = haright,
        #heatmap_width = unit(12, "cm"),heatmap_height =unit(15, "cm")
        )
pdf('D:/LC autoantibodies/LC_diagnosis autoantibodies/112samples/heatmap/IgG heatmap in HuProt.pdf',height = 10,width =8)
draw(ht)
dev.off()



#####3.2IgM Heatmap####
f = colorRamp2(seq(min(mat.IgMs), max(mat.IgMs), length = 3), c("blue", "#EEEEEE", "red"), space = "RGB")
ht_opt$TITLE_PADDING = unit(c(1.5, 1.5), "points")## 默认情况下标题顶部没有空间，现在我们增加4pt的空间

colGender=brewer.pal(3,'Set2')[c(1,2)]
names(colGender)=c('Female','Male')
colStage=brewer.pal(7,'Set1')[c(2,5,1)]
names(colStage)=c('Stage 0','Stage IA','Stage IB')

####top annototation####
hatop=HeatmapAnnotation(
  df=data.frame(Gender=clininfo.IgM$Gender,
                Stage=clininfo.IgM$Stage),
  col = list(Gender = colGender,
             Stage = colStage),
  #Age=anno_points(clininfo$Age,pch = 16,gp = gpar(col= 'red'),size = unit(1.0, "mm")),
  Age=anno_barplot(clininfo.IgM$Age,gp=gpar(fill='#CDB7B5'),height = unit(0.4, "cm")),
  Smoke=anno_text(Smoke.IgM2,location = 0.5,just = "center",rot =360,
                  gp = gpar(fill=ifelse(Smoke.IgM2=='Y','#F08080',ifelse(Smoke.IgM2=='N','#B4CDCD','lightgrey')),
                            fontsize =5), 
                  width = max_text_width(Smoke.IgM2) * 1.5,
                  height = unit(0.28, "cm")),
  Alcohol=anno_text(Alcohol.IgM2,location = 0.5,just = "center",rot =360,
                    gp = gpar(fill=ifelse(Alcohol.IgM2=='Y','#F08080',ifelse(Alcohol.IgM2=='N','#B4CDCD','lightgrey')),
                              fontsize =5), 
                    width = max_text_width(Alcohol.IgM) * 1.5,
                    height = unit(0.28, "cm")),
  simple_anno_size = unit(0.15, "cm"), 
  #gp = gpar(col = "black"),
  height = unit(1.3, "cm"),
  border = c(Gender= F,Stage=F))
draw(hatop)

####Bottom annototation####
# fP = colorRampPalette(c('green','orange','red'))(88)
# fI = colorRampPalette(c('green','orange','red'))(88)
# fL = colorRampPalette(c('blue','yellow','red'))(88)
# 
# habottom=columnAnnotation(PMD = anno_barplot(clininfo.IgM$PMD,gp=gpar(fill=fP[rank(clininfo.IgM$PMD[1:88])])),
#                           IMD=anno_barplot(clininfo.IgM$IMD,gp=gpar(fill=fI[rank(clininfo.IgM$IMD[1:88])])),
#                           LN=anno_barplot(clininfo.IgM$LN,gp=gpar(fill=fL[rank(clininfo.IgM$LN[1:88])])),
#                           height = unit(2.5, "cm"))
# draw(habottom)

####Right annototation####
Se.IgM<-c(Loessfilt$Se[228:296],sf.IgM$Se,Loessfilt$Se[297:387])
FC.IgM<-c(Loessfilt$logFC[228:296],sf.IgM$logFC,Loessfilt$logFC[297:387])

fSe = colorRamp2(seq(min(Se.IgM), max(Se.IgM), length = 3),  c('#F5FFFA','#EEA2AD','#F08080'), space = "RGB")
haright=rowAnnotation(#logFC=anno_points(Se.IgM,pch = 20,gp = gpar(col= rainbow(12))),
                      Se=FC.IgM,col=list(Se=fSe),border=T,
                      simple_anno_size = unit(0.3, "cm"),
                      width = unit(0.85, "cm"))
draw(haright)

#########3.14 IgM final Heatmap ####
dev.new()
ht<-Heatmap(mat.IgMs,col = f,border_gp = gpar(lty = 2, col = "red"),
            rect_gp = gpar(col = "white", lty = 0.5, lwd = 1),name = "Fluorescene",
            column_title = c('EarlyLUAD','BLD','NHC'),
            row_title= c('EarlyLUAD vs BLD','EarlyLUAD vs NHC','EarlyLUAD vs Control'),
            column_title_gp = gpar(fill = c("#CD5555", "#3A5FCD", "#3CB371"), font = 1),
            row_title_gp = gpar(font=0.5, fill = c("#CD5555", "#3A5FCD", "#3CB371")),
            row_names_gp = gpar(fontsize =2.5),
            column_names_gp = gpar(fontsize =3.5,col = c(rep("#CD5555", 58), rep("#3A5FCD", 30),rep('#3CB371',24))),
            cluster_rows = F, cluster_columns = F,
            row_names_centered = TRUE, column_names_centered = TRUE,column_names_rot = 45,
            row_split = factor(c(rep("EarlyLUAD vs BLD", 69),rep("EarlyLUAD vs NHC", 72),rep("EarlyLUAD vs Control", 21)),#可以查看dat.IgM数据看分组数量
                               levels = c("EarlyLUAD vs BLD","EarlyLUAD vs NHC","EarlyLUAD vs Control")),
            column_split = factor(c(rep("EarlyLUAD", 38),rep("BLD", 30),rep("NHC", 24)),
                                  levels = c('EarlyLUAD','BLD','NHC')),
            top_annotation =hatop,
            #bottom_annotation = habottom,
            right_annotation = haright,
            #heatmap_width = unit(12, "cm"),heatmap_height =unit(15, "cm")
)
pdf('D:/LC autoantibodies/LC_diagnosis autoantibodies/112samples/heatmap/IgM heatmap in HuProt.pdf',height = 10,width =8)
draw(ht)
dev.off()
