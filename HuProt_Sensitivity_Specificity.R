###Calculate Sensitivity,Specificity,AUC based on the Loess Normalized matrix
library(impute)
normdata0<-read.table('K:/LC_diagnosis autoantibodies/112samples/Loess/originalData/highdensity/normalization.IgM/normalizedData_fastCyclicloess.IgM.csv',sep=',',header = T)
normdata1<-impute.knn(as.matrix(normdata0[,-c(1:6)])) #impute NA 
normdata2<-cbind(normdata0[,1:6],as.data.frame(normdata1[1]))
colnames(normdata2)<-colnames(normdata0)
write.table(normdata2,'./normalization.IgM/normdata2del.csv',sep = ',')
##extract NHC data
normdat<-as.data.frame(t(normdata2[,7:98])) #7:98 belong to NHC
colnames(normdat)<-normdata2$ID
normdat<-cbind(Sample=colnames(normdata2[7:98]),normdat)
###add group info
normdat1<-cbind(group=c(rep('EarlyLUAD',38),rep('BLD',30),rep('NHC',24)),normdat)
normdat1$group<-as.factor(normdat1$group)
###calculate the threshold of mean+2SD in NHC
fun<-function(x)mean(x)+2*sd(x)
cutoff<-apply(normdat1[normdat1$group=='NHC',-c(1:2)],2,FUN = fun)
#1.transform the expression data into 0 or 1 compared to the threshold
normdat1.t<-t(normdat1[,-c(1:2)])
compare<-ifelse(normdat1.t>cutoff,1,0)
compare.t<-t(compare)
dat.trans<-cbind(normdat1[1:2],as.data.frame(compare.t))
dat.trans.t<-cbind(normdata2[,4:5],as.data.frame(t(dat.trans[,-(1:2)]))) ###row as the proteins, column as samples
write.table(dat.trans.t,'K:/LC_diagnosis autoantibodies/112samples/Loess/originalData/highdensity/Se.Sp.IgM/dat.pos or neg.IgM.csv',sep = ',')
#2.calculate sensitivity and specificity
library(dplyr)
dat.trans[,-2]%>%group_by(group)%>%summarise_if(.,is.numeric,sum)->pos.count
Se<-pos.count[2,-1]/38 ##EarlyLUAD
Se1<-pos.count[1,-1]/30##BLD
Sp1<-(30-pos.count[1,-1])/30 ###EarlyLUAD vs BLD
Sp2<-(24-pos.count[3,-1])/24 ###EarlyLUAD vs NHC
Sp3<-(54-(pos.count[1,-1]+pos.count[3,-1]))/54
Se.Sp<-as.data.frame(t(rbind(Se,Se1,Sp1,Sp2,Sp3)))
colnames(Se.Sp)<-c('Se','Se1','Sp1','Sp2','Sp3')
write.table(Se.Sp,'./Se.Sp.IgM/Se.Sp_IgM.csv',sep = ',')


#3.Combine the filter criteria of logFC, P value and sensitivity,specificity
##3.1 read the differential result of LoessNormalized data£¬filter FC>1 and P<0.05
library(dplyr)
con<-'D:/LC autoantibodies/LC_diagnosis autoantibodies/Loess/originalData/highdensity/differentialAnalysis.IgG'
Ben_Nor<-read.table(paste0(con,'/significantDifferentialTest_BenignvsNormal.csv'),sep=',',header = T)
Ben_Nor<-rename(Ben_Nor,ID=X)
ELC_Ben<-read.table(paste0(con,'/significantDifferentialTest_Early_LCvsBenign.csv'),sep=',',header = T)
ELC_Ben<-rename(ELC_Ben,ID=X)
ELC_Ctl<-read.table(paste0(con,'/significantDifferentialTest_Early_LCvsControl.csv'),sep=',',header = T)
ELC_Ctl<-rename(ELC_Ctl,ID=X)
ELC_Nor<-read.table(paste0(con,'/significantDifferentialTest_Early_LCvsNormal.csv'),sep=',',header = T)
ELC_Nor<-rename(ELC_Nor,ID=X)
Se.Sp.dat<-read.table('D:/LC autoantibodies/LC_diagnosis autoantibodies/Loess/originalData/highdensity/Se.Sp.revised.txt',sep = '\t',header = T)

##column no of three groups
col.LC<-colnames(normdata2[,7:66]) #EarlyLUAD
col.Benign<-colnames(normdata2[,67:96])#BLD
col.Normal<-colnames(normdata2[,97:120])#NHC

#1.EarlyLUADC_BLD FC¡Ý1£¬loess normalized data£¬group means£¬categorical 0,1 data (mean+2SD as threshold) 
ELC_Ben%>%filter(P.Value<0.05)%>%filter(logFC>log2(1.0))%>%left_join(normdata2[,c('Name','ID',col.LC,col.Benign)],by='ID')->ELC_Ben.Loess.FC1.0
ELC_Ben.Loess.FC1.0%>%mutate(mean_EarlyLC=apply(ELC_Ben.Loess.FC1.0[,col.LC],1,mean))%>%
  mutate(mean_Benign=apply(ELC_Ben.Loess.FC1.0[,col.Benign],1,mean))%>%left_join(Se.Sp.dat,by='ID')%>%filter(Se>=0.1)->ELC_Ben.Loess.FC1.0
#2. EarlyLUADC_NHC
ELC_Nor%>%filter(P.Value<0.05)%>%filter(logFC>log2(1))%>%left_join(normdata2[,c('Name','ID',col.LC,col.Normal)],by='ID')->ELC_Nor.Loess.FC1.0
ELC_Nor.Loess.FC1.0%>%mutate(mean_EarlyLC=apply(ELC_Nor.Loess.FC1.0[,col.LC],1,mean))%>%
  mutate(mean_Normal=apply(ELC_Nor.Loess.FC1.0[,col.Normal],1,mean))%>%left_join(Se.Sp.dat,by='ID')%>%filter(Se>=0.1)->ELC_Nor.Loess.FC1.0
#3. EarlyLUADC_Control(BLD+NHC)
ELC_Ctl%>%filter(P.Value<0.05)%>%filter(logFC>log2(1.0))%>%left_join(normdata2[,c('Name','ID',col.LC,col.Benign,col.Normal)],by='ID')->ELC_Ctl.Loess.FC1.0
ELC_Ctl.Loess.FC1.0%>%mutate(mean_EarlyLC=apply(ELC_Ctl.Loess.FC1.0[,col.LC],1,mean))%>%
  mutate(mean_Ctl=apply(ELC_Ctl.Loess.FC1.0[,c(col.Benign,col.Normal)],1,mean))%>%left_join(Se.Sp.dat,by='ID')%>%filter(Se>=0.1)->ELC_Ctl.Loess.FC1.0


####output
output.con<-'D:/LC autoantibodies/LC_diagnosis autoantibodies/Loess selection'
write.table(ELC_Ben.Loess.FC1.0,paste0(output.con,'/EarlyLC_Benign/ELC_Ben.Loess.FC1.0.csv'),sep=',')
write.table(ELC_Nor.Loess.FC1.0,paste0(output.con,'/EarlyLC_Normal/ELC_Nor.Loess.FC1.0.csv'),sep=',')
write.table(ELC_Ctl.Loess.FC1.0,paste0(output.con,'/EarlyLC_Benign+Normal/ELC_Ctl.Loess.FC1.0.csv'),sep=',')



###heatmap outline for primary check#######################################
library(pheatmap)
#EarlyLUAD_BLDign FC=1.0
EarlyLUAD.BLD.clus<-EarlyLUAD_BLD.Loess.FC1.0[,-c(3:9,100:105)]
group1<-c(rep('Early_LC',60),rep('BLDign',30))
annotation_c<-data.frame(group1)
rownames(annotation_c)<-colnames(EarlyLUAD.BLD.clus[,-c(1:2)])
p1<-pheatmap(EarlyLUAD.BLD.clus[,3:92],scale = 'row',clustering_method = 'ward.D2',annotation_col = annotation_c,annotation_legend = T,frontsize_row=7,main='EarlyLUAD_BLD.FC1.0_Se1.0')

#p1.1<-pheatmap(EarlyLUAD.BLD.clus[3:92],scale = 'row',cluster_cols = F,clustering_method = 'ward.D',annotation_col = annotation_c,annotation_legend = T,frontsize_row=7)
##EarlyLUAD_BLD FC=1.3
EarlyLUAD.Nor.clus<-EarlyLUAD_Nor.Loess.FC1.0[,-c(3:9,94:99)]
group2<-c(rep('Early_LC',60),rep('BLD',24))
annotation_c2<-data.frame(group2)
rownames(annotation_c2)<-colnames(EarlyLUAD.Nor.clus[,-c(1:2)])
p2<-pheatmap(EarlyLUAD.Nor.clus[3:86],scale = 'row',clustering_method = 'ward.D2',annotation_col = annotation_c2,annotation_legend = T,frontsize_row=7,main = 'EarlyLUAD_BLD.FC1.0_Se1.0')
##EarlyLUAD_BLDign+BLD FC=1.0
EarlyLUAD.Ctl.clus<-EarlyLUAD_Ctl.Loess.FC1.0[,-c(3:9,124:129)]
group3<-c(rep('Early_LC',60),rep('BLDign+BLD',54))
annotation_c3<-data.frame(group3)
rownames(annotation_c3)<-colnames(EarlyLUAD.Ctl.clus[,-c(1:2)])
p3<-pheatmap(EarlyLUAD.Ctl.clus[3:116],scale = 'row',clustering_method = 'ward.D2',annotation_col = annotation_c3,annotation_legend = T,frontsize_row=7,main = 'EarlyLUAD_Ctrl.FC1.0_Se1.0')

##EarlyLUAD_BLDign+BLD FC=1.0
EarlyLUAD.Ctl.clus<-EarlyLUAD_Ctl.Loess.FC1.0[,-c(3:9,124:129)]
group3<-c(rep('Early_LC',60),rep('BLDign+BLD',54))
annotation_c3<-data.frame(group3)
rownames(annotation_c3)<-colnames(EarlyLUAD.Ctl.clus[,-c(1:2)])
p3<-pheatmap(EarlyLUAD.Ctl.clus[3:116],scale = 'row',clustering_method = 'ward.D2',annotation_col = annotation_c3,annotation_legend = T,frontsize_row=7,main = 'EarlyLUAD_Ctrl.FC1.0_Se1.0')

##BLDign vs BLD FC=1.0
BLD.Nor.clus<-EarlyLUAD_Ctl.Loess.FC1.0[,-c(3:9,124:129)]
group3<-c(rep('Early_LC',60),rep('BLDign+BLD',54))
annotation_c3<-data.frame(group3)
rownames(annotation_c3)<-colnames(EarlyLUAD.Ctl.clus[,-c(1:2)])
p3<-pheatmap(EarlyLUAD.Ctl.clus[3:116],scale = 'row',clustering_method = 'ward.D2',annotation_col = annotation_c3,annotation_legend = T,frontsize_row=7,main = 'EarlyLUAD_Ctrl.FC1.0_Se1.0')


