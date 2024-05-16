
library(dplyr)
library(purrr)
library(rstatix)
library(reshape2)
library(caret)
library(pROC)
library(smotefamily)
library(adabag)
library(ROCR)
library(rpart)
library(randomForest)
library(mlr)
library(e1071)
library(smotefamily)
library(tidyr)
library(stringr)
library(tidyverse)
library(ggpubr)
#####1.read data####
setwd('./LC_diagnosis autoantibodies/ELISA')
ELISA<-read.csv('./combineddata/ELISA combination widedat include Control11.csv')
ELISA<-ELISA%>%select(-c('UCHL1_EDIT','UCHL1.IgG.re'))%>%filter(SampleID!='LC_B43')
ELISA$Group<-factor(ELISA$Group,levels = c('EarlyLUAD','BLD','NHC','Control'))
LUAD_NHC<-ELISA[ELISA$Group%in%c("EarlyLUAD","NHC"),]%>%
  select(-c(1,3))
LUAD_BLD<-ELISA[ELISA$Group%in%c("EarlyLUAD","BLD"),]%>%
  select(-c(1,3))

LUAD_Control<-ELISA[ELISA$Group%in%c("EarlyLUAD","Control"),]%>%
  select(-c(1:3))
LUAD_NHC<-ELISA[ELISA$Group%in%c("EarlyLUAD","NHC"),]%>%
  filter(!SampleID%in%c('M14','M39','NHC06','NHC20','NHC41'))%>%
  select(-c(1:3))


#####2.batch test####
#define batch test function 
ttest=function(dat,group){
  dat1=melt(dat[,1:length(colnames(dat))])%>%group_by(variable) %>%
    t_test(as.formula(paste0('value ~', group)))%>%as.data.frame
}



#######2.2 split data into ten-folds#####
dataOvS<-function(dat){
  dat$y<-factor(ifelse(dat$Group=="EarlyLUAD",1,0))
  dat<-dat[,c(ncol(dat),1,4:ncol(dat)-1)]
  set.seed(2022)
  dat1<-SMOTE(dat[,-c(1:2)],dat[,1],K=5) 
  print(table(dat1$data$class))
  ID.small=dat$SampleID[match(as.numeric(rownames(dat1[["orig_P"]])),rownames(dat))]
  SampleID=c(ID.small,paste0('syn',substr(ID.small[1],1,1),1:nrow(dat1[["syn_data"]])),dat$SampleID[startsWith(dat$SampleID,'LC_')])
  data<-data.frame(y=factor(dat1$data$class,levels = c(0,1)),SampleID,dat1$data[-length(colnames(dat1$data))])
}

LUAD_NHCwilogi<-dataOvS(LUAD_NHC)
LUAD_BLDwilogi<-dataOvS(LUAD_BLD)


set.seed(2022)
kfolds<-createFolds(LUAD_NHCwilogi$y,k=10)
kfolds1<-createFolds(LUAD_BLDwilogi$y,k=10)


#######2.3 ten-folds ttest####
funttest<-function(dat){
  vart<-list()
  for (i in 1:10) {
    test=dat[kfolds[[i]],]
    train=dat[-kfolds[[i]],]
    var=as.data.frame(ttest(train,'y'))
    varfinal=as.character(var[var$p<0.05,1])
    vart[[i]]=varfinal
  }
  vart<-Reduce(intersect,vart)
}
varttLUAD_BLD<-funttest(LUAD_BLDwilogi)
varttLUAD_NHC<-funttest(LUAD_NHCwilogi)
unionvartt<-union(varttLUAD_NHC,varttLUAD_BLD)

#######2.4. ten-folds ttest boxplot####

LUAD_NHCwilogi1=LUAD_NHCwilogi
LUAD_NHCwilogi1$y<-factor(LUAD_NHCwilogi1$y,levels = c(1,0),labels = c('EarlyLUAD','NHC'))
names(LUAD_NHCwilogi1)[1]='Group'
LUAD_BLDwilogi1=LUAD_BLDwilogi
LUAD_BLDwilogi1$y<-factor(LUAD_BLDwilogi1$y,levels = c(1,0),labels = c('EarlyLUAD','BLD'))
names(LUAD_BLDwilogi1)[1]='Group'


foldsdat=function(dat){
  foldsdat<-list()
  foldsdat0<-list()
  vart<-list()
  for (i in 1:10) {
    test=dat[kfolds[[i]],]
    train=dat[-kfolds[[i]],]
    var=as.data.frame(ttest(train,'Group'))
    varfinal=as.character(var[var$p<0.05,1])
    vart[[i]]=varfinal
    foldsdat0[[i]]=train
    foldsdat[[i]]=train[,c('Group',varfinal)]
  }
  vart1<-Reduce(intersect,vart)
  foldsdat1<-lapply(foldsdat,function(x)x[,c('Group',vart1)])##每个fold的差异交集
  foldsdat2<-lapply(foldsdat0,function(x)x[,c('Group',unionvartt)])#每个fold的最终两个比较的10个差异
  foldsdat3=list(vart1,vart,foldsdat,foldsdat1,foldsdat2)
}
viobox=function(dat){
  ggplot(dat, aes(AAbs, OD,color=Group))+ 
    geom_jitter(alpha=0.35,
                position=position_jitterdodge(jitter.width = 0.2, 
                                              jitter.height = 0, 
                                              dodge.width = 0.8))+
    geom_boxplot(alpha=0.35,width=0.2,
                 position=position_dodge(width=0.8),
                 size=0.6,outlier.colour = NA)+
    geom_violin(alpha=0.35,width=1.0,
                position=position_dodge(width=0.8),
                size=0.75)+
    scale_color_manual(values = unique(pal_nejm("default")(8)))+
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

foldsttLUAD_NHC<-foldsdat(LUAD_NHCwilogi1)
foldsttLUAD_BLD<-foldsdat(LUAD_BLDwilogi1)


########2.4.14 final boxplot####
vioboxdatLN<-lapply(foldsttLUAD_NHC[[5]],function(x)gather(x,key='AAbs',value = 'OD',-'Group'))
vioboxdatLB<-lapply(foldsttLUAD_BLD[[5]],function(x)gather(x,key='AAbs',value = 'OD',-'Group'))

files = str_c(str_c('./vioboxplot/vioboxplot of folds varfinal10 in LUADvsNHC in ELISA',1:10), ".pdf")
files1 = str_c(str_c('./vioboxplot/vioboxplot of folds varfinal10 in LUADvsBLD in ELISA',1:10), ".pdf")
walk2(files,lapply(vioboxdatLN,viobox),~ggsave(.x,.y,width=13,height = 6,units = 'in',dpi = 300))
walk2(files1,lapply(vioboxdatLB,viobox),~ggsave(.x,.y,width=10,height = 6,units = 'in',dpi = 300))


##Adaboost modeling ##########################################################################

#####3.1 define function####
adaboost<-function(dat,n,seedmax){
  seed_perform<-data.frame(seed=numeric(),auc=numeric())
  seed_perform1<-data.frame(seed=numeric(),auc=numeric())
  importance<-data.frame()
  probts<-list()
  probtr<-list()
  proball<-list()
  testdat=list()
  traindat=list()
  model<-list()
  for (i in seedmax) {
    set.seed(i)
    #split ELISA data into train and test
    #id=sample(1:nrow(dat), round(n*nrow(dat)))
    id=createDataPartition(dat$y,p=n,list=F)
    tr.dat=dat[-id,]
    ts.dat=dat[id,]
    
    ####Adaboost####
    #hyperparameter tuning
    m<- seq(100,500,1000)
    acc.dat=0
    
    for (j in 1:length(m))
    {
      boo.dat<-boosting(y~., tr.dat[-2], mfinal=m[j],v=10) #tenfolds
      dat.bt<-predict(boo.dat, tr.dat[-2])
      acc.dat[j]<- sum(as.numeric(dat.bt$class == tr.dat[,1]) )/ nrow(tr.dat)
    }
    #plot(m, acc.dat, xlab="mfinal", ylab="acuracy")
    best_j<- which.max(acc.dat)#optimal parameter
    #retrain the model 
    boo.bk<-boosting(y~., tr.dat[-2], mfinal=m[best_j])
    importancei<-boo.bk$importance 
    
    bt.bk<-predict(boo.bk, ts.dat[-2]) #test
    names(bt.bk)
    bt.bk$confusion
    sum(as.numeric(bt.bk$class == ts.dat[, 1]) )/ nrow(ts.dat)
    
    bt.bk1<-predict(boo.bk, tr.dat[-2]) #train
    names(bt.bk1)
    bt.bk1$confusion
    sum(as.numeric(bt.bk1$class == tr.dat[, 1]) )/ nrow(tr.dat)
    
    #ROC
    pred_boo <- prediction(bt.bk$prob[,2], ts.dat$y) 
    pred_boo1 <- prediction(bt.bk1$prob[,2], tr.dat$y) 

    auc_boo <- ROCR::performance(pred_boo, "auc")
    auc_boo1<-ROCR::performance(pred_boo1, "auc")
 
    auc_value_boo = round(auc_boo@y.values[[1]], 3)

    auc_value_boo1 = round(auc_boo1@y.values[[1]], 3)
    
    seed_perform=rbind(seed_perform,cbind(seed=i,auc=auc_value_boo))
    seed_perform1=rbind(seed_perform,cbind(seed=i,auc=auc_value_boo1))
    importance=rbind(importance,cbind(seed=i,t(data.frame(importancei))))
    print(c(which(seedmax==i),auc_value_boo))
    probtr[[which(seedmax==i)]]=data.frame(SampleID=tr.dat$SampleID,group=tr.dat$y,probtr=unlist(pred_boo1@predictions))##traindat
    probts[[which(seedmax==i)]]=data.frame(SampleID=ts.dat$SampleID,group=ts.dat$y,probts=unlist(pred_boo@predictions))##testdat
    proball[[which(seedmax==i)]]=data.frame(SampleID=c(tr.dat$SampleID,SampleID=ts.dat$SampleID),group=c(tr.dat$y,ts.dat$y),
                                            proball=c(unlist(pred_boo1@predictions),unlist(pred_boo@predictions)),
                                            class=c(rep('trainset',nrow(tr.dat)),rep('testset',nrow(ts.dat))))##train+testdat
    traindat[[which(seedmax==i)]]=tr.dat
    testdat[[which(seedmax==i)]]=ts.dat
    model[[which(seedmax==i)]]=boo.bk
  }
  names(probts)=paste0('seed',seedmax)
  names(probtr)=paste0('seed',seedmax)
  names(proball)=paste0('seed',seedmax)
  seed_perform=list(seed_perform,seed_perform1,probtr,probts,proball,traindat,testdat,model,importance)
}
set.seed(2022)
seedmax<-sample(1:10000,size=50)


###modeling ###
set.seed(2022)
#oversampling
dataOvS<-function(dat){
  dat$y<-factor(ifelse(dat$Group=="EarlyLUAD",1,0))
  dat<-dat[,c(ncol(dat),1,4:ncol(dat)-1)]
  set.seed(2022)
  dat1<-SMOTE(dat[,-c(1:2)],dat[,1],K=5) 
  print(table(dat1$data$class))
  ID.small=dat$SampleID[match(as.numeric(rownames(dat1[["orig_P"]])),rownames(dat))]
  SampleID=c(ID.small,paste0('syn',substr(ID.small[1],1,1),1:nrow(dat1[["syn_data"]])),dat$SampleID[startsWith(dat$SampleID,'LC_')])
  data<-data.frame(y=factor(dat1$data$class,levels = c(0,1)),SampleID,dat1$data[-length(colnames(dat1$data))])
}
LUAD_NHCmdl4.1<-dataOvS(LUAD_NHC)[,c('y','SampleID',unionvartt)] 
LUAD_BLDmdl4.1<-dataOvS(LUAD_BLD)[,c('y','SampleID',unionvartt)]

PfmLUAD_NHC4.1<-adaboost(dat=LUAD_NHCmdl4.1,n=0.25,seedmax = seedmax)
aucmeanLN4.1<-mean(PfmLUAD_NHC4.1[[1]]$auc)
PfmLUAD_BLD4.1<-adaboost(dat=LUAD_BLDmdl4.1,n=0.25,seedmax = seedmax)
aucmeanLB4.1<-mean(PfmLUAD_BLD4.1[[1]]$auc)

write.csv(LUAD_BLDmdl4.1,'./model/LUAD_BLD oversampledat in ELISA.csv',row.names = F)
write.csv(LUAD_NHCmdl4.1,'./model/LUAD_NHC oversampledat in ELISA.csv',row.names = F)

########3.3 plot ROC 50 curves#####
library(pROC)
library(ggsci)
colpalettes<-unique(c(pal_npg("nrc")(10),pal_aaas("default")(10),pal_nejm("default")(8),pal_lancet("lanonc")(9),
                      pal_jama("default")(7),pal_jco("default")(10)))
#########50 ROC#### 
roc0=roc(PfmLUAD_NHC4.1[[5]][[1]]$y,PfmLUAD_NHC4.1[[3]][[1]])
plot(smooth(roc0),legacy.axes=T,col=colpalettes[m],print.auc=F)
rocLN<-list()
for (m in 1:50) {
  rocLN[[m]]=roc(PfmLUAD_NHC4.1[[5]][[m]]$y,PfmLUAD_NHC4.1[[3]][[m]])
  roc=roc(PfmLUAD_NHC4.1[[5]][[m]]$y,PfmLUAD_NHC4.1[[3]][[m]])
  plot(smooth(roc),legacy.axes=F,col=colpalettes[m],add=TRUE,print.auc=F)
  legend("bottom", legend=paste0('roc',1:50),col=colpalettes[1:50],lty=1,cex=0.4,ncol = 5)
}

roc0.1=roc(PfmLUAD_BLD4.1[[5]][[1]]$y,PfmLUAD_BLD4.1[[3]][[1]])
plot(smooth(roc0.1),legacy.axes=T,col=colpalettes[m])
rocLN<-list()
for (m in 1:50) {
  rocLN[[m]]=roc(PfmLUAD_BLD4.1[[5]][[m]]$y,PfmLUAD_BLD4.1[[3]][[m]])
  roc=roc(PfmLUAD_BLD4.1[[5]][[m]]$y,PfmLUAD_BLD4.1[[3]][[m]])
  plot(smooth(roc),legacy.axes=F,col=colpalettes[m],add=TRUE)
  legend("bottomright", legend=paste0('roc',1:50),col=colpalettes[1:50],lty=1,cex=0.4,ncol = 5)
}



#### Autoantibodies combined with LDCT####
library(mice)
####6.1 read clinical information and imputation of missing value#####
clinProb<-read.csv('./model/ClinicInfo and probobility in LUAD_BLD.csv',header = T)

BLDmice<-clinProb%>%filter(ModelClass==0)%>%select(Age,PMD,IMD,LN)%>%as.matrix()%>%
  mice(m=5, maxit = 50, method = 'pmm', seed = 500)
BLDcomp<- complete(BLDmice,2)
summary(BLDcomp)


LUADmice<-clinProb%>%filter(ModelClass==1)%>%select(PMD,IMD,LN)%>%
 as.matrix()%>%mice(m=5, maxit = 50, method = 'pmm', seed = 500)
LUADcomp<- complete(LUADmice,2)
summary(LUADcomp)

####integrte data####
clinProball<-clinProb%>%mutate(Age.mice=c(clinProb$Age[1:234],BLDcomp$Age),
                               PMD.mice=c(LUADcomp$PMD,BLDcomp$PMD),
                               IMD.mice=c(LUADcomp$IMD,BLDcomp$IMD),
                               LN.mice=c(LUADcomp$LN,BLDcomp$LN),.before='mean.probts.')%>%
                       mutate(PMD.mice1=c(LUADcomp$PMD,BLDcomp1$PMD,rep('',100)),
                              IMD.mice1=c(LUADcomp$IMD,BLDcomp1$IMD,rep('',100)),
                              LN.mice1=c(LUADcomp$LN,BLDcomp1$LN,rep('',100)),.before='mean.probts.')


#####PPV calculation ####
#total PPV
CT.PPV=(234/(234+100))*100 
CT.PPV1=(234/(234+200))*100

#PPV for LDCT in nodules with 0.8cm,,0.9-20cm,>20mm

ppv0=function(x) x%>%summarise(CTPPV.nd=length(.$group[.$group=='EarlyLUAD'])/n())
CT.PPVnd<-clinProball%>%select(2,21,14:20)%>%
          mutate(IMDgroup=ifelse(IMD.mice<=0.8,1,ifelse(IMD.mice>2,3,2)))%>%
          split(.$IMDgroup)%>%map_dfr(ppv0,.id='IMDgroup')
CT.PPVndall<-clinProball%>%select(2,21,14:20)%>%
  mutate(IMDgroup=ifelse(IMD.mice<=0.8,1,ifelse(IMD.mice>2,3,2)))%>%ppv0 


####PPV for Autoantibodies-LDCT in nodules with 0.8cm,,0.9-20cm,>20mm
ppv=function(x)x%>%filter(PosNeg==1)%>%
  summarise(PPV.nd=length(.$group[.$group=='EarlyLUAD'])/n())

CT_AAbs.PPVnd<-clinProball%>%select(2,21,14:20)%>%mutate(PosNeg=ifelse(mean.probts.>=0.5,1,0))%>%
  mutate(IMDgroup=ifelse(IMD.mice<=0.8,1,ifelse(IMD.mice>2,3,2)))%>%
  split(.$IMDgroup)%>%map_dfr(ppv,.id='IMDgroup')###79%,71%,88%
CT_AAbs.PPVndall<-clinProball%>%select(2,21,14:20)%>%mutate(PosNeg=ifelse(mean.probts.>=0.5,1,0))%>%
  mutate(IMDgroup=ifelse(IMD.mice<=0.8,1,ifelse(IMD.mice>2,3,2)))%>%ppv #78.2%


####calculate RR for CRS####
####imputation data####
clinProbamc<-clinProb%>%mutate(Age.mice=c(clinProb$Age[1:234],BLDcomp$Age),
                               PMD.mice=c(LUADcomp$PMD,BLDcomp$PMD),
                               IMD.mice=c(LUADcomp$IMD,BLDcomp$IMD),
                               LN.mice=c(LUADcomp$LN,BLDcomp$LN),.before='mean.probts.')%>%
  mutate(Age.miceclass=case_when(Age.mice<=40~1,
                                 Age.mice>40&Age.mice<=50~2,
                                 Age.mice>50&Age.mice<=60~3,
                                 Age.mice>60~4),.after='Age')%>%
  mutate(Gender.mice=c(clinProb$Gender[1:234],ifelse(BLDcompgs$Gender01==1,'Male','Female')),
         Smoke.mice=c(clinProb$Smoke[1:234],ifelse(BLDcompgs$Smoke01==1,'Yes','No')),
         Alcohol.mice=c(clinProb$Alcohol[1:234],ifelse(BLDcompgs$Alcohol01==1,'Yes','No')),
         .before='GeneTest')

  write.csv(clinProbamc,'./model/PPV_SeSp in LUAD_BLD/clinic after interpolation and probability.csv',row.names = F)

#categorical variables
datlog<-read.csv('./model/PPV_SeSp in LUAD_BLD/clinic after interpolation and probability for logistic.csv',header = T)
datlog$Age.miceclass<-factor(datlog$Age.miceclass,levels=c(1,2,3,4))
datlog$Gender.mice<-factor(datlog$Gender.mice,levels=c(0,1))  
datlog$Smoke.mice<-factor(datlog$Smoke.mice,levels=c(1,2)) 
datlog$Alcohol.mice<-factor(datlog$Alcohol.mice,levels=c(1,2)) 
datlog$AAbsprobs<-factor(datlog$AAbsprobs,levels=c(0,1))
datlog$IMD.miceclass<-factor(datlog$IMD.miceclass,levels = c(1,2,3))

datlog$ModelClass<-factor(datlog$ModelClass,levels = c(0,1))
fit<-glm(ModelClass~Age.mice+Gender.mice+Smoke.mice+Alcohol.mice+IMD.miceclass+AAbsprobs,
         data=datlog,family=binomial())
coeffs=summary(fit)[["coefficients"]]

library(rio)
write.csv(coeffs,'./model/PPV_SeSp in LUAD_BLD/coeffs.csv',row.names = F)

###CRS was built in SAS with IMD,Gender and AAbs