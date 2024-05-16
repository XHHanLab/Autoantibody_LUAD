library(ggplot2)
library(dplyr)
library(tidyr)
library(purrr)
dat<-read.csv('D:/LC autoantibodies/LC_diagnosis autoantibodies/ELISA/dat for control plot/CVdat for plot del.csv', header=T)

CVmean<-dat%>%group_by(AAbs)%>%na.omit()%>%summarise_if(is.numeric,mean)

mean(CVmean$CV)

###generate labels####
data_text<-data.frame(label=paste0('CV','=',round(CVmean$CV,3),'%'),
                      
                      AAbs=CVmean$AAbs,
                      
                      x=rep(200,22),
                      
                      y=rep(45,22))



repnum=do.call(c,map(table(dat$AAbs),~seq(1,.x,by=1)))

p=ggplot(dat,aes(repnum,CV))+
  geom_point(stat = 'identity',position = 'identity')+
  facet_wrap(~AAbs, scales = "free",nrow = 2)+ylim(c(0,50))+
  geom_hline(aes(yintercept=20),colour='red',linetype='dashed',size=1.0)+
  scale_x_continuous(labels=c("0","100","200","300","400","500"))+
  theme_bw() + theme(panel.grid.minor =element_blank())+
  labs(x="n(samples)") +
  geom_text(data=data_text,mapping=aes(x=x,y=y,label=label),size=5,nudge_x=0.1,nudge_y=0.1)

ggsave('D:/LC autoantibodies/LC_diagnosis autoantibodies/ELISA/dat for control plot/CV plot.pdf',plot = last_plot(),height = 5,width =17,units = 'in')