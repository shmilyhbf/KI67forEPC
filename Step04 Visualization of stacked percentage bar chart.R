rm(list = ls())

StepVar="Step04 "
SavePath<-paste0(getwd(),"/Rdata/")

################################################################
load(paste0(SavePath, "Step02 matchdata.Rdata"))

###########################################
data.Tstage=as.data.frame(table(mydata$Tstage, mydata$MKI67Status))

ka.Tstage<-xtabs(~mydata$Tstage+mydata$MKI67Status,data=mydata)
ka.Tstage
chi.res<-chisq.test(ka.Tstage)
chi.res[["statistic"]]
Pvalue=ifelse(chi.res[["p.value"]]<0.001,"<0.001",
              paste0("=",sprintf("%.3f",round(chi.res[["p.value"]],3))))
text.x=paste0("Chi-squared Value= ",
              sprintf("%.3f",round(chi.res[["statistic"]],3)),
              ", P value",
              Pvalue)

library(scales)
library(ggplot2)
library(RColorBrewer)
library(dplyr)
library(plyr)

################********Calculate percentage**********
colnames(data.Tstage)
labels=levels(mydata$MKI67Status)
data3 = ddply(data.Tstage,'Var1',transform,percentage=Freq/sum(Freq),
              percent=paste0(sprintf("%.2f",round(Freq/sum(Freq)*100,2)),"%"))
data3=data3[data3$percentage!=0,]
p1.Tstage=ggplot(data3,aes(x=Var1,y=percentage,fill=Var2))+
  geom_bar(stat ="identity",alpha=0.9)+
  geom_text(aes(label=percent),######Count
            position = position_stack(vjust =0.5),size=3)+
  scale_fill_brewer(palette = "Set2",name = "Ki-67 status",
                    labels=labels)+
  #scale_fill_discrete(name = "Dose", labels = c("A", "B", "C"))+ 
  scale_y_continuous(expand = c(0,0),label=scales::percent)+
  labs(y="Percentage")+
  labs(x="Different T stage")+
  theme_bw()+
  theme(#axis.title.x =element_blank(),
    legend.key.size = unit(0.5,'cm'))
p1.Tstage=p1.Tstage+theme(plot.margin = margin(t = 1, r = 1, b = 1, l = 1, unit = "cm"))
p1.Tstage=p1.Tstage+ labs(title = text.x) 
p1.Tstage

##########################################
data.Nstage=as.data.frame(table(mydata$Nstage, mydata$MKI67Status))
ka.Nstage<-xtabs(~mydata$Nstage+mydata$MKI67Status,data=mydata)
ka.Nstage
chi.res<-chisq.test(ka.Nstage)
chi.res[["statistic"]]
Pvalue=ifelse(chi.res[["p.value"]]<0.001,"<0.001",
              paste0("=",sprintf("%.3f",round(chi.res[["p.value"]],3))))
text.x=paste0("Chi-squared Value= ",
              sprintf("%.3f",round(chi.res[["statistic"]],3)),
              ", P value",
              Pvalue)

colnames(data.Nstage)
labels=levels(mydata$MKI67Status)
data3 = ddply(data.Nstage,'Var1',transform,percentage=Freq/sum(Freq),
              percent=paste0(sprintf("%.2f",round(Freq/sum(Freq)*100,2)),"%"))
data3=data3[data3$percentage!=0,]
p1.Nstage=ggplot(data3,aes(x=Var1,y=percentage,fill=Var2))+
  geom_bar(stat ="identity",alpha=0.9)+
  geom_text(aes(label=percent),######Count
            position = position_stack(vjust =0.5),size=3)+
  scale_fill_brewer(palette = "Set2",name = "Ki-67 status",
                    labels=labels)+#+ scale_fill_discrete(name = "Dose", labels = c("A", "B", "C"))
  scale_y_continuous(expand = c(0,0),label=scales::percent)+
  labs(y="Percentage")+
  labs(x="Different lymph node metastasis status")+
  theme_bw()+
  theme(#axis.title.x =element_blank(),
    legend.key.size = unit(0.5,'cm'))
p1.Nstage=p1.Nstage+theme(plot.margin = margin(t = 1, r = 1, b = 1, l = 1, unit = "cm"))
p1.Nstage=p1.Nstage+ labs(title = text.x) 
p1.Nstage
###########################################
data.Stage=as.data.frame(table(mydata$Stage, mydata$MKI67Status))

ka.Stage<-xtabs(~mydata$Stage+mydata$MKI67Status,data=mydata)
ka.Stage
chi.res<-chisq.test(ka.Stage)
chi.res[["statistic"]]
Pvalue=ifelse(chi.res[["p.value"]]<0.001,"<0.001",
              paste0("=",sprintf("%.3f",round(chi.res[["p.value"]],3))))
text.x=paste0("Chi-squared Value= ",
              sprintf("%.3f",round(chi.res[["statistic"]],3)),
              ", P value",
              Pvalue)

colnames(data.Stage)
labels=levels(mydata$MKI67Status)
data3 = ddply(data.Stage,'Var1',transform,percentage=Freq/sum(Freq),
              percent=paste0(sprintf("%.2f",round(Freq/sum(Freq)*100,2)),"%"))
data3=data3[data3$percentage!=0,]
p1.Stage=ggplot(data3,aes(x=Var1,y=percentage,fill=Var2))+
  geom_bar(stat ="identity",alpha=0.9)+
  geom_text(aes(label=percent),######Count
            position = position_stack(vjust =0.5),size=3)+
  scale_fill_brewer(palette = "Set2",name = "Ki-67 status",
                    labels=labels)+
  scale_y_continuous(expand = c(0,0),label=scales::percent)+
  labs(y="Percentage")+
  labs(x="Different TNM stage")+
  theme_bw()+
  theme(#axis.title.x =element_blank(),
    legend.key.size = unit(0.5,'cm'))
p1.Stage=p1.Stage+theme(plot.margin = margin(t = 1, r = 1, b = 1, l = 1, unit = "cm"))
p1.Stage=p1.Stage+ labs(title = text.x) 
p1.Stage

############

###########################################
data.VascularInvasion=as.data.frame(table(mydata$VascularInvasion, mydata$MKI67Status))

ka.VascularInvasion<-xtabs(~mydata$VascularInvasion+mydata$MKI67Status,data=mydata)
ka.VascularInvasion
chi.res<-chisq.test(ka.VascularInvasion)
chi.res[["statistic"]]
Pvalue=ifelse(chi.res[["p.value"]]<0.001,"<0.001",
              paste0("=",sprintf("%.3f",round(chi.res[["p.value"]],3))))
text.x=paste0("Chi-squared Value= ",
              sprintf("%.3f",round(chi.res[["statistic"]],3)),
              ", P value",
              Pvalue)

colnames(data.VascularInvasion)
labels=levels(mydata$MKI67Status)
data3 = ddply(data.VascularInvasion,'Var1',transform,percentage=Freq/sum(Freq),
              percent=paste0(sprintf("%.2f",round(Freq/sum(Freq)*100,2)),"%"))
data3=data3[data3$percentage!=0,]
p1.VascularInvasion=ggplot(data3,aes(x=Var1,y=percentage,fill=Var2))+
  geom_bar(stat ="identity",alpha=0.9)+
  geom_text(aes(label=percent),######Count
            position = position_stack(vjust =0.5),size=3)+
  scale_fill_brewer(palette = "Set2",name = "Ki-67 status",
                    labels=labels)+
  scale_y_continuous(expand = c(0,0),label=scales::percent)+
  labs(y="Percentage")+
  labs(x="Different vascular invasion status")+
  theme_bw()+
  theme(#axis.title.x =element_blank(),
    legend.key.size = unit(0.5,'cm'))
p1.VascularInvasion=p1.VascularInvasion+theme(plot.margin = margin(t = 1, r = 1, b = 1, l = 1, unit = "cm"))
p1.VascularInvasion=p1.VascularInvasion+ labs(title = text.x) 
p1.VascularInvasion

###########################################
data.PerineuralInvasion=as.data.frame(table(mydata$PerineuralInvasion, mydata$MKI67Status))

ka.PerineuralInvasion<-xtabs(~mydata$PerineuralInvasion+mydata$MKI67Status,data=mydata)
ka.PerineuralInvasion
chi.res<-chisq.test(ka.PerineuralInvasion)
chi.res[["statistic"]]
Pvalue=ifelse(chi.res[["p.value"]]<0.001,"<0.001",
              paste0("=",sprintf("%.3f",round(chi.res[["p.value"]],3))))
text.x=paste0("Chi-squared Value= ",
              sprintf("%.3f",round(chi.res[["statistic"]],3)),
              ", P value",
              Pvalue)


colnames(data.PerineuralInvasion)
labels=levels(mydata$MKI67Status)
data3 = ddply(data.PerineuralInvasion,'Var1',transform,percentage=Freq/sum(Freq),
              percent=paste0(sprintf("%.2f",round(Freq/sum(Freq)*100,2)),"%"))
data3=data3[data3$percentage!=0,]
p1.PerineuralInvasion=ggplot(data3,aes(x=Var1,y=percentage,fill=Var2))+
  geom_bar(stat ="identity",alpha=0.9)+
  geom_text(aes(label=percent),######Count
            position = position_stack(vjust =0.5),size=3)+
  scale_fill_brewer(palette = "Set2",name = "Ki-67 status",
                    labels=labels)+#+ scale_fill_discrete(name = "Dose", labels = c("A", "B", "C"))
  scale_y_continuous(expand = c(0,0),label=scales::percent)+
  labs(y="Percentage")+
  labs(x="Different perineural invasion status")+
  theme_bw()+
  theme(#axis.title.x =element_blank(),
    legend.key.size = unit(0.5,'cm'))
p1.PerineuralInvasion=p1.PerineuralInvasion+theme(plot.margin = margin(t = 1, r = 1, b = 1, l = 1, unit = "cm"))
p1.PerineuralInvasion=p1.PerineuralInvasion+ labs(title = text.x) 
p1.PerineuralInvasion


library(patchwork)
p.merge=p1.Stage+p1.Nstage+  #p1.Tstage+
  p1.VascularInvasion+p1.PerineuralInvasion +
  plot_layout(guides = 'collect',nrow = 2)+plot_annotation(tag_levels = "A")

ggsave(filename =paste0(SavePath,StepVar,"Figure 4.pdf"),p.merge,
       width =12,height = 6, dpi = 300, units = "in", device='pdf')
ggsave(filename =paste0(SavePath,StepVar,"Figure 4.tif"),p.merge,
       width =12,height = 6, dpi = 300, units = "in", device='tiff')
ggsave(filename =paste0(SavePath,StepVar,"Figure 4.png"),p.merge,
       width =12,height = 6, dpi = 300, units = "in", device='png')
ggsave(filename =paste0(SavePath,StepVar,"Figure 4.jpg"),p.merge,
       width =12,height = 6, dpi = 300, units = "in", device='jpg')

