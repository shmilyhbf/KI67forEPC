rm(list = ls())

StepVar="Step03 "
SavePath<-paste0(getwd(),"/Rdata/")

################################################################
load(paste0(SavePath, "Step02 matchdata.Rdata"))
if(F){
  openxlsx::write.xlsx(matchdata,colNames = TRUE,rowNames = TRUE,
                       file =paste0("./Rdata/","Supplementary Table S2.xlsx"))
  openxlsx::write.xlsx(mydata,colNames = TRUE,rowNames = TRUE,
                       file =paste0("./Rdata/","Supplementary Table S1.xlsx"))
  
}

dput(colnames(mydata))
myVars <- c( "Tstage", "Nstage","Stage", "VascularInvasion", "PerineuralInvasion", "PostTreatment", 
             "Gender", "Age", "AgeGroup", "TumorLocation", "Differentiation","Coloring.strength", "Positive.percentage", 
             "GeneScore")
surVars<-c("Futime.OS", "Fustat.OS")
strataVar<-"MKI67"

contVars=c("Age","Coloring.strength", "Positive.percentage", 
           "GeneScore")
catVars <- myVars[!myVars%in%c(contVars,surVars)]

if (strataVar=="") {
  for (i in catVars) {
    if(class(mydata[,i])!="factor"){
      mydata[,i]<-factor(mydata[,i])
    }
  } 
} else {
  for (i in c(strataVar,catVars)) {
    if(class(mydata[,i])!="factor"){
      mydata[,i]<-factor(mydata[,i])
    }
  }
}
str( mydata)

testVars<-c(catVars,contVars)

##########################################################################
###########################Loading Fnction######################################
functionpath<-"F:/Bioinformatics/myfunction/GeneralStatistics/"
source(paste0(functionpath,"myStatisticalFunction.R")) 

ttest.list<-list()
for (i in strataVar) {
  ttest.list[[i]]<-myStatisticalFunction(mydata=mydata, testvar.name=testVars,
                                         Catvar=catVars,Rankvar="",TotalShow=1,
                                         groupvar.name=i,forcedMethod=1,
                                         tablename=paste0(SavePath,StepVar,i))
}

ttest.list<-list()
for (i in strataVar) {
  ttest.list[[i]]<-myStatisticalFunction(mydata=matchdata, testvar.name=testVars,
                                         Catvar=catVars,Rankvar="",TotalShow=1,
                                         groupvar.name=i,forcedMethod=1,
                                         tablename=paste0(SavePath,StepVar,"match",i))
}
#######################################################################
median(mydata$Futime.OS)
range(mydata$Futime.OS)
table(mydata$Fustat.OS)
#################################################
###**************survival analysis******
###########################################################
dput(colnames(mydata))
targetGene="MKI67"
varsU=c(targetGene,"Stage", "VascularInvasion", "PerineuralInvasion", "PostTreatment", 
        "Gender", "AgeGroup", "TumorLocation", "Differentiation")

SurSet1=c("Futime.OS", "Fustat.OS")
Uregdata.OS=mydata[,colnames(mydata)%in%c(SurSet1,varsU)]
######################################################################
functionpath<-"F:/Bioinformatics/myfunction/bioinformatics/"
source(paste0(functionpath,"Step07 UniGene.KaplanMeier.R")) 

KM.res.OS<-UniGene.KaplanMeier(mydata=Uregdata.OS,groupVar=targetGene,TimeType="months",
                               legendlabs=levels(Uregdata.OS[,targetGene]),
                               plot.filename =paste0(SavePath,StepVar,"01 survivalplot.OS.png"),
                               saveplot=1,#
                               title.curves="A",method="AutoGroup")
if(F){
  openxlsx::write.xlsx(KM.res.OS[["SurvivalProportion"]],colNames = TRUE,rowNames = TRUE,
                       file =paste0(SavePath,StepVar,"01 SurvivalProportion.OS at diff timepoints.xlsx"))
}
KM.curve.OS=KM.res.OS[["ggsurvplot"]]

Uregdata.PFS=matchdata[,colnames(mydata)%in%c(SurSet1,varsU)]

#dd=survdiff(Surv(Futime.OS,Fustat.OS) ~ MKI67,data=matchdata)


KM.res.PFS<-UniGene.KaplanMeier(mydata=Uregdata.PFS,groupVar=targetGene,TimeType="months",
                                legendlabs=levels(Uregdata.PFS[,targetGene]),
                                plot.filename =paste0(SavePath,StepVar,"02 survivalplot.os match.png"),
                                saveplot=1,#
                                title.curves="B",method="AutoGroup")
if(F){
  openxlsx::write.xlsx(KM.res.PFS[["SurvivalProportion"]],colNames = TRUE,rowNames = TRUE,
                       file =paste0(SavePath,StepVar,"02 SurvivalProportion.PFS at diff timepoints.xlsx"))
}
KM.curve.PFS=KM.res.PFS[["ggsurvplot"]]
class(KM.curve.PFS)

surPlot.list<-list(KM.curve.OS, KM.curve.PFS)
surPlot<-arrange_ggsurvplots(surPlot.list, print = TRUE,labels = c("A", "B"),  
                             ncol = 1, nrow = length(surPlot.list)) 
class(surPlot)

ggsave(filename =paste0(SavePath,StepVar,"Figure 4.png"),surPlot,width =5, 
       height = 7, dpi = 300, units = "in", device='png')
ggsave(filename =paste0(SavePath,StepVar,"Figure 4.tif"),surPlot,width =5, 
       height = 7, dpi = 300, units = "in", device='tiff')
ggsave(filename =paste0(SavePath,StepVar,"Figure 4.jpg"),surPlot,width =5, 
       height = 7, dpi = 300, units = "in", device='jpg')


#########################################################
functionpath<-"F:/Bioinformatics/myfunction/GeneralStatistics/"
source(paste0(functionpath,"UniReg.cox.R")) 

UniReg.Res.OS<-UniReg.cox(varsU=varsU,Uregdata=Uregdata.OS,Pvalue=0.05,
                          colnames.Result=c("Characteristics", "Levels", "Beta", "SE", "HR (95% CI for HR)", 
                                            "Statistics(Z value)", "P"))

if(F){
  openxlsx::write.xlsx(UniReg.Res.OS[["Result"]],colNames = TRUE,rowNames = TRUE,
                       file =paste0(SavePath,StepVar,"Table 4.xlsx"))
}

source(paste0(functionpath,"MultiReg.cox.R"))
MultiReg.Res.OS<-MultiReg.cox(UniReg.Res=UniReg.Res.OS,y.var=c("futime","fustat"),
                              colnames.Result=c("Characteristics", "Levels", "Beta", "SE", "HR (95% CI for HR)", 
                                                "Statistics(Z value)", "P"))

if(F){
  openxlsx::write.xlsx(MultiReg.Res.OS[["Result.SigU"]],colNames = TRUE,rowNames = TRUE,
                       file =paste0(SavePath,StepVar,"Table 5.xlsx"))
  
  openxlsx::write.xlsx(MultiReg.Res.OS[["Result.beta.fitStep"]],colNames = TRUE,rowNames = TRUE,
                       file =paste0(SavePath,StepVar,"Table 6.xlsx"))
}

####################################

source(paste0(functionpath,"myforestplot.R"))
forest.unicox.OS=myforestplot(inputData=UniReg.Res.OS[["Result"]],effectindex="HR")
forest.Multicox.OS=myforestplot(inputData=MultiReg.Res.OS[["Result.SigU"]],effectindex="HR")
forest.Multicox.fitStep.OS=myforestplot(inputData=MultiReg.Res.OS[["Result.beta.fitStep"]],effectindex="HR")
p3.OS<-plot_grid(forest.unicox.OS, forest.Multicox.OS, forest.Multicox.fitStep.OS,
                 labels = "AUTO",ncol = 1,
                 rel_heights=c(6,3,2))
ggsave(filename =paste0(SavePath,StepVar,"Figure 5.png"),p3.OS,width =10, 
       height = 15, dpi = 300, units = "in", device='png')
ggsave(filename =paste0(SavePath,StepVar,"Figure 5.tif"),p3.OS,width =10, 
       height = 15, dpi = 300, units = "in", device='tiff')
ggsave(filename =paste0(SavePath,StepVar,"Figure 5.jpg"),p3.OS,width =10, 
       height = 15, dpi = 300, units = "in", device='jpg')

