rm(list = ls())
load(paste0("./Rdata/","Step01 Raw survival data.Rdata"))
#save(mydata,file = paste0("./Rdata/","Step01 Raw survival data.Rdata"))
############################# psm###############################
package.need<-c('tableone', 'Matching','MatchIt',
                'survey','reshape2','ggplot2',"dplyr",
                "plyr","survival","survminer")

for (package in package.need) {
  if (!require(package, character.only=T, quietly=T)) {
    install.packages(package)
    suppressMessages(library(package, character.only=T))
  }
}
dput(names(mydata))

vars <- c("Tstage", "Nstage",# "Mstage", 
          "Stage", "MKI67", "VascularInvasion", "PerineuralInvasion", "PostTreatment", 
          "Gender", "Age", "AgeGroup", "TumorLocation", "Differentiation")

#vars <- c("sex","age","phone","city","freq")
tab <- CreateTableOne(vars=vars, strata = "MKI67" , data = mydata)
print(tab, showAllLevels = TRUE, smd = TRUE)

# psm
formula<-"MKI67~"
for (var in vars) {
  formula<-paste0(formula,var,"+")
}
formula
psm <- matchit(MKI67~Stage+PostTreatment+Gender+AgeGroup+TumorLocation+Differentiation,
               data=mydata,method="nearest",ratio = 1,
                 distance = "logit",
                 replace = FALSE,
                 caliper = 0.05)

summary(psm)
matchdata <- match.data(psm)
save(matchdata,mydata,file = paste0(SavePath,StepVar,"matchdata.Rdata"))
