if(F){
  BiocManager::install("BioinformaticsFMRP/GEObiolinksGUI.data")
  BiocManager::install("BioinformaticsFMRP/GEObiolinks")
  #安装github包GEOmirror和AnnoProbe
  remotes::install_github("jmzeng1314/GEOmirror")
  remotes::install_github("jmzeng1314/AnnoProbe")
  #三个包同时加载
}
# install cran packages
cran_pkgs <- c('tidyr',"plyr","ggfortify","rtracklayer", 'factoextra',  'devtools',
               'tibble','dplyr','stringr','ggplot2','ggpubr', 'FactoMineR',
               'cowplot', 'patchwork','basetheme','paletteer','AnnoProbe',
               'ggthemes','VennDiagram','tinyarray') 
bioc_pkgs <- c('GEOquery',"pheatmap",'hgu133plus2.db','ggnewscale',"limma",
               "impute","GSEABase","GSVA","clusterProfiler","org.Hs.eg.db",
               "preprocessCore","enrichplot","ggplotify")

for (pkg in bioc_pkgs){
  if (! require(pkg,character.only=T) ) {
    if(!require("BiocManager")) {
      install.packages("BiocManager",update = F,ask = F)
    }
    BiocManager::install(pkg,ask = F,update = F)
  }
  require(pkg,character.only=T)
}


for (pkg in cran_pkgs){
  if (!require(pkg,character.only=T) ) {
    install.packages(pkg,ask = F,update = F)
  }
  require(pkg,character.only=T)
}

###########################creat folders################################
rm(list = ls())
StepVar="Step05 "

###########创建项目及目录
project<- "ESCA"

#############################################GEOID#################################
geoID<-""

#R script保存名字
paste0("Step02 GEOdata calling.",project)
savename<-paste0(project,geoID)
#################################01 dir finding###########################
functionpath<-"F:/Bioinformatics/myfunction/bioinformatics/"
source(paste0(functionpath,"FileCreate.function.R")) 

ReadPath<-FileCreate.function(Type="Path",StaMethod="TCGA",
                              project=project,geoID=geoID,
                              Filename="")
ReadPath
SavePath<-FileCreate.function(Type="Path",StaMethod="DDGEO",
                              project=project,geoID=geoID,
                              Filename="Rdata")

SavePath="F:/桌面文件/食管癌Ki67/RawdataandRscript/Rdata/"

######################################GEOdata calling######################

GEOfilename<-list.files(ReadPath)
Fileselect<-GEOfilename[grepl("exprSet.annotated.Rds", GEOfilename)]
#*********注意查看文件情况************
Fileselect

data.list=readRDS(paste0(ReadPath,Fileselect))

######################################GEOdata calling######################

exprSet.log2=data.list[["exprSet.annotated"]][["exprSet.log2"]]
group_list=data.list[["group_list"]]


TargetGene="MKI67"

assign(paste0("exprSet.",TargetGene),as.vector(t(exprSet.log2[rownames(exprSet.log2)%in%c("MKI67"),grepl("Tumor",group_list)])))


#########################
IDtransDir<-"F:/Bioinformatics/idTrans/Rdata/"
idTrans.biomaRt.GTF.TCGA<-readRDS(paste0(IDtransDir,"000 idTrans.biomaRt.GTF.TCGA.V44.Rds"))
unique(idTrans.biomaRt.GTF.TCGA$gene_biotype)
#lncRNA.list<-unique(idTrans.biomaRt.GTF.TCGA$symbol[idTrans.biomaRt.GTF.TCGA$gene_biotype=="lncRNA"])
mRNA.list<-unique(idTrans.biomaRt.GTF.TCGA$symbol[idTrans.biomaRt.GTF.TCGA$gene_biotype=="protein_coding"])

#group.select=group_list[grepl("Tumor",group_list)]

Gene.exp=exprSet.log2[rownames(exprSet.log2)%in%mRNA.list,grepl("Tumor",group_list)]

functionpath<-"F:/Bioinformatics/myfunction/bioinformatics/"
source(paste0(functionpath,"Gene.CorAnalysis.R"))

TargetGene2=rownames(Gene.exp)
Exp_GS=as.data.frame(t(exprSet.log2[rownames(exprSet.log2)%in%c(TargetGene,TargetGene2),]))
Toal1 <- Gene.CorAnalysis(TargetGene=TargetGene,
                         TargetGene2=TargetGene2,
                         Exp_GS=Exp_GS)

Toal2 <- Gene.CorAnalysis2(TargetGene=TargetGene,
                         TargetGene2=TargetGene2,
                         Exp_GS=Exp_GS)

Toal3 <- filter(Toal2,abs(cor)>0.4 & FDR<0.05)
Toal3<-Toal3[Toal3$ToGene!=TargetGene,]

Toal3<-Toal3[order(Toal3$cor),]
Cor.GeneList<-Toal3$ToGene
Gene.neg<-as.character(head(Toal3$ToGene,20))
Gene.pos<-as.character(tail(Toal3$ToGene,20))
Gene.select=c(TargetGene,Gene.pos,Gene.neg)

####################################################################
library(corrplot)
library(rstatix) 
Exp_cor=as.data.frame(t(exprSet.log2[rownames(exprSet.log2)%in%Gene.select,]))
resmcor <- cor(Exp_cor)

resmorp <- cor.mtest(resmcor, conf.level = .95) 
p.mat <- resmorp$p

dfLong.resmcor = resmcor %>%as.data.frame(.) %>%
  rownames_to_column("FromGene") %>%  #("Gene")
  pivot_longer(-1,names_to = "ToGene",values_to = "CorrValue")
dfLong.resmcor

dfLong.p.mat = p.mat %>%as.data.frame(.) %>%
  rownames_to_column("Cell Type 1") %>%  #("Gene")
  pivot_longer(-1,names_to = "Cell Type 2",values_to = "PValue")
dfLong.p.mat
dfLong.resmcor$PValue=dfLong.p.mat$PValue

dfLong.resmcor$stars=ifelse(dfLong.resmcor$PValue<0.001,"***",
                            ifelse(dfLong.resmcor$PValue<0.01,"**",
                                   ifelse(dfLong.resmcor$PValue<0.05,"*","")))
level.cor=Gene.select
dfLong.resmcor$FromGene<-factor(dfLong.resmcor$FromGene,
                                     levels = level.cor)
dfLong.resmcor$ToGene<-factor(dfLong.resmcor$ToGene,
                                     levels =level.cor)

library(ggplot2)

corplot<-ggplot(dfLong.resmcor, aes(FromGene,ToGene)) +
  geom_tile(aes(fill=CorrValue)) +
  #geom_text(aes(label=stars), color="black", size=4) + 
  scale_fill_gradient2(low='blue', high='red',mid = 'white', 
                       #name=paste0("*    p < 0.05","\n\n","**  p < 0.01","\n\n","*** p < 0.001","\n\n","Correlation") ,
                       limit=c(-1,1)) + 
  labs(x=NULL,y=NULL) + 
  theme(axis.text.x = element_text(size=8,angle = 30,hjust = 1,color = "black"),
        axis.text.y = element_text(size=8,color = "black"),
        axis.ticks.y = element_blank(),
        panel.background=element_blank()) # 做一些简单的美化

corplot

if(F){
  ggsave(filename =paste0(SavePath,StepVar,"Figure 7.png"),corplot,width =10, 
         height = 10, dpi = 300, units = "in", device='png')
  
  ggsave(filename =paste0(SavePath,StepVar,"Figure 7.tif"),corplot,width =10, 
         height = 10, dpi = 300, units = "in", device='tiff')
  
  ggsave(filename =paste0(SavePath,StepVar,"Figure 7.jpg"),corplot,width =10, 
         height = 10, dpi = 300, units = "in", device='jpg')
  
  row.names(Toal3)=NULL
 
  openxlsx::write.xlsx(Toal3,colNames = TRUE,rowNames = TRUE,
                       file =paste0(SavePath,StepVar,"Supplementary Table S3.xlsx"))
}


########################################
###########################Loading function######################################
functionpath<-"F:/Bioinformatics/myfunction/bioinformatics/"
source(paste0(functionpath,"Step04 DEA.function.R")) 

###############################################################
group_list<-ifelse(get(paste0("exprSet.",TargetGene))<=median(get(paste0("exprSet.",TargetGene))),"Low","High")#data.list[["group_list"]]
group_list=factor(group_list,levels = c("Low","High"))
group.level<-levels(group_list)
expr<-Gene.exp #data.list[["exprSet.annotated"]][["exprSet.log2"]]


#trt="non_Responder";ctl="Responder"

if(length(levels(group_list))>2){
  print("group_list levels more than 2")
  print("Please set trt and ctl!")
}else{
  trt=levels(group_list)[2];ctl=levels(group_list)[1]
  Group<-group_list
}


if(length(group.level)>2 & trt!=""){
  Group<-group_list[grepl(paste0(trt,"|",ctl),group_list)]
  Group<-as.character(Group)
  Group<-factor(Group,levels = c(ctl,trt))
  expr<-expr[,grepl(paste0(trt,"|",ctl),group_list)]
}

########***查看group情况**** 
Group

###################################################
IDtransDir<-"F:/Bioinformatics/idTrans/Rdata/"
idTrans.biomaRt.GTF.TCGA<-readRDS(paste0(IDtransDir,"000 idTrans.biomaRt.GTF.TCGA.V44.Rds"))
unique(idTrans.biomaRt.GTF.TCGA$gene_biotype)
#lncRNA.list<-unique(idTrans.biomaRt.GTF.TCGA$symbol[idTrans.biomaRt.GTF.TCGA$gene_biotype=="lncRNA"])
mRNA.list<-unique(idTrans.biomaRt.GTF.TCGA$symbol[idTrans.biomaRt.GTF.TCGA$gene_biotype=="protein_coding"])

#expr.lncRNA<-expr[rownames(expr)%in%lncRNA.list,]
expr.mRNA<-expr[rownames(expr)%in%mRNA.list,]

###################################################
if(F){
exprSet.DEA.lncRNA<-Exp.filtered_genes(expr=expr.lncRNA,group_list=group_list, 
                                       ctl.trt=c(1,2),Filter.Method="Mean",DataType="Array")


DEA.limma.lncRNA<-DEA.function.limma(exprSet.DEA=exprSet.DEA.lncRNA,DEAindex="FDR",##"PValue"
                                     LogFC.method="Fix1",Colname="EN",voomaM=1)
}
###################################################

exprSet.DEA.mRNA<-Exp.filtered_genes(expr=expr.mRNA,group_list=group_list, 
                                     ctl.trt=c(1,2),Filter.Method="Mean",DataType="Array")

DEA.limma.mRNA<-DEA.function.limma(exprSet.DEA=exprSet.DEA.mRNA,DEAindex="FDR",##"PValue"
                                   LogFC.method="Fix1",Colname="EN",voomaM=1)

DEA.list.mRNA<-list()
DEA.list.mRNA[["group_list"]]<-group_list
DEA.list.mRNA[["exprSet.DEA.mRNA"]]<-exprSet.DEA.mRNA
DEA.list.mRNA[["DEA.limma.mRNA"]]<-DEA.limma.mRNA


##########################################################################
table(DEA.limma.mRNA$diffsig)

DEA.GeneList<-DEA.limma.mRNA$symbol[grepl("Down|Up",DEA.limma.mRNA$diffsig)]

DEA.DF=DEA.limma.mRNA[grepl("Down|Up",DEA.limma.mRNA$diffsig),]
DEA.DF=DEA.DF[order(DEA.DF$logFC),]

DEA.list.mRNA[["DEA.GeneList"]]<-DEA.GeneList
DEA.list.mRNA[["DEA.DF"]]<-DEA.DF


if(F){
  openxlsx::write.xlsx(DEA.DF,colNames = TRUE,rowNames = TRUE,
                       file =paste0(SavePath,StepVar,"Supplementary Table S4.xlsx"))
  saveRDS(DEA.list.mRNA,file=paste0(SavePath,StepVar,savename," DEA.list.mRNA.Rds"))
  
}

GeneList.select<-intersect(DEA.GeneList,Cor.GeneList)

DEA.DF2=DEA.limma.mRNA[GeneList.select,]
DEA.DF2=DEA.DF2[order(DEA.DF2$logFC),]
if(F){
  openxlsx::write.xlsx(DEA.DF2,colNames = TRUE,rowNames = TRUE,
                       file =paste0(SavePath,StepVar,"Table 8.xlsx"))
  
}


GeneList.select<-DEA.GeneList

###########################################################
table(group_list)

#######################################################################
###################*******volcano plot**********************
##***********************Data Prepairing***************
DEA.limma=DEA.limma.mRNA
###########################Loading Fnction######################################
if(F){
  functionpath<-"F:/Bioinformatics/myfunction/bioinformatics/"
  source(paste0(functionpath,"Step05 VolcanoHeatmapPlot.R")) 
  
  volcano.PlotData<-Volcano.function(DEA.limma=DEA.limma,volcanoMethod=1,
                                     Labelshow=1,top.gene="",
                                     DEAindex="FDR",GeneNo=10,
                                     width =5.5,height = 5.5,
                                     SavePath=SavePath,FigureType="png",
                                     FigureName=paste0(StepVar,"Figure 8"))
  
  volcano.PlotData
  Heatmap.PlotData<-Heatmap.function(exprSet=expr.mRNA,group_list=Group,
                                     DEAindex="FDR",GeneNo=10,DEA.limma=DEA.limma,
                                     top.gene="",
                                     pheatmapMethod=3,#1,2
                                     width =10,height = 5.5,SavePath=SavePath,
                                     FigureType="png",
                                     FigureName=paste0(StepVar,"Figure 8"))
  
  Heatmap.PlotData
  pMerge<-pheatmapMergevolcano(pheatmap.grid=Heatmap.PlotData,
                               volcano= volcano.PlotData,
                               width =10,height = 5.5,SavePath=SavePath,
                               FigureType="png",FigureName=paste0(StepVar,"Figure 8"))
  pMerge<-pheatmapMergevolcano(pheatmap.grid=Heatmap.PlotData,
                               volcano= volcano.PlotData,
                               width =10,height = 5.5,SavePath=SavePath,
                               FigureType="tiff",FigureName=paste0(StepVar,"Figure 8"))
  pMerge<-pheatmapMergevolcano(pheatmap.grid=Heatmap.PlotData,
                               volcano= volcano.PlotData,
                               width =10,height = 5.5,SavePath=SavePath,
                               FigureType="jpg",FigureName=paste0(StepVar,"Figure 8"))
  ##############################Saving data############################################
  save(DEA.limma,volcano.PlotData,
       Heatmap.PlotData,Group,
       pheatmapMergevolcano,
       file=paste0(SavePath,StepVar,"volcano.Heatmap.PlotData.Rdata"))
  
  load(paste0(SavePath,StepVar,"volcano.Heatmap.PlotData.Rdata"))
}



