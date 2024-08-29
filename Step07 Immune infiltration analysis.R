rm(list=ls())
StepVar="Step07 "

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
group_list=data.list[["group_list"]]

exprSet.TPM <-data.list[["exprSet.annotated"]][["exprSet.TPM"]]

#########################
IDtransDir<-"F:/Bioinformatics/idTrans/Rdata/"
idTrans.biomaRt.GTF.TCGA<-readRDS(paste0(IDtransDir,"000 idTrans.biomaRt.GTF.TCGA.V44.Rds"))
unique(idTrans.biomaRt.GTF.TCGA$gene_biotype)
#lncRNA.list<-unique(idTrans.biomaRt.GTF.TCGA$symbol[idTrans.biomaRt.GTF.TCGA$gene_biotype=="lncRNA"])
mRNA.list<-unique(idTrans.biomaRt.GTF.TCGA$symbol[idTrans.biomaRt.GTF.TCGA$gene_biotype=="protein_coding"])

exprSet.TPM <-exprSet.TPM[rownames(exprSet.TPM)%in%mRNA.list,
                          grepl("Tumor",group_list)]
exprSet.TPM <- log2(exprSet.TPM + 1)

TargetGene="MKI67"
assign(paste0("exprSet.",TargetGene),
       as.vector(t(exprSet.TPM[rownames(exprSet.TPM)%in%c("MKI67"),grepl("Tumor",group_list)])))
group_list<-ifelse(get(paste0("exprSet.",TargetGene))<=median(get(paste0("exprSet.",TargetGene))),"Low","High")#data.list[["group_list"]]
group_list=factor(group_list,levels = c("Low","High"))
group.level<-levels(group_list)
########regroup with group_list
grp1=group_list[grepl(group.level[1],group_list)]
grp2=group_list[grepl(group.level[2],group_list)]

exprSet.TPM1=exprSet.TPM[,grepl(group.level[1],group_list)]
exprSet.TPM2=exprSet.TPM[,grepl(group.level[2],group_list)]
group_list=c(grp1,grp2)
exprSet.TPM=cbind(exprSet.TPM1,exprSet.TPM2)


library(GSVA)

#https://doi.org/10.1016/j.celrep.2016.12.019)
geneset<-openxlsx::read.xlsx(paste0(SavePath,"mmc3.xlsx"),
                             colNames = TRUE,rowNames = F)
#*
geneset2 <- split(geneset$Metagene,geneset$Cell.type)

geneset2[1:2]

#*Immune cell infiltration analysis
res.gsva <- gsva(
  as.matrix(exprSet.TPM), 
  geneset2, 
  method = "ssgsea", #"gsva", "ssgsea", "zscore" and "plage"
  kcdf = "Gaussian",
  mx.diff = F,
  verbose   = F 
)

#Min-Max scaled
res.scaled <- res.gsva
for (i in colnames(res.gsva)) {
  res.scaled[,i] <- (res.gsva[,i] -min(res.gsva[,i]))/(max(res.gsva[,i] )-min(res.gsva[,i] ))
}

if(F){
  res.gsva.df=as.data.frame(res.gsva)
  openxlsx::write.xlsx(res.gsva.df,colNames = TRUE,rowNames = TRUE,
                       file =paste0(SavePath,StepVar,project," Supplementary Table S6.xlsx"))
}
gsvaList=list(exprSet.TPM,geneset2,res.gsva,res.scaled)
names(gsvaList)=c("exprSet.TPM","geneset2","res.gsva","res.scaled")

save(gsvaList,file = paste0(SavePath,StepVar,project,"-immune_cell_result.Rdata")) #保存一下结果，方便下次使用#*

#*

options(stringsAsFactors = F)
library(pacman)
#chooseCRANmirror()
#install.packages("ggvenn")
p_load(tidyverse,ggpubr,paletteer,corrplot,rstatix,ggheatmap,aplot,gridExtra,
       pheatmap,ggplot2,stringr,ggplotify)

#########################################################
##################************heatplot**************
#First Method
if(F){
  pheatmap(res.gsva, show_colnames = F)
  #group_list <- ifelse(str_sub(colnames(res.gsva),14,15) < 10,"tumor","normal")
  annotation <- data.frame(group_list)
  rownames(annotation) <- colnames(res.gsva)
  table(annotation)
  head(annotation)
  pheatmap(res.gsva,show_colnames = F,annotation_col = annotation,fontsize = 10)
}



########Second Method
ncol(res.gsva)
length(group_list)
length(colnames(res.gsva))


col_metaData <- data.frame(Group=group_list)
rownames(col_metaData) <- colnames(res.gsva)


Groupcol <- c("#98D352","#FF7F0E")
names(Groupcol) <- levels(group_list)
col <- list(Group=Groupcol)
text_columns <- sample(colnames(res.gsva),0)

#########wide to long
dfLong = res.gsva %>%as.data.frame(.) %>%
  rownames_to_column("Cell Type") %>%  #("Gene")
  pivot_longer(-1,names_to = "Sample",values_to = "Value")
dfLong

dfLong$Sample<-factor(dfLong$Sample,levels = unique(dfLong$Sample))


p = ggplot(dfLong,aes(x=Sample,y=`Cell Type`,fill=Value))+
  geom_raster()+
  scale_fill_gradientn(colours = c("blue", "green", "yellow", "red"))+
  #scale_fill_gradient2(low="#0000ff", high="#ff0000", mid="#ffffff")+
  scale_y_discrete(position="right") +
  theme_minimal()+
  theme(panel.grid.major=element_blank(),
        axis.text.x = element_blank())
p
#*
#*
pGroup = col_metaData %>%
  mutate(X = factor(colnames(res.gsva),levels=unique(colnames(res.gsva))),Y="Group") %>%
  ggplot(aes(x=X,y=Y,fill=Group))+
  geom_tile() +
  theme_void()+
  labs(fill = "Group")+
  theme(panel.grid.major=element_blank(),
        axis.text.x = element_blank())
pGroup

heatplot=p %>%  insert_top(pGroup, height = .05)

heatplot

#########################################################
##################************boxplot**************
dt <- res.scaled %>% t() %>% as.data.frame() %>%
  rownames_to_column("sample") %>%
  gather(key = cell_type,
         value = value, -sample)
head(dt)
#*
dtt <- dt %>%
  group_by(sample) %>%
  mutate(proportion = round(value/sum(value),3))
head(dtt)
#*
dtt$cell_type <- factor(dtt$cell_type,levels = unique(rownames(res.gsva)))

#Custom Theme:
mytheme <- theme(axis.title = element_text(size = 12),
                 axis.text.x = element_blank(),
                 axis.ticks.x = element_blank(),
                 plot.title = element_text(size = 13,
                                           hjust = 0.5,
                                           face = "bold"),
                 legend.text = element_text(size = 10),
                 legend.position = "bottom")

d_palettes <- palettes_d_names #(paletteer package)
col <- paletteer_d("Redmonder::dPBIPuOr", 28, type = "continuous")

p <- ggplot(dtt,
            aes(x = cell_type,y = proportion,fill = cell_type)) +
  geom_boxplot(color = "black",alpha = 0.6,outlier.shape = 21,outlier.size = 1.2) +
  scale_fill_manual(values = col) +
  labs(x = "cell type", y = "proportion") +
  theme_bw() + mytheme
p
#*
#*
#*#reodering
dtt_arrange <- dtt %>%
  group_by(cell_type) %>%
  summarise(de = median(proportion)) %>%
  arrange(desc(de)) %>%
  pull(cell_type)

dtt$cell_type <- factor(dtt$cell_type,levels = unique(dtt_arrange))

# replotting
p1 <- ggplot(dtt,
             aes(x = cell_type,y = proportion,fill = cell_type)) +
  geom_boxplot(color = "black",alpha = 0.6,outlier.shape = 21,outlier.size = 1.2) +
  scale_fill_manual(values = col) +
  labs(x = "cell type", y = "proportion") +
  theme_bw() + mytheme
p1

##################************boxplot**************
##**********3. Group box plot+significance markers***********
dtt$group <-rep(group_list,length.out=nrow(dtt))

p2 <- ggplot(dtt,
             aes(x = cell_type,y = proportion,fill = group)) +
  geom_boxplot(color = "black",alpha = 0.6,outlier.shape = 21,outlier.size = 1.2) +
  scale_fill_manual(values = c("#4979b6","#d9352a")) +
  labs(x = "cell type", y = "proportion") +
  theme_bw() + mytheme + theme(axis.text.x = element_text(angle = 45))
p2
#*t test or wilcox test to test
t <- t_test(group_by(dtt, cell_type), proportion ~ group)
tj <- adjust_pvalue(t, method = 'fdr') 
tj
#* sig index
tj <- add_significance(tj, 'p.adj')
tj
#class(tj)
#*
#*adding sig index
lab <- add_xy_position(tj, x = 'cell_type', dodge = 0.65)
lab$p.adj.signif=ifelse(lab$p.adj.signif=="****","***",lab$p.adj.signif)
if(F){
  lab.df=as.data.frame(lab)
  lab.df=lab.df[order(lab.df$p.adj.signif),]
  openxlsx::write.xlsx(lab.df,colNames = TRUE,rowNames = TRUE,
                       file =paste0(SavePath,StepVar,project," Table 9.xlsx"))
}

#ggpubr plotting
p3 <- ggboxplot(dtt, x = "cell_type", y = "proportion",
                fill = "group", color = "black",
                width=0.7, alpha = 0.6,
                outlier.shape = 21, outlier.size = 1.2) +
  scale_fill_manual(values = c("#4979b6","#d9352a")) +
  labs(x = "cell type", y = "proportion") +
  theme_bw() + mytheme + theme(axis.text.x = element_text(angle = 45)) +
  stat_pvalue_manual(lab, label = 'p.adj.signif', label.size=3, bracket.size=0.5, tip.length = 0.02)
p3

Boxplot=p3+theme(plot.margin = margin(t = 1, r = 2, b = 1, l = 1, unit = "cm"))
Boxplot

###################################################################
###############*********Immune cell related heatmap***************
#Calculate the correlation between immune cells:
res.scaled.cor <- cor(t(res.scaled))
#View(res.scaled.cor)
#Add significance markers;Use cor.mtest for significance testing
resmorp <- cor.mtest(res.scaled.cor, conf.level = .95) 

p.mat <- resmorp$p
################wide to long
dfLong.resmcor = res.scaled.cor %>%as.data.frame(.) %>%
  rownames_to_column("Cell Type 1") %>%  #("Gene")
  pivot_longer(-1,names_to = "Cell Type 2",values_to = "CorrValue")
dfLong.resmcor

dfLong.p.mat = p.mat %>%as.data.frame(.) %>%
  rownames_to_column("Cell Type 1") %>%  #("Gene")
  pivot_longer(-1,names_to = "Cell Type 2",values_to = "PValue")
dfLong.p.mat
dfLong.resmcor$PValue=dfLong.p.mat$PValue

dfLong.resmcor$stars=ifelse(dfLong.resmcor$PValue<0.001,"***",
                            ifelse(dfLong.resmcor$PValue<0.01,"**",
                                   ifelse(dfLong.resmcor$PValue<0.05,"*"," ")))

dfLong.resmcor$`Cell Type 1`<-factor(dfLong.resmcor$`Cell Type 1`,
                                     levels = unique(dfLong.resmcor$`Cell Type 1`))
dfLong.resmcor$`Cell Type 2`<-factor(dfLong.resmcor$`Cell Type 2`,
                                     levels = unique(dfLong.resmcor$`Cell Type 2`))

## Drawing with geom_tile in ggplot2
corplot<-ggplot(dfLong.resmcor, aes(`Cell Type 1`,`Cell Type 2`)) +
  geom_tile(aes(fill=CorrValue)) +
  geom_text(aes(label=stars), color="black", size=4) +
  scale_fill_gradient2(low='blue', high='red',mid = 'white', 
                       limit=c(-1,1),
                       name=paste0("*    p < 0.05","\n\n","**  p < 0.01","\n\n","*** p < 0.001","\n\n","Correlation")) +
  labs(x=NULL,y=NULL) +
  theme(axis.text.x = element_text(size=8,angle = 30,hjust = 1,color = "black"),
        axis.text.y = element_text(size=8,color = "black"),
        axis.ticks.y = element_blank(),
        panel.background=element_blank()) 

corplot
#######################################################################
#####*****Heat map of the correlation between immune cells and target genes****

identical(colnames(res.scaled),colnames(exprSet.TPM))

DEA.list.mRNA=readRDS("F:/桌面文件/食管癌Ki67/RawdataandRscript/Rdata/Step05 ESCA DEA.list.mRNA.Rds")
DEA.GeneList=DEA.list.mRNA[["DEA.GeneList"]]

load(paste0(SavePath,"immuGene.DF.Rdata"))

genes <-DEA.GeneList[DEA.GeneList%in%immuGene.list]
exp_genes <- exprSet.TPM[genes,]

rb <- rbind(res.scaled,exp_genes)

rownames(rb)

rbcor <- cor(t(rb))
rbcorp <- cor.mtest(rbcor, conf.level = .95)
p.mat2 <- rbcorp$p

split <- rbcor[1:nrow(res.scaled), 
               (nrow(rb)-length(genes)+1):nrow(rb)] 
#View(split)
#*
#*Cutting p-value matrix:
splitp <- p.mat2[1:nrow(res.scaled), 
                 (nrow(rb)-length(genes)+1):nrow(rb)]

dfLong.split = split %>%as.data.frame(.) %>%
  rownames_to_column("Cell Type") %>%  #("Gene")
  pivot_longer(-1,names_to = "Gene",values_to = "CorrValue")
dfLong.split

dfLong.splitp = splitp %>%as.data.frame(.) %>%
  rownames_to_column("Cell Type") %>%  #("Gene")
  pivot_longer(-1,names_to = "Gene",values_to = "PValue")
dfLong.splitp
dfLong.split$PValue=dfLong.splitp$PValue

#dfLong.split$stars=ifelse(dfLong.split$PValue<0.001,"***",
#                          ifelse(dfLong.split$PValue<0.01,"**",
#                                 ifelse(dfLong.split$PValue<0.05,"*"," ")))
#
dfLong.split$stars=ifelse(dfLong.split$PValue<0.05,"*"," ")

dfLong.split$`Cell Type`<-factor(dfLong.split$`Cell Type`,
                                 levels = unique(dfLong.split$`Cell Type`))
dfLong.split$`Gene`<-factor(dfLong.split$`Gene`,
                            levels = unique(dfLong.split$`Gene`))

#*plotting
corplot2<-ggplot(dfLong.split, aes(`Cell Type`,`Gene`)) +
  geom_tile(aes(fill=CorrValue)) +
  geom_text(aes(label=stars), color="black", size=4) +
  scale_fill_gradient2(low='blue', high='red',mid = 'white', 
                       name=paste0("*    P < 0.05","\n\n","Correlation"), 
                       #name=paste0("*    p < 0.05","\n\n","**  p < 0.01","\n\n","*** p < 0.001","\n\n","Correlation"), 
                       limit=c(-1,1)) + 
  labs(x=NULL,y=NULL) + 
  theme(axis.text.x = element_text(size=8,angle = 30,hjust = 1,color = "black"),
        axis.text.y = element_text(size=8,color = "black"),
        axis.ticks.y = element_blank(),
        panel.background=element_blank()) 
corplot2
#library(gridExtra)
#combined_plot <- grid.arrange(heatplot,Boxplot.LUAD,corplot,corplot2, 
#                              nrow = 4, ncol = 1)
pheatmap.grid <- ggplotify::grid2grob(print(heatplot))

pMerge<-cowplot::plot_grid(pheatmap.grid,Boxplot,corplot,corplot2,
                           labels = "AUTO",
                           #label_size= 12,
                           #label_fontfamily= "serif",
                           #label_fontface= "plain",
                           #label_colour= "black",
                           rel_widths = c(1.5, 1),
                           rel_heights=c(1.0,1),
                           ncol=2)
pMerge

if(F){
  ggsave(filename =paste0(SavePath,StepVar,"Figure 10.png"),pMerge,width =16, 
         height = 12, dpi = 300, units = "in", device='png')
  
  ggsave(filename =paste0(SavePath,StepVar,"Figure 10.tiff"),pMerge,width =16, 
         height = 12, dpi = 300, units = "in", device='tiff')
  
  ggsave(filename =paste0(SavePath,StepVar,"Figure 10.jpg"),pMerge,width =16, 
         height = 12, dpi = 300, units = "in", device='jpg')
}
