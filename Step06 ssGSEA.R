#################################package downloading###################
bioc_pkgs <- c("ReactomePA","reactome.db","org.Hs.eg.db",
               "clusterProfiler","biomaRt","enrichplot",
               "DO.db","DOSE","pathview","ggstatsplot")
#BiocManager::install("ggstatsplot",ask = F,update = F)
for (pkg in bioc_pkgs){
  if (! require(pkg,character.only=T) ) {
    BiocManager::install(pkg,ask = F,update = F)
  }
  require(pkg,character.only=T)
}

package.list=c("tidyverse","data.table","ggridges",
               "dplyr","GOplot","ggplot2","forcats","cowplot",
               "ggnewscale","ggupset","ggstance","patchwork")
for (package in package.list) {
  if (!require(package,character.only=T, quietly=T)) {
    install.packages(package)
    library(package, character.only=T)
  }
}
rm(list = ls())
StepVar="Step06 "
SavePath="F:/桌面文件/食管癌Ki67/RawdataandRscript/Rdata/"


#load(paste0("F:/Bioinformatics/重要基因集/免疫基因集/immuGene.DF.Rdata"))
DEA.list.mRNA=readRDS("F:/桌面文件/食管癌Ki67/RawdataandRscript/Rdata/Step05 ESCA DEA.list.mRNA.Rds")
project="ESCA"
DEA.res=DEA.list.mRNA[["DEA.limma.mRNA"]]
#DEA.res=DEA.res[rownames(DEA.res)%in%GeneList.select,]
df <- data.frame(DEA.res$logFC,rownames(DEA.res))
colnames(df) <- c("logFC","SYMBOL")
head(df)#查看前面几行
dim(df)#数据总共几行几列

# 转换基因ID
df_id<-bitr(df$SYMBOL, 
            fromType = "SYMBOL",
            toType = "ENTREZID",
            OrgDb = "org.Hs.eg.db")#对应的物种，小鼠的是org.Mm.eg.db


df_all<-merge(df,df_id,by="SYMBOL",all=F)
head(df_all) 
dim(df_all) 

#df_geneselect<-df_all[df_all$SYMBOL%in%GeneList.select,]
df_geneselect<-df_all
#排序
df_all_sort <- df_geneselect[order(df_geneselect$logFC, decreasing = T),]#先按照logFC降序排序
gene_fc = df_all_sort$logFC #把foldchange按照从大到小提取出来
head(gene_fc)
names(gene_fc) <- df_all_sort$ENTREZID 
head(gene_fc)


#GO富集
GO_GSEA <- gseGO(pvalueCutoff = 1,
  gene_fc, #gene_fc
  ont = "ALL",# "BP"、"MF"和"CC"或"ALL"
  OrgDb = org.Hs.eg.db,#人类注释基因
  keyType = "ENTREZID")
p1 <- dotplot(GO_GSEA,showCategory = 5, title = "GO",
              split= "ONTOLOGY", font.size = 10)+
  facet_grid(ONTOLOGY~.,scales = "free")
p1
GO_GSEA.result=as.data.frame(GO_GSEA@result)

GO_GSEA.result<- GO_GSEA.result[order(GO_GSEA.result$ONTOLOGY,GO_GSEA.result$p.adjust, decreasing = F),]#先按照logFC降序排序

openxlsx::write.xlsx(GO_GSEA.result,colNames = TRUE,rowNames = TRUE,
                     file =paste0(SavePath,StepVar,project,"-GO_GSEA.result.xlsx"))
#ggsave(paste0(path_result,"GO.pdf"),p1,width = 6,height = 8)


#KEGG
KEGG_GSEA <- gseKEGG(gene_fc,pvalueCutoff = 1, organism = "hsa")
# 可视化
p2<- dotplot(KEGG_GSEA, showCategory=5,title=paste0("GSEA_KEGG")) 
p2
KEGG_GSEA.result=as.data.frame(KEGG_GSEA@result)

KEGG_GSEA.result<- KEGG_GSEA.result[order(KEGG_GSEA.result$p.adjust, decreasing = F),]#先按照logFC降序排序

openxlsx::write.xlsx(KEGG_GSEA.result,colNames = TRUE,rowNames = TRUE,
                     file =paste0(SavePath,StepVar,project," KEGG_GSEA.result.xlsx"))

# theother: groupGO, enrichGO, enrichKEGG, enrichMKEGG, enrichWP ,
# enricher ,enrichPathway, enrichDO, enrichNCG, enrichDGN , enrichMeSH

##wiki
Wiki_GSEA <- gseWP(gene_fc,pvalueCutoff = 1,organism = "Homo sapiens")
p3 <- dotplot(Wiki_GSEA, showCategory=5,title=paste0("GSEA_WikiPathway"))
p3
Wiki_GSEA.result=as.data.frame(Wiki_GSEA@result)

Wiki_GSEA.result<- Wiki_GSEA.result[order(Wiki_GSEA.result$p.adjust, decreasing = F),]#先按照logFC降序排序

openxlsx::write.xlsx(Wiki_GSEA.result,colNames = TRUE,rowNames = TRUE,
                     file =paste0(SavePath,StepVar,project," Wiki_GSEA.result.xlsx"))


##Reactome
reactome_GSEA <- gsePathway(gene_fc, pvalueCutoff = 1, verbose = FALSE)
p4 <- dotplot(reactome_GSEA, showCategory=5,title=paste0("GSEA_Reactome"))
p4
reactome_GSEA.result=as.data.frame(reactome_GSEA@result)

reactome_GSEA.result<- reactome_GSEA.result[order(reactome_GSEA.result$p.adjust, decreasing = F),]#先按照logFC降序排序

openxlsx::write.xlsx(reactome_GSEA.result,colNames = TRUE,rowNames = TRUE,
                     file =paste0(SavePath,StepVar,project," reactome_GSEA.result.xlsx"))


#结果放在一张图中
p5 <- p1/p2/p3/p4

class(p1)
pMerge <-p5 + plot_annotation(tag_levels = "A") +     ## 大写字母编号
              plot_layout(ncol = 1,#图形设置为两列，默认按行填充，
              #widths = c(2, 1),
              heights=c(2,1,1,1.5))#两列之间相对宽度比为2：1
ggsave(filename =paste0(SavePath,StepVar,project," Figure 9.tif"),pMerge,width =9, 
       height =20, dpi = 300, units = "in", device='tiff')
ggsave(filename =paste0(SavePath,StepVar,project," Figure 9.png"),pMerge,width =9, 
       height =20, dpi = 300, units = "in", device='png')
ggsave(filename =paste0(SavePath,StepVar,project," Figure 9.jpg"),pMerge,width =9, 
       height =20, dpi = 300, units = "in", device='jpg')

##保存数据
save(GO_GSEA,KEGG_GSEA,Wiki_GSEA,reactome_GSEA,
     file = paste0(SavePath,StepVar,project," enrichment_GSEA.Rdata"))
#"F:/Bioinformatics/TCGA突变数据分析/LUAD/Rdata/enrichment_GSEA.Rdata"


#gseaplot2用法 ##用misgdb得到的结果,直接按照序号出图，可输入"hsa04621"替换为3
View(KEGG_GSEA@result) #查看需要绘制的通路ID号
# 绘制"Nod-like receptor signaling pathway"即hsa04621

p6 <- gseaplot2(
  KEGG_GSEA, #gseaResult object，即GSEA结果
  "hsa01521",#富集的ID编号
  title = "EGFR tyrosine kinase inhibitor resistance", #标题
  color = "green",#GSEA线条颜色
  base_size = 11,#基础字体大小
  rel_heights = c(1.5, 0.5, 1),#副图的相对高度
  subplots = 1:3, #要显示哪些副图 如subplots=c(1,3) #只要第一和第三个图，subplots=1#只要第一个图
  pvalue_table = FALSE, #是否添加 pvalue table
  ES_geom = "line" #running enrichment score用先还是用点ES_geom = "dot"
)
p6
ggsave(paste0(SavePath,StepVar,project," GSEA_KEGG2.png"),p6,width = 6,height = 4)


################################################################
