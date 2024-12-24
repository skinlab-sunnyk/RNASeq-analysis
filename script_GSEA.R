library(tidyverse)
library(limma)
library(gplots) #for heatmaps
library(DT) #interactive and searchable tables of our GSEA results
library(GSEABase) #functions and methods for Gene Set Enrichment Analysis
library(Biobase) #base functions for bioconductor; required by GSEABase
library(GSVA) #Gene Set Variation Analysis, a non-parametric and unsupervised method for estimating variation of gene set enrichment across samples.
library(gprofiler2) #tools for accessing the GO enrichment results using g:Profiler web resources
library(patchwork)
library(clusterProfiler) # provides a suite of tools for functional enrichment analysis
library(msigdbr) # access to msigdb collections directly within R
library(enrichplot) # great for making the standard GSEA enrichment plots
library(readr)

myData_q <- subset(myTopHits_0, myTopHits_0$adj.P.Val<0.05)
myData_q.ord.byFC <- myData_q[order(myData_q$logFC, decreasing = T),]
myData_q.ord.byFC.selectbyfc <- subset(myData_q.ord.byFC, abs(myData_q.ord.byFC$logFC)>log2(1.5))
dim(myData_q)
dim(myData_q.ord.byFC)
dim(myData_q.ord.byFC.selectbyfc)
mydata_UP<-subset(myData_q, myData_q$logFC>0)
dim(mydata_UP)
mydata_DOWN<-subset(myData_q, myData_q$logFC<0)
dim(mydata_DOWN)
mydata_UP_FC1.5<-subset(mydata_UP, mydata_UP$logFC>log2(1.5))
dim(mydata_UP_FC1.5)
mydata_DOWN_FC1.5<-subset(mydata_DOWN, mydata_DOWN$logFC < -log2(1.5))
dim(mydata_DOWN_FC1.5)


myData_q.ord.byFC.tibble <- as_tibble(myData_q.ord.byFC, rownames = "GeneID")
myData_q.ord.byFC.selectbyfc.tibble <- as_tibble(myData_q.ord.byFC.selectbyfc, rownames = "GeneID")

mydata.df.sub <- dplyr::select(myData_q.ord.byFC.tibble, GeneID, logFC)
mydata.gsea <- mydata.df.sub$logFC
names(mydata.gsea) <- as.character(mydata.df.sub$GeneID)

mydata.df.sub2 <- dplyr::select(myData_q.ord.byFC.selectbyfc.tibble, GeneID, logFC)
mydata.gsea2 <- mydata.df.sub2$logFC
names(mydata.gsea2) <- as.character(mydata.df.sub2$GeneID)

C2CP <- read.gmt("m5.go.v2023.2.Mm.symbols.gmt")
set.seed(123) #set a random seed so that we can reproducible ordering for our GSEA results below
myGSEA.res <- GSEA(mydata.gsea, TERM2GENE=C2CP, verbose=FALSE) #could replace C2CP with hs_gsea_c2 object you retrieved from msigdb above
beep(sound = 1)

myGSEA.df <- as_tibble(myGSEA.res@result)


set.seed(124) #set a random seed so that we can reproducible ordering for our GSEA results below
myGSEA.res2 <- GSEA(mydata.gsea2, TERM2GENE=C2CP, verbose=FALSE) #could replace C2CP with hs_gsea_c2 object you retrieved from msigdb above
beep(sound = 1)

myGSEA.df2 <- as_tibble(myGSEA.res2@result)

write.csv(myGSEA.df2, "LimmabyFC_mGO_mSigDB_GSEA.csv")

gseaplot2(myGSEA.res,
          geneSetID = c("GOBP_CELL_JUNCTION_ASSEMBLY", "GOBP_CELL_JUNCTION_ORGANIZATION", "GOBP_REGULATION_OF_CELL_JUNCTION_ASSEMBLY", "GOCC_ANCHORING_JUNCTION"), #can choose multiple signatures to overlay in this plot
          pvalue_table = FALSE, #can set this to FALSE for a cleaner plot
          #title = myGSEA.res$Description[14]
) #can also turn off this title

A#Atry GSEA with up and down
mydata.up.all.ord <- mydata_UP[order(mydata_UP$logFC, decreasing = T),]
mydata.down.all.ord <- mydata_DOWN[order(mydata_DOWN$logFC, decreasing = T),]
mydata.up.fc1.5.ord <- mydata_UP_FC1.5[order(mydata_UP_FC1.5$logFC, decreasing = T),]
mydata.down.fc1.5.ord <- mydata_DOWN_FC1.5[order(mydata_DOWN_FC1.5$logFC, decreasing = T),]


mydata.up.all.ord.tibble <- as_tibble(mydata.up.all.ord, rownames = "GeneID")
mydata.df.sub3 <- dplyr::select(mydata.up.all.ord.tibble, GeneID, logFC)
mydata.gsea3 <- mydata.df.sub3$logFC
names(mydata.gsea3) <- as.character(mydata.df.sub3$GeneID)


set.seed(125) #set a random seed so that we can reproducible ordering for our GSEA results below
myGSEA.res3 <- GSEA(mydata.gsea3, TERM2GENE=C2CP, verbose=FALSE) #could replace C2CP with hs_gsea_c2 object you retrieved from msigdb above
beep(sound = 1)

myGSEA.res.df3 <- as_tibble(myGSEA.res3@result)
rm(myGSEA.res.df3)
rm(myGSEA.res3)
rm(mydata.gsea3)
rm(mydata.df.sub3)
rm(mydata.up.all.ord.tibble)


myGSEA.df <- myGSEA.df %>%
  mutate(phenotype = case_when(
    NES > 0 ~ "Treated",
    NES < 0 ~ "Control"))

GOBP_VASCULATURE_DEVELOPMENT
GOBP_CORONARY_VASCULATURE_MORPHOGENESIS
GOBP_VASCULAR_ENDOTHELIAL_CELL_PROLIFERATION
GOBP_CELLULAR_RESPONSE_TO_VASCULAR_ENDOTHELIAL_GROWTH_FACTOR_STIMULUS
GOBP_VASOCONSTRICTION
GOBP_VASCULAR_ENDOTHELIAL_GROWTH_FACTOR_SIGNALING_PATHWAY
GOBP_VASCULATURE_DEVELOPMENT
GOBP_VASCULAR_ENDOTHELIAL_GROWTH_FACTOR_RECEPTOR_SIGNALING_PATHWAY
GOBP_REGULATION_OF_CELL_JUNCTION_ASSEMBLY
GOBP_CELL_JUNCTION_ASSEMBLY
GOBP_CELL_JUNCTION_ORGANIZATION
GOCC_ANCHORING_JUNCTION
GOBP_MESENCHYMAL_CELL_DIFFERENTIATION
GOBP_MESENCHYME_DEVELOPMENT
GOMF_EXTRACELLULAR_MATRIX_STRUCTURAL_CONSTITUENT
GOCC_COLLAGEN_CONTAINING_EXTRACELLULAR_MATRIX
GOMF_EXTRACELLULAR_MATRIX_STRUCTURAL_CONSTITUENT_CONFERRING_TENSILE_STRENGTH
GOBP_EXTRACELLULAR_MATRIX_ASSEMBLY
GOBP_REGULATION_OF_EXTRACELLULAR_MATRIX_ORGANIZATION
GOBP_COLLAGEN_METABOLIC_PROCESS
GOCC_COLLAGEN_TRIMER
GOMF_COLLAGEN_BINDING
GOBP_COLLAGEN_FIBRIL_ORGANIZATION
GOBP_COLLAGEN_CATABOLIC_PROCESS
GOBP_REGULATION_OF_COLLAGEN_METABOLIC_PROCESS
GOBP_COLLAGEN_BIOSYNTHETIC_PROCESS

rownames(myGSEA.df)
myGSEA.df[1,]

ggplot(myGSEA.df[1:25,], aes(x=phenotype, y=ID)) +
  geom_point(aes(size=setSize, color = NES, alpha=-log10(p.adjust))) +
  scale_color_gradient(low="blue", high="red") +
  theme_bw()
