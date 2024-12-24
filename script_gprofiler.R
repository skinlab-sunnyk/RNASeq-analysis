library(tidyverse)
library(limma)
library(gplots) #for heatmaps
library(DT) #interactive and searchable tables of our GSEA results
library(GSEABase) #functions and methods for Gene Set Enrichment Analysis
library(Biobase) #base functions for bioconductor; required by GSEABase
library(GSVA) #Gene Set Variation Analysis, a non-parametric and unsupervised method for estimating variation of gene set enrichment across samples.
library(gprofiler2) #tools for accessing the GO enrichment results using g:Profiler web resources
library(clusterProfiler) # provides a suite of tools for functional enrichment analysis
library(msigdbr) # access to msigdb collections directly within R
library(enrichplot) # great for making the standard GSEA enrichment plots
library(beepr)
# Carry out GO enrichment using gProfiler2 ----
# use topTable result to pick the top genes for carrying out a Gene Ontology (GO) enrichment analysis
#myTopHits_0 <- topTable(ebFit, adjust ="BH", coef=1, number=50, sort.by="logFC")
# use the 'gost' function from the gprofiler2 package to run GO enrichment analysis

mydata.up.all.ord <- mydata_UP[order(mydata_UP$logFC, decreasing = T),]
mydata.down.all.ord <- mydata_DOWN[order(mydata_DOWN$logFC, decreasing = T),]
mydata.up.fc1.5.ord <- mydata_UP_FC1.5[order(mydata_UP_FC1.5$logFC, decreasing = T),]
mydata.down.fc1.5.ord <- mydata_DOWN_FC1.5[order(mydata_DOWN_FC1.5$logFC, decreasing = T),]

dim(mydata.down.fc1.5.ord)
rownames(mydata.down.all.ord)

gost.res.up.all <- gost(rownames(mydata.up.all.ord), organism = "mmusculus", correction_method = "fdr")
beep(sound =1)
?gost()

gost.res.down.all <- gost(rownames(mydata.down.all.ord), organism = "mmusculus", correction_method = "fdr")
beep(sound =1)

gost.res.up.fc1.5 <- gost(rownames(mydata.up.fc1.5.ord), organism = "mmusculus", correction_method = "fdr")
beep(sound =1)

gost.res.down.fc1.5 <- gost(rownames(mydata.down.fc1.5.ord), organism = "mmusculus", correction_method = "fdr")
beep(sound =1)

dim(mydata.down.fc1.5.ord)
Gprofiler.res.down.1.5 <- gost.res.down.fc1.5$result

# produce an interactive manhattan plot of enriched GO terms
gostplot(gost.res, interactive = TRUE, capped = TRUE)
#set interactive=FALSE to get plot for publications
mygostplot <- gostplot(gost.res, interactive = FALSE, capped = TRUE)
# produce a publication quality static manhattan plot with specific GO terms highlighted
publish_gostplot(
  mygostplot, #your static gostplot from above
  highlight_terms = c("REAC:R-HSA-9662851", "REAC:R-HSA-9824443", "REAC:R-HSA-9664407", "REAC:R-HSA-9658195", "REAC:R-HSA-9664417"),
  filename = NULL,
  width = NA,
  height = NA)
