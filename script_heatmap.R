library(tidyverse)
library(limma) #we only use limma in this script for the 'avearrays' function
library(RColorBrewer) #need colors to make heatmaps
library(gplots) #the heatmap2 function in this package is a primary tool for making heatmaps
library(gameofthrones) #because...why not.  Install using 'devtools::install_github("aljrico/gameofthrones")'
library(d3heatmap) # install from github with: devtools::install_github("talgalili/d3heatmap")
library(heatmaply)
#library(tblhelpr)
#library(magrittr)

myheatcolors1 <- redblue(75) #this is from the 'colorpanel' function in gplots (same package that heatmap.2 comes from)


results <- decideTests(ebFit_0, method="global", adjust.method="BH", p.value=0.05, lfc=log2(1.5))
results.all.sig <- decideTests(ebFit_0, method="global", adjust.method="BH", p.value=0.05)
colnames(v.DList.filtered.norm_0$E) <- sampleLabels
diffGenes <- v.DList.filtered.norm_0$E[results[,1] !=0,]
dim(diffGenes)
view(diffGenes)
view(ebFit_0)
view(v.DList.filtered.norm_0)
dim(x)
view(results[.1] == 1)
upgenesFC1.5 <- v.DList.filtered.norm_0$E[results[,1] == 1,]
dim(upgenesFC1.5)
downgenesFC1.5 <- v.DList.filtered.norm_0$E[results[,1] == -1,]
dim(downgenesFC1.5)
nsFCtgenes <- v.DList.filtered.norm_0$E[results[,1] == 0,]
dim(nsFCtgenes)
All.UP <- v.DList.filtered.norm_0$E[results.all.sig[,1] == 1,] 
All.DOWN <- v.DList.filtered.norm_0$E[results.all.sig[,1] == -1,] 
All.nonsig <- v.DList.filtered.norm_0$E[results.all.sig[,1] == 0,] 
dim(All.UP)
dim(All.DOWN)
dim(All.nonsig)


clustRows <- hclust(as.dist(1-cor(t(diffGenes), method="pearson")), method="complete") 
clustColumns <- hclust(as.dist(1-cor(diffGenes, method="spearman")), method="complete") #cluster columns by spearman correlation
module.assign <- cutree(clustRows, k=2)

module.color <- rainbow(length(unique(module.assign)), start=0.1, end=0.9) 
module.color <- module.color[as.vector(module.assign)] 


?ggplot

heatmap.2(diffGenes, 
          Rowv=as.dendrogram(clustRows), 
          Colv=as.dendrogram(clustColumns),
          RowSideColors=module.color,
          col=rev(myheatcolors1), scale='row',
          density.info="none", trace="none",  
          cexRow=0.05, cexCol=1, margins=c(8,20) )
          
#diffGenes <- v.DList.filtered.norm_0$E
#diffGenes <- data.matrix(diffGenes)
#rownames(diffGenes)
#diffGenes[,]
 
d3heatmap(diffGenes,
          colors = myheatcolors1,
          Rowv=as.dendrogram(clustRows),
          row_side_colors = module.color,
          scale='row')

names(module.color) <- names(module.assign) 

module.assign.df <- as_tibble(as.list(module.assign))
module.assign.pivot <- pivot_longer(module.assign.df, # dataframe to be pivoted
                                    cols = everything(), # column names to be stored as a SINGLE variable
                                    names_to = "geneID", # name of that new variable (column)
                                    values_to = "module") # name of new variable (column) storing all the values (data)

view(module.assign)
