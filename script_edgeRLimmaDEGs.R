# Load packages -----
library(tidyverse) # you know it well by now!
library(limma) # venerable package for differential gene expression using linear modeling
library(edgeR)
library(gt)
library(DT)
library(plotly)
library(ggrepel)

# Set up your design matrix ----
group <- factor(targets$Condition)
design <- model.matrix(~0 + group)
colnames(design) <- levels(group)

# NOTE: if you need a paired analysis (a.k.a.'blocking' design) or have a batch effect, the following design is useful
# design <- model.matrix(~block + treatment)
# this is just an example. 'block' and 'treatment' would need to be objects in your environment

# Model mean-variance trend and fit linear model to data ----
# Use VOOM function from Limma package to model the mean-variance relationship
v.DList.filtered.norm_0 <- voom(DList.filtered.norm_0, design, plot = TRUE)
v.DList.filtered.norm_1 <- voom(DList.filtered.norm_1, design, plot = TRUE)
v.DList.filtered.norm_10 <- voom(DList.filtered.norm_10, design, plot = TRUE)
view(v.DList.filtered.norm_0)
view(DList.filtered.norm_0)
DList.filtered.norm_0.df <- v.DList.filtered.norm_0
colnames(DList.filtered.norm_0.df) <- c("geneID", sampleLabels)
write.csv(v.DList.filtered.norm_0, "DList.filtered.norm_0.df.csv")

# fit a linear model to your data
fit_0 <- lmFit(v.DList.filtered.norm_0, design)
fit_1 <- lmFit(v.DList.filtered.norm_1, design)
fit_10 <- lmFit(v.DList.filtered.norm_10, design)
view(fit_0)
contrast.matrix <- makeContrasts(Treated-Control,
                                 levels=design)

fits_0 <- contrasts.fit(fit_0, contrast.matrix)
fits_1 <- contrasts.fit(fit_1, contrast.matrix)
fits_10 <- contrasts.fit(fit_10, contrast.matrix)
view(fits_0)
ebFit_0 <- eBayes(fits_0)
ebFit_1 <- eBayes(fits_1)
ebFit_10 <- eBayes(fits_10)
view(ebFit_0)
myTopHits_0 <- topTable(ebFit_0, adjust ="BH", coef=1, number=40000, sort.by="logFC")
myTopHits_0.tibble <- as_tibble(myTopHits_0, rownames = "Gene_ID")
linfit.tibble_0 <- as_tibble(v.DList.filtered.norm_0$E, rownames = "Gene_ID")

?eBayes

myTopHits_1 <- topTable(ebFit_1, adjust ="BH", coef=1, number=40000, sort.by="logFC")
myTopHits_1.tibble <- as_tibble(myTopHits_1, rownames = "Gene_ID")
linfit.tibble_1 <- as_tibble(v.DList.filtered.norm_1$E, rownames = "Gene_ID")

myTopHits_10 <- topTable(ebFit_10, adjust ="BH", coef=1, number=40000, sort.by="logFC")
myTopHits_10.tibble <- as_tibble(myTopHits_10, rownames = "Gene_ID")
linfit.tibble_10 <- as_tibble(v.DList.filtered.norm_10$E, rownames = "Gene_ID")

myTopHits.df <- myTopHits_0 %>%
  as_tibble(rownames = "geneID")

view(v.Dlist.df)
view(myTopHits.df)
view(myTopHits_0)
write.csv(myTopHits.df, "myTopHits_ER_Limma.csv")
write.csv(v.DList.filtered.norm_0 , "v.DList.filtered.norm_0.csv")

myTopHits.df2 <- myTopHits.df
myTopHits.df2$diffexp <- "NO"
myTopHits.df2$diffexp[myTopHits.df2$adj.P.Val<0.05 & myTopHits.df2$logFC>log2(1.5)] <- "UP"
myTopHits.df2$diffexp[myTopHits.df2$adj.P.Val<0.05 & myTopHits.df2$logFC < -log2(1.5)] <- "DOWN"
myTopHits.df2$diffexp[myTopHits.df2$adj.P.Val<0.05 & myTopHits.df2$logFC>0 & myTopHits.df2$logFC<=log2(1.5)] <- "UPcut"
myTopHits.df2$diffexp[myTopHits.df2$adj.P.Val<0.05 & myTopHits.df2$logFC<0 & myTopHits.df2$logFC >= -log2(1.5)] <- "DOWNCut"
top20degs <- head(myTopHits.df2[order(myTopHits.df2$adj.P.Val), 'geneID'], 20)
GenesOfChoice <- data.frame(geneID = c("Tjp1", 
                                       "Cldn5", 
                                       "Cldn25", 
                                       "Pecam1", 
                                       "Hmox1",
                                       "Ahr",
                                       "Cyp1b1",
                                       "Foxf2",
                                       "Tal1",
                                       "Col1a1",
                                       "Col1a2",
                                       "Col3a1",
                                       "Col5a1",
                                       "Col6a1",
                                       "Col6a2",
                                       "Col6a3",
                                       "Col6a6",
                                       "Axin2",
                                       "Tcf7",
                                       "Angptl7",
                                       "Flrt1",
                                       "Dock8",
                                       "Plekha4",
                                       "Gm14567",
                                       "Cyp2c55",
                                       "Sprr2a2",
                                       "Fxyd3",
                                       "Sez6l"
                                       ))
colnames(GenesOfChoice) <- "geneID"
gList <- rbind(top20degs,GenesOfChoice)
?ggplot
which(myTopHits.df2$geneID,)
myTopHits.df2$dlabel <- ifelse(myTopHits.df2$geneID %in% GenesOfChoice$geneID, myTopHits.df2$geneID, NA)
dim(gList)
view(myTopHits.df2)
view(top20degs)
view(GenesOfChoice)
view(gList)

vplot <- ggplot(myTopHits.df2) +
  aes(y=-log10(adj.P.Val), x=logFC, text = paste("Symbol:", geneID), col = diffexp, label = dlabel) +
  geom_point(size=2) +
  geom_hline(yintercept = -log10(0.05), linetype="longdash", colour="grey", linewidth=1) +
  geom_vline(xintercept = log2(1.5), linetype="longdash", colour="#BE684D", linewidth=1) +
  geom_vline(xintercept = -log2(1.5), linetype="longdash", colour="#BE684D", linewidth=1) +
  annotate("rect", xmin = log2(1.5), xmax = 12, ymin = -log10(0.05), ymax = 12, alpha=.2, fill="#BE684D") +
  annotate("rect", xmin = -log2(1.5), xmax = -12, ymin = -log10(0.05), ymax = 12, alpha=.2, fill="#2C467A") +
  labs(title="Volcano plot",
       subtitle = "UAS treatment",
       caption=paste0("q<0.05 and -1.5<FC<1.5 - EdgeR Limma filter (cpm>0)>=3")) +
  #theme_bw()
  scale_color_manual(values = c("blue","cyan","black", "red", "darkgoldenrod1"))
  geom_text_repel(max.overlaps = inf)
  
  
#2C467A

vplot <- ggplot(myTopHits.df2) +
    aes(y=-log10(adj.P.Val), x=logFC, col = diffexp, label = dlabel) +
    geom_point(size=2) +
    geom_hline(yintercept = -log10(0.05), linetype="longdash", colour="grey", linewidth=1) +
    geom_vline(xintercept = log2(1.5), linetype="longdash", colour="#BE684D", linewidth=1) +
    geom_vline(xintercept = -log2(1.5), linetype="longdash", colour="#BE684D", linewidth=1) +
    annotate("rect", xmin = log2(1.5), xmax = 12, ymin = -log10(0.05), ymax = 12, alpha=.2, fill="#BE684D") +
    annotate("rect", xmin = -log2(1.5), xmax = -12, ymin = -log10(0.05), ymax = 12, alpha=.2, fill="#2C467A") +
    #geom_text(color="black")+
    #labs(title="Volcano plot",
         #subtitle = "UAS treatment",
         #caption=paste0("q<0.05 and -1.5<FC<1.5 - EdgeR Limma filter (cpm>0)>=3")) +
    #theme_bw()
    scale_color_manual(values = c("royalblue2","lightskyblue2","grey", "indianred2", "darkgoldenrod1"))+
    geom_text_repel(max.overlaps = Inf, color = "black")


ggplot(myTopHits.df2) +
  aes(y=-log10(adj.P.Val), x=logFC, col = diffexp) +
  geom_point(size=2) +
  geom_hline(yintercept = -log10(0.05), linetype="longdash", colour="grey", linewidth=1) +
  geom_vline(xintercept = log2(1.5), linetype="longdash", colour="#BE684D", linewidth=1) +
  geom_vline(xintercept = -log2(1.5), linetype="longdash", colour="#BE684D", linewidth=1) +
  annotate("rect", xmin = log2(1.5), xmax = 12, ymin = -log10(0.05), ymax = 12, alpha=.2, fill="#BE684D") +
  annotate("rect", xmin = -log2(1.5), xmax = -12, ymin = -log10(0.05), ymax = 12, alpha=.2, fill="#2C467A") +
  #geom_text(color="black")+
  #labs(title="Volcano plot",
  #subtitle = "UAS treatment",
  #caption=paste0("q<0.05 and -1.5<FC<1.5 - EdgeR Limma filter (cpm>0)>=3")) +
  #theme_bw()
  scale_color_manual(values = c("royalblue2","lightskyblue2","grey", "indianred2", "darkgoldenrod1"))+
  #geom_text_repel(max.overlaps = Inf, color = "black")