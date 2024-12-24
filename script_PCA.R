# Load packages ------
library(tidyverse) # you're familiar with this fromt the past two lectures
library(DT) # for making interactive tables
library(plotly) # for making interactive plots
library(gt) # A layered 'grammar of tables' - think ggplot, but for tables
# Identify variables of interest in study design file ----
targets
group <- targets$Condition
group <- factor(group)

log2.cpm.filtered.norm.df_0
log2.cpm.filtered.norm.df_1
log2.cpm.filtered.norm.df_10

# Hierarchical clustering ---------------
distance_0 <- dist(t(log2.cpm.filtered.norm_0), method = "maximum") #other distance methods are "euclidean", maximum", "manhattan", "canberra", "binary" or "minkowski"
distance_1 <- dist(t(log2.cpm.filtered.norm_1), method = "maximum") #other distance methods are "euclidean", maximum", "manhattan", "canberra", "binary" or "minkowski"
distance_10 <- dist(t(log2.cpm.filtered.norm_10), method = "maximum") #other distance methods are "euclidean", maximum", "manhattan", "canberra", "binary" or "minkowski"

clusters_0 <- hclust(distance_0, method = "average") #other agglomeration methods are "ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median", or "centroid"
clusters_1 <- hclust(distance_1, method = "average") #other agglomeration methods are "ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median", or "centroid"
clusters_10 <- hclust(distance_10, method = "average") #other agglomeration methods are "ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median", or "centroid"

plot(clusters_0, labels=sampleLabels)
plot(clusters_1, labels=sampleLabels)
plot(clusters_10, labels=sampleLabels)

# Principal component analysis (PCA) -------------
pca.res_0 <- prcomp(t(log2.cpm.filtered.norm_0), scale.=F, retx=T)
pca.res_1 <- prcomp(t(log2.cpm.filtered.norm_1), scale.=F, retx=T)
pca.res_10 <- prcomp(t(log2.cpm.filtered.norm_10), scale.=F, retx=T)

#look at the PCA result (pca.res) that you just created
ls(pca.res_0)
ls(pca.res_1)
ls(pca.res_10)

summary(pca.res_0) # Prints variance summary for all principal components.
summary(pca.res_1) # Prints variance summary for all principal components.
summary(pca.res_10) # Prints variance summary for all principal components.

#pca.res_0$rotation #$rotation shows you how much each gene influenced each PC (called 'scores')
#pca.res$x # 'x' shows you how much each sample influenced each PC (called 'loadings')
#note that these have a magnitude and a direction (this is the basis for making a PCA plot)
screeplot(pca.res_0) # A screeplot is a standard way to view eigenvalues for each PCA
screeplot(pca.res_1) # A screeplot is a standard way to view eigenvalues for each PCA
screeplot(pca.res_10) # A screeplot is a standard way to view eigenvalues for each PCA


pc.var_0<-pca.res_0$sdev^2 # sdev^2 captures these eigenvalues from the PCA result
pc.per_0<-round(pc.var_0/sum(pc.var_0)*100, 1) # we can then use these eigenvalues to calculate the percentage variance explained by each PC
pc.per_0

pc.var_1<-pca.res_1$sdev^2 # sdev^2 captures these eigenvalues from the PCA result
pc.per_1<-round(pc.var_1/sum(pc.var_1)*100, 1) # we can then use these eigenvalues to calculate the percentage variance explained by each PC
pc.per_1

pc.var_10<-pca.res_10$sdev^2 # sdev^2 captures these eigenvalues from the PCA result
pc.per_10<-round(pc.var_10/sum(pc.var_10)*100, 1) # we can then use these eigenvalues to calculate the percentage variance explained by each PC
pc.per_10


# Visualize your PCA result ------------------
#lets first plot any two PCs against each other
#We know how much each sample contributes to each PC (loadings), so let's plot
pca.res.df_0 <- as_tibble(pca.res_0$x)
pca.res.df_1 <- as_tibble(pca.res_1$x)
pca.res.df_10 <- as_tibble(pca.res_10$x)


ggplot(pca.res.df_0) +
  aes(x=PC1, y=PC2, label=sampleLabels, color=group) +
  geom_point(size=4) +
  #geom_label() +
  #stat_ellipse() +
  xlab(paste0("PC1 (",pc.per_0[1],"%",")")) + 
  ylab(paste0("PC2 (",pc.per_0[2],"%",")")) +
  labs(title="PCA plot",
       caption=paste0("Based on filter (cpm>0)>=3 ", Sys.time())) +
  coord_fixed() +
  theme_bw()

ggplot(pca.res.df_1) +
  aes(x=PC1, y=PC2, label=sampleLabels, color=group) +
  geom_point(size=4) +
  #geom_label() +
  #stat_ellipse() +
  xlab(paste0("PC1 (",pc.per_1[1],"%",")")) + 
  ylab(paste0("PC2 (",pc.per_1[2],"%",")")) +
  labs(title="PCA plot",
       caption=paste0("Based on filter (cpm>1)>=3 ", Sys.time())) +
  coord_fixed() +
  theme_bw()

ggplot(pca.res.df_10) +
  aes(x=PC1, y=PC2, label=sampleLabels, color=group) +
  geom_point(size=4) +
  #geom_label() +
  #stat_ellipse() +
  xlab(paste0("PC1 (",pc.per_10[1],"%",")")) + 
  ylab(paste0("PC2 (",pc.per_10[2],"%",")")) +
  labs(title="PCA plot",
       caption=paste0("Based on filter (cpm>10)>=3 ", Sys.time())) +
  coord_fixed() +
  theme_bw()

# Create a PCA 'small multiples' chart for cpm>0----

pca.res.df_0 <- pca.res_0$x[,1:6] %>% # note that this is the first time you've seen the 'pipe' operator from the magrittr package
  as_tibble() %>%
  add_column(sample = sampleLabels,
             group = group)

pca.pivot_0 <- pivot_longer(pca.res.df_0, # dataframe to be pivoted
                          cols = PC1:PC6, # column names to be stored as a SINGLE variable
                          names_to = "PC", # name of that new variable (column)
                          values_to = "loadings") # name of new variable (column) storing all the values (data)

ggplot(pca.pivot_0) +
  aes(x=sample, y=loadings, fill=group) + # you could iteratively 'paint' different covariates onto this plot using the 'fill' aes
  geom_bar(stat="identity") +
  facet_wrap(~PC) +
  labs(title="PCA 'small multiples' plot",
       caption=paste0("based on (cpm>0)>=3 ", Sys.time())) +
  theme_bw() +
  coord_flip()

# Create a PCA 'small multiples' chart for cpm>1----

pca.res.df_1 <- pca.res_1$x[,1:6] %>% # note that this is the first time you've seen the 'pipe' operator from the magrittr package
  as_tibble() %>%
  add_column(sample = sampleLabels,
             group = group)

pca.pivot_1 <- pivot_longer(pca.res.df_1, # dataframe to be pivoted
                            cols = PC1:PC6, # column names to be stored as a SINGLE variable
                            names_to = "PC", # name of that new variable (column)
                            values_to = "loadings") # name of new variable (column) storing all the values (data)

ggplot(pca.pivot_1) +
  aes(x=sample, y=loadings, fill=group) + # you could iteratively 'paint' different covariates onto this plot using the 'fill' aes
  geom_bar(stat="identity") +
  facet_wrap(~PC) +
  labs(title="PCA 'small multiples' plot",
       caption=paste0("based on (cpm>1)>=3 ", Sys.time())) +
  theme_bw() +
  coord_flip()

# Create a PCA 'small multiples' chart for cpm>1----

pca.res.df_10 <- pca.res_10$x[,1:6] %>% # note that this is the first time you've seen the 'pipe' operator from the magrittr package
  as_tibble() %>%
  add_column(sample = sampleLabels,
             group = group)

pca.pivot_10 <- pivot_longer(pca.res.df_10, # dataframe to be pivoted
                            cols = PC1:PC6, # column names to be stored as a SINGLE variable
                            names_to = "PC", # name of that new variable (column)
                            values_to = "loadings") # name of new variable (column) storing all the values (data)

ggplot(pca.pivot_10) +
  aes(x=sample, y=loadings, fill=group) + # you could iteratively 'paint' different covariates onto this plot using the 'fill' aes
  geom_bar(stat="identity") +
  facet_wrap(~PC) +
  labs(title="PCA 'small multiples' plot",
       caption=paste0("based on (cpm>10)>=3 ", Sys.time())) +
  theme_bw() +
  coord_flip()

beep(sound = 1)

