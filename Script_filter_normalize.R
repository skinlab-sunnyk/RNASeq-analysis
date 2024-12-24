#Introduction----
#script 2 
#Filter and normalize the data
#Visulize the effect of filtering and normalization

# Load packages -----
library(tidyverse) # already know about this from Step 1 script
library(edgeR) # well known package for differential expression analysis, but we only use for the DGEList object and for normalization methods
library(matrixStats) # let's us easily calculate stats on rows or columns of a data matrix
library(cowplot) # allows you to combine multiple plots in one figure
library(ggplot2)

#Examine data----
TPM <- Txi_gene$abundance
Counts <- Txi_gene$counts
colSums(TPM)
colSums(Counts)

# capture sample labels from the study design file that you worked with and saved as 'targets' in step 1
targets
sampleLabels <- targets$Sample_type

TPM.stats <- transform(TPM, 
                         SD=rowSds(TPM), 
                         AVG=rowMeans(TPM),
                         MED=rowMedians(TPM))

Counts.stats <- transform(Counts, 
                         SD=rowSds(Counts), 
                         AVG=rowMeans(Counts),
                         MED=rowMedians(Counts))



head(TPM.stats)
head(Counts.stats)

Counts.stats.tibble <- as_tibble(Counts.stats, rownames = "GeneID")
view(Counts.stats.tibble)
write.csv(Counts.stats.tibble, "counts.stats.tibble.")

TPM.stats.tibble <- as_tibble(TPM.stats, rownames = "GeneID")
view(TPM.stats.tibble)
write.csv(TPM.stats.tibble, "TPM.stats.tibble.")
#visualize----
p2.raw.tpm <- ggplot(TPM.stats) +
  aes(x=SD, y=MED) +
  geom_point(shape = 16, size=2) +
  geom_smooth(method=lm) +
  labs(y="Median", x = "SD",
       title="Transcripts per million (TPM; Abundance)",
       subtitle="Raw abundance:unfiltered, non-normalized data",
       caption = paste0("Binita's data ", Sys.time())) +
  theme_bw()

p2.raw.counts <- ggplot(Counts.stats) +
  aes(x=SD, y=MED) +
  geom_point(shape = 16, size=2) +
  geom_smooth(method=lm) +
  labs(y="Median", x = "SD",
       title="Raw counts",
       subtitle="unfiltered, non-normalized data",
       caption = paste0("Binita's data ", Sys.time())) +
  theme_bw()  

#Use edgeR to create DGElist

DList <- DGEList(Counts)
save(DList, file = "DGEList")
view(DList)
cpm <- cpm(DList) 
colSums(cpm)

log2.cpm <- cpm(DList, log=TRUE)
colSums(log2.cpm)
?cpm()

log2.cpm.df <- as_tibble(log2.cpm, rownames = "geneID")
cpm.df <- as_tibble(cpm, rownames = "geneID")
view(cpm.df)

colnames(log2.cpm.df) <- c("geneID", sampleLabels)
colnames(cpm.df) <- c("geneID", sampleLabels)

log2.cpm.df.pivot <- pivot_longer(log2.cpm.df, # dataframe to be pivoted
                                  cols = C1:T3, # column names to be stored as a SINGLE variable
                                  names_to = "samples", # name of that new variable (column)
                                  values_to = "expression") # name of new variable (column) storing all the values (data)

cpm.df.pivot <- pivot_longer(cpm.df, # dataframe to be pivoted
                             cols = C1:T3, # column names to be stored as a SINGLE variable
                             names_to = "samples", # name of that new variable (column)
                             values_to = "expression") # name of new variable (column) storing all the values (data)

DList$counts

p.raw.log2cpm <- ggplot(log2.cpm.df.pivot) +
  aes(x=samples, y=expression, fill=samples) +
  geom_violin(trim = FALSE, show.legend = FALSE) +
  stat_summary(fun = "median", 
               geom = "point", 
               shape = 95, 
               size = 10, 
               color = "black", 
               show.legend = FALSE) +
  labs(y="log2 expression", x = "sample",
       title="Log2 Counts per Million (CPM)",
       subtitle="unfiltered, non-normalized",
       caption=paste0("produced on ", Sys.time())) +
  theme_bw()

p.raw.cpm <- ggplot(cpm.df.pivot) +
  aes(x=samples, y=expression, fill=samples) +
  geom_violin(trim = FALSE, show.legend = FALSE) +
  stat_summary(fun = "median", 
               geom = "point", 
               shape = 95, 
               size = 10, 
               color = "black", 
               show.legend = FALSE) +
  labs(y="expression", x = "sample",
       title="Raw Counts per Million (CPM)",
       subtitle="unfiltered, non-normalized",
       caption=paste0("produced on ", Sys.time())) +
  theme_bw()

#filter the data----
table(rowSums(DList$counts==0)==6)
keepers1 <- rowSums(cpm>1)>=3
summary(keepers1)
keepers0 <- rowSums(cpm>0)>=3
summary(keepers0)
keepers10 <- rowSums(cpm>10)>=3
summary(keepers0)


#keepers0_2 <- rowSums(cpm>0)>2
#summary(keepers0_2)
#keepers0_4 <- rowSums(cpm>0)>4
#summary(keepers0_2)

DList.filtered_1 <- DList[keepers1,]
DList.filtered_0 <- DList[keepers0,]
DList.filtered_10 <- DList[keepers10,]

#DList.filtered_0_2 <- DList[keepers0_2,]
#DList.filtered_0_4 <- DList[keepers0_4,]

dim(DList.filtered_1)
dim(DList.filtered_0)
dim(DList.filtered_10)
#dim(DList.filtered_0_2)
#dim(DList.filtered_0_4)

log2.cpm.filtered_1 <- cpm(DList.filtered_1, log=TRUE)
log2.cpm.filtered_0 <- cpm(DList.filtered_0, log=TRUE)
log2.cpm.filtered_10 <- cpm(DList.filtered_10, log=TRUE)
view(DList.filtered_0)

#visulize filtered data----

log2.cpm.filtered.df_1 <- as_tibble(log2.cpm.filtered_1, rownames = "geneID")
log2.cpm.filtered.df_0 <- as_tibble(log2.cpm.filtered_0, rownames = "geneID")
log2.cpm.filtered.df_10 <- as_tibble(log2.cpm.filtered_10, rownames = "geneID")

view(log2.cpm.filtered.df_0)

colnames(log2.cpm.filtered.df_1) <- c("geneID", sampleLabels)
colnames(log2.cpm.filtered.df_0) <- c("geneID", sampleLabels)
colnames(log2.cpm.filtered.df_10) <- c("geneID", sampleLabels)


log2.cpm.filtered.df.pivot_1 <- pivot_longer(log2.cpm.filtered.df_1, # dataframe to be pivoted
                                           cols = C1:T3, # column names to be stored as a SINGLE variable
                                           names_to = "samples", # name of that new variable (column)
                                           values_to = "expression") # name of new variable (column) storing all the values (data)


log2.cpm.filtered.df.pivot_0 <- pivot_longer(log2.cpm.filtered.df_0, # dataframe to be pivoted
                                             cols = C1:T3, # column names to be stored as a SINGLE variable
                                             names_to = "samples", # name of that new variable (column)
                                             values_to = "expression") # name of new variable (column) storing all the values (data)

log2.cpm.filtered.df.pivot_10 <- pivot_longer(log2.cpm.filtered.df_10, # dataframe to be pivoted
                                             cols = C1:T3, # column names to be stored as a SINGLE variable
                                             names_to = "samples", # name of that new variable (column)
                                             values_to = "expression") # name of new variable (column) storing all the values (data)


p.log2cpmfilted_1 <- ggplot(log2.cpm.filtered.df.pivot_1) +
  aes(x=samples, y=expression, fill=samples) +
  geom_violin(trim = FALSE, show.legend = FALSE) +
  stat_summary(fun = "median", 
               geom = "point", 
               shape = 95, 
               size = 10, 
               color = "black", 
               show.legend = FALSE) +
  labs(y="log2 expression", x = "sample",
       title="Log2 Counts per Million (CPM)",
       subtitle="filtered for cpm>1, non-normalized",
       caption=paste0("produced on ", Sys.time())) +
  theme_bw()

p.log2cpmfilted_0 <- ggplot(log2.cpm.filtered.df.pivot_0) +
  aes(x=samples, y=expression, fill=samples) +
  geom_violin(trim = FALSE, show.legend = FALSE) +
  stat_summary(fun = "median", 
               geom = "point", 
               shape = 95, 
               size = 10, 
               color = "black", 
               show.legend = FALSE) +
  labs(y="log2 expression", x = "sample",
       title="Log2 Counts per Million (CPM)",
       subtitle="filtered for cpm>0, non-normalized",
       caption=paste0("produced on ", Sys.time())) +
  theme_bw()

p.log2cpmfilted_10 <- ggplot(log2.cpm.filtered.df.pivot_10) +
  aes(x=samples, y=expression, fill=samples) +
  geom_violin(trim = FALSE, show.legend = FALSE) +
  stat_summary(fun = "median", 
               geom = "point", 
               shape = 95, 
               size = 10, 
               color = "black", 
               show.legend = FALSE) +
  labs(y="log2 expression", x = "sample",
       title="Log2 Counts per Million (CPM)",
       subtitle="filtered for cpm>10, non-normalized",
       caption=paste0("produced on ", Sys.time())) +
  theme_bw()



# Normalize your data ----
DList.filtered.norm_1 <- calcNormFactors(DList.filtered_1, method = "TMM")
DList.filtered.norm_0 <- calcNormFactors(DList.filtered_0, method = "TMM")
DList.filtered.norm_10 <- calcNormFactors(DList.filtered_10, method = "TMM")

?calcNormFactors()
DList.filtered.norm_0.df <- as_tibble(DList.filtered.norm_0, rownames = "geneID")

log2.cpm.filtered.norm_1 <- cpm(DList.filtered.norm_1, log=TRUE)
log2.cpm.filtered.norm_0 <- cpm(DList.filtered.norm_0, log=TRUE)
log2.cpm.filtered.norm_10 <- cpm(DList.filtered.norm_10, log=TRUE)

#Visulaize the filtered and normalized data----

log2.cpm.filtered.norm.df_1 <- as_tibble(log2.cpm.filtered.norm_1, rownames = "geneID")
log2.cpm.filtered.norm.df_0 <- as_tibble(log2.cpm.filtered.norm_0, rownames = "geneID")
log2.cpm.filtered.norm.df_10 <- as_tibble(log2.cpm.filtered.norm_10, rownames = "geneID")


colnames(log2.cpm.filtered.norm.df_1) <- c("geneID", sampleLabels)
colnames(log2.cpm.filtered.norm.df_0) <- c("geneID", sampleLabels)
colnames(log2.cpm.filtered.norm.df_10) <- c("geneID", sampleLabels)

view(log2.cpm.filtered.norm.df_0)
log2.cpm.filtered.norm.df.pivot_1 <- pivot_longer(log2.cpm.filtered.norm.df_1, # dataframe to be pivoted
                                                cols = C1:T3, # column names to be stored as a SINGLE variable
                                                names_to = "samples", # name of that new variable (column)
                                                values_to = "expression") # name of new variable (column) storing all the values (data)

log2.cpm.filtered.norm.df.pivot_0 <- pivot_longer(log2.cpm.filtered.norm.df_0, # dataframe to be pivoted
                                                  cols = C1:T3, # column names to be stored as a SINGLE variable
                                                  names_to = "samples", # name of that new variable (column)
                                                  values_to = "expression") # name of new variable (column) storing all the values (data)

log2.cpm.filtered.norm.df.pivot_10 <- pivot_longer(log2.cpm.filtered.norm.df_10, # dataframe to be pivoted
                                                  cols = C1:T3, # column names to be stored as a SINGLE variable
                                                  names_to = "samples", # name of that new variable (column)
                                                  values_to = "expression") # name of new variable (column) storing all the values (data)


p.norm.filtered_1 <- ggplot(log2.cpm.filtered.norm.df.pivot_1) +
  aes(x=samples, y=expression, fill=samples) +
  geom_violin(trim = FALSE, show.legend = FALSE) +
  stat_summary(fun = "median", 
               geom = "point", 
               shape = 95, 
               size = 10, 
               color = "black", 
               show.legend = FALSE) +
  labs(y="log2 expression", x = "sample",
       title="Log2 Counts per Million",
       subtitle="filtered cpm>1, TMM normalized",
       caption=paste0("produced on ", Sys.time())) +
  theme_bw()


p.norm.filtered_0 <- ggplot(log2.cpm.filtered.norm.df.pivot_0) +
  aes(x=samples, y=expression, fill=samples) +
  geom_violin(trim = FALSE, show.legend = FALSE) +
  stat_summary(fun = "median", 
               geom = "point", 
               shape = 95, 
               size = 10, 
               color = "black", 
               show.legend = FALSE) +
  labs(y="log2 expression", x = "sample",
       title="Log2 Counts per Million",
       subtitle="filtered cpm>0, TMM normalized",
       caption=paste0("produced on ", Sys.time())) +
  theme_bw()

p.norm.filtered_10 <- ggplot(log2.cpm.filtered.norm.df.pivot_10) +
  aes(x=samples, y=expression, fill=samples) +
  geom_violin(trim = FALSE, show.legend = FALSE) +
  stat_summary(fun = "median", 
               geom = "point", 
               shape = 95, 
               size = 10, 
               color = "black", 
               show.legend = FALSE) +
  labs(y="log2 expression", x = "sample",
       title="Log2 Counts per Million",
       subtitle="filtered cpm>10, TMM normalized",
       caption=paste0("produced on ", Sys.time())) +
  theme_bw()

beep(sound = 1)


plot_grid(p.raw.cpm, p.raw.log2cpm, p.log2cpmfilted_1,p.norm.filtered_1, labels = c('A', 'B', 'C', 'D'), label_size = 12)
plot_grid(p.raw.cpm, p.raw.log2cpm, p.log2cpmfilted_0,p.norm.filtered_0, labels = c('A', 'B', 'C', 'D'), label_size = 12)
plot_grid(p.raw.cpm, p.raw.log2cpm, p.log2cpmfilted_10,p.norm.filtered_10, labels = c('A', 'B', 'C', 'D'), label_size = 12)

