#Introduction----
#This script will pull the transcript IDs and genes from EnsDB of mouse 

#Libraries----
library(rhdf5) #provides functions for handling hdf5 file formats (kallisto outputs bootstraps in this format)
library(tidyverse) # provides access to Hadley Wickham's collection of R packages for data science, which we will use throughout the course
library(tximport) # package for getting Kallisto results into R
library(ensembldb) #helps deal with ensembl
library(EnsDb.Mmusculus.v79) #replace with your organism-specific database package
library(beepr) #just for fun
#library(biomaRt) #an alternative for annotation
#library(datapasta) #for copy/paste data into the R environment

#Read study design file----
#read_tsv function from readr package to read .txt study design file
targets <- read_tsv("metadata.txt")
#?read_tsv

#Create paths for abundance files to be read 
path <- file.path(targets$Sample_ID, "abundance.tsv") 
#check if paths exist
all(file.exists(path)) 
which(file.exists(path)) 
#any(file.exists(path)) 

#annotation using ensDB.MMusculus.v79 database
Tx <- transcripts(EnsDb.Mmusculus.v79, columns=c("tx_id", "gene_name")) %>%
      as_tibble() %>% #converts the large Granges Tx object to tibble 
      dplyr::rename(target_id = tx_id) %>%
      dplyr::select("target_id", "gene_name")

#TX <- transcripts(EnsDb.Mmusculus.v79, columns=c("tx_id", "gene_name"))
#TX <-  as_tibble(TX) #converts the large Granges Tx object to tibble 
#TX <-  dplyr::rename(TX, target_id = tx_id)
#TX <-  dplyr::select(TX, "target_id", "gene_name")

Txi_gene <- tximport(path, 
                     type = "kallisto", 
                     tx2gene = Tx, 
                     txOut = FALSE, #FALSE value gets gene level data
                     countsFromAbundance = "lengthScaledTPM",
                     ignoreTxVersion = TRUE)

beep(sound = 1)
summary(Txi_gene)
class(Txi_gene)
names(Txi_gene)
